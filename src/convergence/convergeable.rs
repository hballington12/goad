/// Trait for types that can be tracked for convergence.
/// Provides operations needed for online mean/variance computation.
pub trait Convergeable: Clone + Sized {
    /// Create a zero/empty version with the same structure
    fn zero_like(&self) -> Self;

    /// Weighted addition: combines self (with weight w1) and other (with weight w2).
    /// For simple quantities: (self * w1 + other * w2) / (w1 + w2)
    /// For derived quantities like asymmetry: uses appropriate weighting (e.g., by ScatCross)
    fn weighted_add(&self, other: &Self, self_weight: f32, other_weight: f32) -> Self;

    /// Element-wise multiplication (for computing x²)
    fn mul_elem(&self, other: &Self) -> Self;

    /// Element-wise division
    fn div_elem(&self, other: &Self) -> Self;

    /// Element-wise addition
    fn add_elem(&self, other: &Self) -> Self;

    /// Element-wise subtraction
    fn sub_elem(&self, other: &Self) -> Self;

    /// Scale by a scalar
    fn scale(&self, scalar: f32) -> Self;

    /// Element-wise square root (for SEM computation)
    fn sqrt_elem(&self) -> Self;

    /// Returns a pre-weighted version for convergence tracking.
    /// e.g., asymmetry becomes asymmetry * scat_cross
    fn to_weighted(&self) -> Self;

    /// Returns the weights for each field.
    /// e.g., asymmetry weight is scat_cross, powers weight is 1.0
    fn weights(&self) -> Self;
}

/// Tracks running statistics for convergence using Welford's online algorithm.
/// Matches the Python implementation in goad/convergence/convergable.py exactly.
/// Computes mean and standard error of the mean (SEM) incrementally.
#[derive(Debug)]
pub struct ConvergenceTracker<T: Convergeable> {
    i: usize, // iteration counter
    m: T,     // running weighted sum (stores value*weight accumulated via Welford)
    s: T,     // sum of squared deltas (for variance)
    w: T,     // running mean weight
}

impl<T: Convergeable> ConvergenceTracker<T> {
    /// Create a new tracker using a template for structure
    pub fn new(template: &T) -> Self {
        Self {
            i: 0,
            m: template.zero_like(),
            s: template.zero_like(),
            w: template.zero_like(),
        }
    }

    /// Update with a new result.
    /// Internally computes weighted value and weight via to_weighted()/weights(),
    /// then applies Welford's algorithm matching Python exactly.
    pub fn update(&mut self, result: &T) {
        self.i += 1;

        // Get pre-weighted value and weights from result
        let value = result.to_weighted();
        let weight = result.weights();

        if self.i == 1 {
            self.m = value;
            self.w = weight;
            // s stays zero
        } else {
            // delta = value - m_old
            let delta = value.sub_elem(&self.m);

            // m = m_old + delta / i
            self.m = self.m.add_elem(&delta.scale(1.0 / self.i as f32));

            // s = s + delta^2 * (i-1)/i
            let delta_sq = delta.mul_elem(&delta);
            let factor = (self.i - 1) as f32 / self.i as f32;
            self.s = self.s.add_elem(&delta_sq.scale(factor));

            // w = w_old + (weight - w_old) / i
            let dw = weight.sub_elem(&self.w);
            self.w = self.w.add_elem(&dw.scale(1.0 / self.i as f32));
        }
    }

    /// Get the current count
    pub fn count(&self) -> usize {
        self.i
    }

    /// Get the running mean: m / w
    pub fn mean(&self) -> T {
        if self.i == 0 {
            return self.m.zero_like();
        }
        self.m.div_elem(&self.w)
    }

    /// Get the standard error of the mean: sqrt(s / (i-1)^2 / w^2)
    pub fn sem(&self) -> T {
        if self.i < 2 {
            return self.m.zero_like();
        }
        // SEM = sqrt(s / (i-1)^2 / w^2)
        let n_minus_1 = (self.i - 1) as f32;
        let w_sq = self.w.mul_elem(&self.w);
        self.s
            .scale(1.0 / (n_minus_1 * n_minus_1))
            .div_elem(&w_sq)
            .sqrt_elem()
    }

    /// Merge another tracker into this one using Chan's parallel algorithm.
    ///
    /// This enables batched/parallel Welford computation. Workers can accumulate
    /// results locally, then the master merges the partial trackers.
    ///
    /// Mathematical basis (Chan et al.):
    ///   n = n_a + n_b
    ///   delta = mean_b - mean_a
    ///   mean = mean_a + delta * (n_b / n)
    ///   M2 = M2_a + M2_b + delta^2 * (n_a * n_b / n)
    pub fn merge(&mut self, other: &Self) {
        if other.i == 0 {
            return;
        }
        if self.i == 0 {
            self.i = other.i;
            self.m = other.m.clone();
            self.s = other.s.clone();
            self.w = other.w.clone();
            return;
        }

        let n_a = self.i as f32;
        let n_b = other.i as f32;
        let n = n_a + n_b;

        // delta = mean_b - mean_a (for weighted values)
        let delta = other.m.sub_elem(&self.m);

        // mean = mean_a + delta * (n_b / n)
        self.m = self.m.add_elem(&delta.scale(n_b / n));

        // M2 = M2_a + M2_b + delta^2 * (n_a * n_b / n)
        let delta_sq = delta.mul_elem(&delta);
        self.s = self
            .s
            .add_elem(&other.s)
            .add_elem(&delta_sq.scale(n_a * n_b / n));

        // Same formula for weights
        let dw = other.w.sub_elem(&self.w);
        self.w = self.w.add_elem(&dw.scale(n_b / n));

        self.i = (n_a + n_b) as usize;
    }
}
