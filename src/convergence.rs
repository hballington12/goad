use std::collections::HashMap;
use std::fmt::Debug;

/// Tracks the convergence of values across multiple iterations of a Monte Carlo simulation.
/// Type T represents the type of data being tracked.
#[derive(Debug, Clone)]
pub struct Convergence<T> {
    /// Maps parameter names to their tracked statistics
    values: HashMap<String, ValueStats>,
    /// Threshold for determining convergence
    threshold: f32,
    /// Minimum number of iterations before checking convergence
    min_iterations: usize,
    /// Associated data (original parameters, context, etc.)
    data: Option<T>,
}

/// Statistics for a single value being tracked for convergence
#[derive(Debug, Clone)]
struct ValueStats {
    /// Number of iterations processed
    count: usize,
    /// Current mean value
    mean: f32,
    /// Sum of squared differences term (used for variance calculation)
    s: f32,
    /// History of raw values (optional)
    history: Vec<f32>,
}

impl<T: Debug + Clone> Convergence<T> {
    /// Create a new Convergence tracker with the specified threshold
    pub fn new(threshold: f32, min_iterations: usize) -> Self {
        Self {
            values: HashMap::new(),
            threshold,
            min_iterations,
            data: None,
        }
    }

    /// Create a new Convergence tracker with associated data
    pub fn with_data(threshold: f32, min_iterations: usize, data: T) -> Self {
        Self {
            values: HashMap::new(),
            threshold,
            min_iterations,
            data: Some(data),
        }
    }

    /// Set associated data
    pub fn set_data(&mut self, data: T) {
        self.data = Some(data);
    }

    /// Get a reference to the associated data
    pub fn data(&self) -> Option<&T> {
        self.data.as_ref()
    }

    /// Get a mutable reference to the associated data
    pub fn data_mut(&mut self) -> Option<&mut T> {
        self.data.as_mut()
    }

    /// Take ownership of the associated data, replacing it with None
    pub fn take_data(&mut self) -> Option<T> {
        self.data.take()
    }

    /// Add a new parameter to track
    pub fn add_parameter(&mut self, name: &str) {
        self.values.insert(
            name.to_string(),
            ValueStats {
                count: 0,
                mean: 0.0,
                s: 0.0,
                history: Vec::new(),
            },
        );
    }

    /// Check if a parameter exists
    pub fn has_parameter(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }

    /// Update a parameter with a new value from the current iteration
    pub fn update(&mut self, name: &str, value: f32) -> Result<(), String> {
        let stats = self
            .values
            .get_mut(name)
            .ok_or_else(|| format!("Parameter '{}' not found", name))?;

        if stats.count == 0 {
            // First iteration
            stats.mean = value;
            stats.count = 1;
            stats.s = 0.0;
        } else {
            // Subsequent iterations
            stats.count += 1;
            let delta = value - stats.mean;
            stats.mean += delta / stats.count as f32;
            stats.s += ((stats.count - 1) as f32 / stats.count as f32) * delta * delta;
        }

        stats.history.push(value);
        Ok(())
    }

    /// Get the current mean value for a parameter
    pub fn mean(&self, name: &str) -> Option<f32> {
        self.values.get(name).map(|stats| stats.mean)
    }

    /// Get the current variance estimate for a parameter
    pub fn variance(&self, name: &str) -> Option<f32> {
        self.values.get(name).and_then(|stats| {
            if stats.count > 1 {
                // This is s²
                Some(stats.s / (stats.count - 1) as f32)
            } else {
                None // Need at least 2 samples for variance
            }
        })
    }

    /// Get the current standard deviation estimate for a parameter
    pub fn std_dev(&self, name: &str) -> Option<f32> {
        self.variance(name).map(|var| var.sqrt())
    }

    /// Check if a specific parameter has converged
    pub fn has_parameter_converged(&self, name: &str) -> Option<bool> {
        let stats = self.values.get(name)?;

        // Don't check convergence until minimum iterations are reached
        if stats.count < self.min_iterations {
            return Some(false);
        }

        let std_dev = self.std_dev(name)?;
        let mean = self.mean(name)?;

        // Avoid division by zero
        if mean.abs() < 1e-10 {
            Some(std_dev < self.threshold)
        } else {
            // Relative error: standard deviation / mean < threshold
            Some((std_dev / mean.abs()) < self.threshold)
        }
    }

    /// Check if all parameters have converged
    pub fn has_converged(&self) -> bool {
        if self.values.is_empty() {
            return false;
        }

        self.values
            .keys()
            .all(|name| self.has_parameter_converged(name).unwrap_or(false))
    }

    /// Get the current iteration count
    pub fn iteration_count(&self) -> Option<usize> {
        self.values.values().next().map(|stats| stats.count)
    }

    /// Get all parameter names
    pub fn parameter_names(&self) -> Vec<String> {
        self.values.keys().cloned().collect()
    }

    /// Reset the convergence tracking
    pub fn reset(&mut self) {
        for stats in self.values.values_mut() {
            stats.count = 0;
            stats.mean = 0.0;
            stats.s = 0.0;
            stats.history.clear();
        }
    }

    /// Get the threshold value
    pub fn threshold(&self) -> f32 {
        self.threshold
    }

    /// Set the threshold value
    pub fn set_threshold(&mut self, threshold: f32) {
        self.threshold = threshold;
    }

    /// Get the minimum iterations
    pub fn min_iterations(&self) -> usize {
        self.min_iterations
    }

    /// Set the minimum iterations
    pub fn set_min_iterations(&mut self, min_iterations: usize) {
        self.min_iterations = min_iterations;
    }
}

/// Trait for types that can be tracked for convergence
pub trait Convergable {
    /// Register the parameters to be tracked in the convergence object
    fn register_parameters(&self, convergence: &mut Convergence<Self>)
    where
        Self: Sized;

    /// Update the convergence with values from this instance
    fn update_convergence(&self, convergence: &mut Convergence<Self>) -> Result<(), String>
    where
        Self: Sized;

    /// Create a new instance from the convergence mean values
    fn from_convergence_means(convergence: &Convergence<Self>) -> Self
    where
        Self: Sized;
}

/// Implementation of Convergable for Params
impl Convergable for crate::params::Params {
    fn register_parameters(&self, convergence: &mut Convergence<Self>) {
        if self.asymettry.is_some() {
            convergence.add_parameter("asymettry");
        }
        if self.scat_cross.is_some() {
            convergence.add_parameter("scat_cross");
        }
        if self.ext_cross.is_some() {
            convergence.add_parameter("ext_cross");
        }
        if self.albedo.is_some() {
            convergence.add_parameter("albedo");
        }
    }

    fn update_convergence(&self, convergence: &mut Convergence<Self>) -> Result<(), String> {
        if let Some(val) = self.asymettry {
            if convergence.has_parameter("asymettry") {
                convergence.update("asymettry", val)?;
            }
        }

        if let Some(val) = self.scat_cross {
            if convergence.has_parameter("scat_cross") {
                convergence.update("scat_cross", val)?;
            }
        }

        if let Some(val) = self.ext_cross {
            if convergence.has_parameter("ext_cross") {
                convergence.update("ext_cross", val)?;
            }
        }

        if let Some(val) = self.albedo {
            if convergence.has_parameter("albedo") {
                convergence.update("albedo", val)?;
            }
        }

        Ok(())
    }

    fn from_convergence_means(convergence: &Convergence<Self>) -> Self {
        Self {
            asymettry: convergence.mean("asymettry"),
            scat_cross: convergence.mean("scat_cross"),
            ext_cross: convergence.mean("ext_cross"),
            albedo: convergence.mean("albedo"),
        }
    }
}
