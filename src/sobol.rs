//! Owen-scrambled Sobol sequence for orientation sampling.
//!
//! Three dimensions with full 32-bit direction numbers, so one seed gives
//! 2^32 distinct points. This replaces `sobol_burley`, whose 16-bit tables
//! limit a sequence to 2^16 points: past that it silently revisits the same
//! 65,536 points, which made long convergence runs re-weight one point set.
//!
//! Construction:
//! - Direction numbers from Joe & Kuo (new-joe-kuo-6.21201), dimensions 1-3.
//! - The point is the XOR of the direction numbers selected by the index bits
//!   (direct form, stateless, so any index can be drawn in any order).
//! - Each dimension is then Owen-scrambled with a seed-dependent hash-based
//!   nested uniform scramble (Laine-Karras form with the constants of Vegdahl's
//!   `sobol_burley`, MIT/Apache-2.0).
//!
//! The index itself is NOT shuffled. Shuffling all 32 index bits would give the
//! first 2^m points random high bits, and they would no longer form a net.
//! Without the shuffle every power-of-two prefix is a scrambled (t,m,3)-net,
//! and different seeds are independent randomisations of the same net.
//!
//! Indices past 2^32 move to a fresh scramble seed per 2^32 block, so the
//! sampler is unbounded.

/// Number of bits per direction number (and the per-seed sequence length 2^BITS).
const BITS: usize = 32;

/// Direction numbers v_k = m_k << (32 - k), k = 1..=32, for one dimension with
/// primitive polynomial degree `s`, interior coefficients `a` and initial `m`.
const fn directions(s: usize, a: u32, m_init: [u32; 2]) -> [u32; BITS] {
    let mut m = [0u32; BITS + 1]; // 1-indexed
    let mut k = 1;
    while k <= BITS {
        if s == 0 {
            m[k] = 1; // van der Corput
        } else if k <= s {
            m[k] = m_init[k - 1];
        } else {
            // m_k = 2 a_1 m_{k-1} ^ 4 a_2 m_{k-2} ^ ... ^ 2^s m_{k-s} ^ m_{k-s}
            let mut x = m[k - s] ^ (m[k - s] << s);
            let mut j = 1;
            while j < s {
                let a_j = (a >> (s - 1 - j)) & 1;
                if a_j == 1 {
                    x ^= m[k - j] << j;
                }
                j += 1;
            }
            m[k] = x;
        }
        k += 1;
    }
    let mut v = [0u32; BITS];
    let mut k = 1;
    while k <= BITS {
        v[k - 1] = m[k] << (BITS - k);
        k += 1;
    }
    v
}

/// Dimension 1 is van der Corput (m_k = 1); dimensions 2 and 3 from Joe & Kuo.
const DIRECTIONS: [[u32; BITS]; 3] = [
    directions(0, 0, [1, 0]),
    directions(1, 0, [1, 0]),
    directions(2, 1, [1, 3]),
];

/// Unscrambled Sobol point (as a 32-bit fixed-point fraction) in `dimension`.
#[inline]
pub fn sobol_u32(index: u32, dimension: usize) -> u32 {
    let v = &DIRECTIONS[dimension];
    let mut x = 0u32;
    let mut i = index;
    while i != 0 {
        x ^= v[i.trailing_zeros() as usize];
        i &= i - 1;
    }
    x
}

/// 32-bit integer hash (so zero does not map to zero).
#[inline]
fn hash(mut n: u32) -> u32 {
    n ^= 0xe6fe3beb;
    n ^= n >> 16;
    n = n.wrapping_mul(0x7feb352d);
    n ^= n >> 15;
    n = n.wrapping_mul(0x846ca68b);
    n ^= n >> 16;
    n
}

/// Hash-based nested uniform (Owen) scramble on a bit-reversed value.
/// Each output bit depends only on more significant input bits.
#[inline]
fn owen_scramble_rev(mut n_rev: u32, scramble: u32) -> u32 {
    n_rev ^= n_rev.wrapping_mul(0x3d20adea);
    n_rev = n_rev.wrapping_add(scramble);
    n_rev = n_rev.wrapping_mul((scramble >> 16) | 1);
    n_rev ^= n_rev.wrapping_mul(0x05526c56);
    n_rev ^= n_rev.wrapping_mul(0x53a22864);
    n_rev
}

/// Per-dimension scramble key for a seed.
#[inline]
fn scramble_key(seed: u32, dimension: usize) -> u32 {
    const DIM_KEYS: [u32; 3] = [0x912f69ba, 0x174f18ab, 0x691e72ca];
    hash(seed.wrapping_mul(0x9c8f2d3b) ^ DIM_KEYS[dimension])
}

/// Effective 32-bit seed for a 64-bit seed and 2^32-block number.
/// Block 0 of a seed below 2^32 uses the seed directly.
#[inline]
fn block_seed(seed: u64, block: u32) -> u32 {
    let lo = seed as u32;
    let hi = (seed >> 32) as u32;
    if hi == 0 && block == 0 {
        lo
    } else {
        hash(lo ^ hash(hi ^ hash(block)))
    }
}

/// Scrambled Sobol sample in [0, 1) for a 64-bit `index`, `dimension` in 0..3.
#[inline]
pub fn sample(index: u64, dimension: usize, seed: u64) -> f32 {
    let s = block_seed(seed, (index >> 32) as u32);
    let x = sobol_u32(index as u32, dimension);
    let y = owen_scramble_rev(x.reverse_bits(), scramble_key(s, dimension)).reverse_bits();
    // Top 24 bits -> exact f32 in [0, 1).
    (y >> 8) as f32 * (1.0 / (1u32 << 24) as f32)
}

/// Uniformly distributed orientation (alpha, beta, gamma) in degrees for `index`.
#[inline]
pub fn euler_angles(index: u64, seed: u64) -> (f32, f32, f32) {
    let u1 = sample(index, 0, seed);
    let u2 = sample(index, 1, seed);
    let u3 = sample(index, 2, seed);
    let alpha = u1 * 360.0;
    let beta = (1.0_f32 - u2 * 2.0).acos().to_degrees();
    let gamma = u3 * 360.0;
    (alpha, beta, gamma)
}
