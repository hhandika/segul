//! Deterministic pseudo-random number generator for tree search and bootstrapping.

/// Deterministic Xoshiro256** PRNG seeded via SplitMix64.
#[derive(Debug, Clone)]
pub(crate) struct SeededRng {
    state: [u64; 4],
}

impl SeededRng {
    /// Create a generator from a 64-bit seed.
    pub(crate) fn from_seed(seed: u64) -> Self {
        let mut sm = seed;
        let s0 = splitmix64(&mut sm);
        let s1 = splitmix64(&mut sm);
        let s2 = splitmix64(&mut sm);
        let s3 = splitmix64(&mut sm);
        let state = if (s0 | s1 | s2 | s3) == 0 {
            [
                0x9E37_79B9_7F4A_7C15,
                0xBF58_476D_1CE4_E5B9,
                0x94D0_49BB_1331_11EB,
                1,
            ]
        } else {
            [s0, s1, s2, s3]
        };
        Self { state }
    }

    /// Derive an independent stream from `base_seed` and `stream` index.
    pub(crate) fn derive(base_seed: u64, stream: u64) -> Self {
        let mut x = stream.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mixed = base_seed ^ splitmix64(&mut x);
        Self::from_seed(mixed)
    }

    /// Generate the next 64-bit unsigned integer.
    #[inline]
    pub(crate) fn next_u64(&mut self) -> u64 {
        let result = (self.state[1].wrapping_mul(5))
            .rotate_left(7)
            .wrapping_mul(9);
        let t = self.state[1] << 17;

        self.state[2] ^= self.state[0];
        self.state[3] ^= self.state[1];
        self.state[1] ^= self.state[2];
        self.state[0] ^= self.state[3];

        self.state[2] ^= t;
        self.state[3] = self.state[3].rotate_left(45);

        result
    }

    /// Generate a uniformly distributed integer in `0..upper`.
    ///
    /// # Panics
    /// Panics if `upper == 0`.
    #[inline]
    pub(crate) fn gen_range(&mut self, upper: usize) -> usize {
        assert!(upper > 0, "Upper bound must be positive");
        let u = upper as u64;
        let threshold = u.wrapping_neg() % u;
        loop {
            let r = self.next_u64();
            if r >= threshold {
                return (r % u) as usize;
            }
        }
    }

    /// Shuffle `slice` in place using Fisher-Yates.
    pub(crate) fn shuffle<T>(&mut self, slice: &mut [T]) {
        let len = slice.len();
        for i in (1..len).rev() {
            let j = self.gen_range(i + 1);
            slice.swap(i, j);
        }
    }
}

#[inline]
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rng_determinism_and_stream_separation() {
        let mut r1 = SeededRng::derive(42, 0);
        let mut r2 = SeededRng::derive(42, 0);
        let mut r3 = SeededRng::derive(42, 1);

        let a: Vec<u64> = (0..10).map(|_| r1.next_u64()).collect();
        let b: Vec<u64> = (0..10).map(|_| r2.next_u64()).collect();
        let c: Vec<u64> = (0..10).map(|_| r3.next_u64()).collect();

        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn test_shuffle_preserves_elements() {
        let mut rng = SeededRng::from_seed(123);
        let mut items: Vec<usize> = (0..12).collect();
        rng.shuffle(&mut items);
        assert_ne!(items, (0..12).collect::<Vec<_>>());
        items.sort_unstable();
        assert_eq!(items, (0..12).collect::<Vec<_>>());
    }
}
