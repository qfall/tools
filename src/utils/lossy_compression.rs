// Copyright © 2025 Niklas Siemer
//
// This file is part of qFALL-tools.
//
// qFALL-tools is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implements lossy (de-)compression as specified in ML-KEM.
//!
//! Reference:
//! - \[1\] National Institute of Standards and Technology (2024).
//!   Module-Lattice-Based Key-Encapsulation Mechanism Standard.
//!   Federal Information Processing Standards Publication (FIPS 203).
//!   <https://doi.org/10.6028/NIST.FIPS.203>

use qfall_math::{
    integer::Z,
    integer_mod_q::{MatPolynomialRingZq, PolynomialRingZq},
    traits::{
        GetCoefficient, MatrixDimensions, MatrixGetEntry, MatrixSetEntry, Pow, SetCoefficient,
    },
};

/// This trait is implemented by data-structures, which may use lossy compression by dropping lower order bits
/// as specified in [\[1\]](<index.html#:~:text=[1]>).
pub trait LossyCompression {
    /// Compresses by keeping only `d` higher-order bits.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Compress_d(x) := ⌈(2^d / q) * x⌋ mod 2^d`.
    ///
    /// Parameters:
    /// - `d`: specifies the number of bits that is kept to represent values
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn compress(&mut self, d: impl Into<Z>);

    /// Decompresses a previously compressed value by mapping it to the closest recoverable value over the ring `Z_q`.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Decompress_d(y) := ⌈(q / 2^d) * y⌋`.
    ///
    /// Parameters:
    /// - `d`: specifies the number of bits that was kept during compression
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn decompress(&mut self, d: impl Into<Z>);
}

impl LossyCompression for PolynomialRingZq {
    /// Compresses by keeping only `d` higher-order bits.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Compress_d(x) := ⌈(2^d / q) * x⌋ mod 2^d`.
    ///
    /// Parameters:
    /// - `d`: specifies the number of bits that is kept to represent each value
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer_mod_q::PolynomialRingZq;
    /// use qfall_tools::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompression};
    ///
    /// let modulus = new_anticyclic(16, 257).unwrap();
    /// let mut poly = PolynomialRingZq::sample_uniform(&modulus);
    ///
    /// poly.compress(4);
    /// ```
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn compress(&mut self, d: impl Into<Z>) {
        let d = d.into();
        assert!(d >= Z::ONE, "Performing this function with d < 1 implies reducing mod 1, leaving no information to recover. Choose a larger parameter d.");
        let two_pow_d = Z::from(2).pow(d).unwrap();
        let q = self.get_mod().get_q();

        for coeff_i in 0..=self.get_degree() {
            let mut coeff: Z = unsafe { self.get_coeff_unchecked(coeff_i) };

            coeff *= &two_pow_d;
            let res = (coeff / &q).round() % &q;

            unsafe { self.set_coeff_unchecked(coeff_i, res) };
        }
    }

    /// Decompresses a previously compressed value by mapping it to the closest recoverable value over the ring `Z_q`.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Decompress_d(y) := ⌈(q / 2^d) * y⌋`.
    ///
    /// Parameters:
    /// - `d`: specifies the number of bits that was kept during compression
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer_mod_q::PolynomialRingZq;
    /// use qfall_tools::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompression};
    ///
    /// let modulus = new_anticyclic(16, 257).unwrap();
    /// let mut poly = PolynomialRingZq::sample_uniform(&modulus);
    ///
    /// poly.compress(4);
    /// poly.decompress(4);
    /// ```
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn decompress(&mut self, d: impl Into<Z>) {
        let d = d.into();
        assert!(d >= Z::ONE, "Performing this function with d < 1 implies reducing mod 1, leaving no information to recover. Choose a larger parameter d.");
        let two_pow_d = Z::from(2).pow(d).unwrap();
        let q = self.get_mod().get_q();

        for coeff_i in 0..=self.get_degree() {
            let mut coeff: Z = unsafe { self.get_coeff_unchecked(coeff_i) };

            coeff *= &q;
            let res = (coeff / &two_pow_d).round();

            unsafe { self.set_coeff_unchecked(coeff_i, res) };
        }
    }
}

impl LossyCompression for MatPolynomialRingZq {
    /// Compresses by keeping only `d` higher-order bits.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Compress_d(x) := ⌈(2^d / q) * x⌋ mod 2^d`.
    ///
    /// Parameters:
    /// - `d`: specifies the number of bits that is kept to represent each value
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer_mod_q::MatPolynomialRingZq;
    /// use qfall_tools::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompression};
    ///
    /// let modulus = new_anticyclic(16, 257).unwrap();
    /// let mut poly = MatPolynomialRingZq::sample_uniform(2, 3, &modulus);
    ///
    /// poly.compress(4);
    /// ```
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn compress(&mut self, d: impl Into<Z>) {
        let d = d.into();

        for row in 0..self.get_num_rows() {
            for col in 0..self.get_num_columns() {
                let mut entry: PolynomialRingZq = unsafe { self.get_entry_unchecked(row, col) };
                entry.compress(&d);
                unsafe { self.set_entry_unchecked(row, col, entry) };
            }
        }
    }

    /// Decompresses a previously compressed value by mapping it to the closest recoverable value over the ring `Z_q`.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Decompress_d(y) := ⌈(q / 2^d) * y⌋`.
    ///
    /// Parameters:
    /// - `d`: specifies the number of bits that was kept during compression
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer_mod_q::MatPolynomialRingZq;
    /// use qfall_tools::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompression};
    ///
    /// let modulus = new_anticyclic(16, 257).unwrap();
    /// let mut poly = MatPolynomialRingZq::sample_uniform(2, 3, &modulus);
    ///
    /// poly.compress(4);
    /// poly.decompress(4);
    /// ```
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn decompress(&mut self, d: impl Into<Z>) {
        let d = d.into();

        for row in 0..self.get_num_rows() {
            for col in 0..self.get_num_columns() {
                let mut entry: PolynomialRingZq = unsafe { self.get_entry_unchecked(row, col) };
                entry.decompress(&d);
                unsafe { self.set_entry_unchecked(row, col, entry) };
            }
        }
    }
}

#[cfg(test)]
mod test_compression_poly {
    use crate::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompression};
    use qfall_math::{
        integer::Z,
        integer_mod_q::PolynomialRingZq,
        traits::{Distance, GetCoefficient, Pow},
    };

    /// Ensures that decompressing compressed values results in values close to the original for small d.
    #[test]
    fn round_trip_small_d() {
        let d = 4;
        let q = Z::from(257);
        let modulus = new_anticyclic(16, &q).unwrap();
        let mut poly = PolynomialRingZq::sample_uniform(&modulus);
        let cmp = poly.clone();

        poly.compress(d);
        poly.decompress(d);

        for i in 0..modulus.get_degree() {
            let orig_coeff: Z = cmp.get_coeff(i).unwrap();
            let coeff: Z = poly.get_coeff(i).unwrap();

            let mut distance = orig_coeff.distance(coeff);
            if distance > &q / 2 {
                distance = &q - distance;
            }

            assert!(distance <= Z::from(2).pow(q.log_ceil(2).unwrap() - d - 1).unwrap());
        }
    }

    /// Ensures that decompressing compressed values results in values close to the original.
    #[test]
    fn round_trip() {
        let d = 11;
        let q = Z::from(3329);
        let modulus = new_anticyclic(16, &q).unwrap();
        let mut poly = PolynomialRingZq::sample_uniform(&modulus);
        let cmp = poly.clone();

        poly.compress(d);
        poly.decompress(d);

        for i in 0..modulus.get_degree() {
            let orig_coeff: Z = cmp.get_coeff(i).unwrap();
            let coeff: Z = poly.get_coeff(i).unwrap();

            let mut distance = orig_coeff.distance(coeff);
            if distance > &q / 2 {
                distance = &q - distance;
            }

            assert!(distance <= Z::from(2).pow(q.log_ceil(2).unwrap() - d - 1).unwrap());
        }
    }

    /// Ensures that the function panics if `d = 0` or smaller.
    #[test]
    #[should_panic]
    fn too_small_d() {
        let d = 0;
        let q = Z::from(3329);
        let modulus = new_anticyclic(16, &q).unwrap();
        let mut poly = PolynomialRingZq::sample_uniform(&modulus);

        poly.compress(d);
    }
}

#[cfg(test)]
mod test_compression_matrix {
    use crate::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompression};
    use qfall_math::{
        integer::Z,
        integer_mod_q::{MatPolynomialRingZq, PolynomialRingZq},
        traits::{Distance, GetCoefficient, MatrixDimensions, MatrixGetEntry, Pow},
    };

    /// Ensures that decompressing compressed values results in values close to the original.
    #[test]
    fn round_trip() {
        let d = 11;
        let q = Z::from(3329);
        let modulus = new_anticyclic(16, &q).unwrap();
        let mut matrix = MatPolynomialRingZq::sample_uniform(2, 2, &modulus);
        let cmp = matrix.clone();

        matrix.compress(d);
        matrix.decompress(d);

        for row in 0..matrix.get_num_rows() {
            for col in 0..matrix.get_num_columns() {
                let orig_entry: PolynomialRingZq = cmp.get_entry(row, col).unwrap();
                let entry: PolynomialRingZq = matrix.get_entry(row, col).unwrap();

                for i in 0..modulus.get_degree() {
                    let orig_coeff: Z = orig_entry.get_coeff(i).unwrap();
                    let coeff: Z = entry.get_coeff(i).unwrap();

                    let mut distance = orig_coeff.distance(coeff);
                    if distance > &q / 2 {
                        distance = &q - distance;
                    }

                    assert!(distance <= Z::from(2).pow(q.log_ceil(2).unwrap() - d - 1).unwrap());
                }
            }
        }
    }

    /// Ensures that the function panics if `d = 0` or smaller.
    #[test]
    #[should_panic]
    fn too_small_d() {
        let d = 0;
        let q = Z::from(3329);
        let modulus = new_anticyclic(16, &q).unwrap();
        let mut poly = MatPolynomialRingZq::sample_uniform(2, 3, &modulus);

        poly.compress(d);
    }
}
