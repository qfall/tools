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

use flint_sys::fmpz_poly::fmpz_poly_set_coeff_fmpz;
use qfall_math::{
    integer::{MatPolyOverZ, PolyOverZ, Z},
    integer_mod_q::{MatPolynomialRingZq, ModulusPolynomialRingZq, PolynomialRingZq},
    traits::{GetCoefficient, MatrixDimensions, MatrixGetEntry, MatrixSetEntry, Pow},
};

/// This trait is implemented by data-structures, which may use lossy compression by dropping lower order bits
/// as specified in [\[1\]](<index.html#:~:text=[1]>).
pub trait LossyCompressionFIPS203 {
    /// Defines the datatype that the compressed value will have.
    type CompressedType;
    /// Defines the type of the modulus object.
    type ModulusType;

    /// Compresses by keeping only `d` higher-order bits.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Compress_d(x) := ⌈(2^d / q) * x⌋ mod 2^d`.
    ///
    /// Parameters:
    /// - `d`: specifies the number of bits that is kept to represent values
    ///
    /// Returns a new instance of type [`Self::CompressedType`] containing the compressed coefficients with a loss-factor
    /// defined by `q` and `d`.
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn lossy_compress(&self, d: impl Into<Z>) -> Self::CompressedType;

    /// Decompresses a previously compressed value by mapping it to the closest recoverable value over the ring `Z_q`.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Decompress_d(y) := ⌈(q / 2^d) * y⌋`.
    ///
    /// Parameters:
    /// - `d`: specifies the number of bits that was kept during compression
    ///
    /// Returns a new instance of type [`Self`] with decompressed values according to the loss-factor
    /// defined by `q` and `d`.
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn lossy_decompress(
        compressed: &Self::CompressedType,
        d: impl Into<Z>,
        modulus: &Self::ModulusType,
    ) -> Self;
}

impl LossyCompressionFIPS203 for PolynomialRingZq {
    type CompressedType = PolyOverZ;
    type ModulusType = ModulusPolynomialRingZq;

    /// Compresses by keeping only `d` higher-order bits.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Compress_d(x) := ⌈(2^d / q) * x⌋ mod 2^d`.
    ///
    /// Parameters:
    /// - `d`: specifies the number of bits that is kept to represent each value
    ///
    /// Returns a [`PolyOverZ`] containing the compressed coefficients with a loss-factor
    /// defined by `q` and `d`.
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer_mod_q::PolynomialRingZq;
    /// use qfall_tools::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompressionFIPS203};
    ///
    /// let modulus = new_anticyclic(16, 257).unwrap();
    /// let mut poly = PolynomialRingZq::sample_uniform(&modulus);
    ///
    /// let compressed = poly.lossy_compress(4);
    /// ```
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn lossy_compress(&self, d: impl Into<Z>) -> Self::CompressedType {
        let d = d.into();
        assert!(d >= Z::ONE, "Performing this function with d < 1 implies reducing mod 1, leaving no information to recover. Choose a larger parameter d.");
        let two_pow_d = Z::from(2).pow(d).unwrap();
        let q = self.get_mod().get_q();
        let q_div_2 = q.div_floor(2);

        let mut out = PolyOverZ::default();

        for coeff_i in 0..=self.get_degree() {
            let mut coeff: Z = unsafe { self.get_coeff_unchecked(coeff_i) };

            coeff *= &two_pow_d;
            coeff += &q_div_2;
            let mut res = coeff.div_floor(&q) % &q;

            unsafe {
                fmpz_poly_set_coeff_fmpz(out.get_fmpz_poly_struct(), coeff_i, res.get_fmpz());
            };
        }

        out
    }

    /// Decompresses a previously compressed value by mapping it to the closest recoverable value over the ring `Z_q`.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Decompress_d(y) := ⌈(q / 2^d) * y⌋`.
    ///
    /// Parameters:
    /// - `compressed`: specifies the compressed value
    /// - `d`: specifies the number of bits that was kept during compression
    /// - `modulus`: specifies the modulus of the returned value
    ///
    /// Returns a [`PolynomialRingZq`] with decompressed values according to the loss-factor
    /// defined by `q` and `d`.
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer_mod_q::PolynomialRingZq;
    /// use qfall_tools::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompressionFIPS203};
    ///
    /// let modulus = new_anticyclic(16, 257).unwrap();
    /// let mut poly = PolynomialRingZq::sample_uniform(&modulus);
    ///
    /// let compressed = poly.lossy_compress(4);
    /// let decompressed = PolynomialRingZq::lossy_decompress(&compressed, 4, &modulus);
    /// ```
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn lossy_decompress(
        compressed: &Self::CompressedType,
        d: impl Into<Z>,
        modulus: &Self::ModulusType,
    ) -> Self {
        let d = d.into();
        assert!(d >= Z::ONE, "Performing this function with d < 1 implies reducing mod 1, leaving no information to recover. Choose a larger parameter d.");
        let two_pow_d_minus_1 = Z::from(2).pow(d - 1).unwrap();
        let two_pow_d = &two_pow_d_minus_1 * 2;
        let q = modulus.get_q();

        let mut out = Self::from(modulus);

        for coeff_i in 0..=compressed.get_degree() {
            let mut coeff: Z = unsafe { compressed.get_coeff_unchecked(coeff_i) };

            coeff *= &q;
            coeff += &two_pow_d_minus_1;
            let mut res = coeff.div_floor(&two_pow_d);

            unsafe {
                fmpz_poly_set_coeff_fmpz(out.get_fmpz_poly_struct(), coeff_i, res.get_fmpz());
            };
        }

        out
    }
}

impl LossyCompressionFIPS203 for MatPolynomialRingZq {
    type CompressedType = MatPolyOverZ;
    type ModulusType = ModulusPolynomialRingZq;

    /// Compresses by keeping only `d` higher-order bits.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Compress_d(x) := ⌈(2^d / q) * x⌋ mod 2^d`.
    ///
    /// Parameters:
    /// - `d`: specifies the number of bits that is kept to represent each value
    ///
    /// Returns a [`MatPolyOverZ`] containing the compressed coefficients with a loss-factor
    /// defined by `q` and `d`.
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer_mod_q::MatPolynomialRingZq;
    /// use qfall_tools::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompressionFIPS203};
    ///
    /// let modulus = new_anticyclic(16, 257).unwrap();
    /// let mut poly = MatPolynomialRingZq::sample_uniform(2, 3, &modulus);
    ///
    /// let compressed = poly.lossy_compress(4);
    /// ```
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn lossy_compress(&self, d: impl Into<Z>) -> Self::CompressedType {
        let d = d.into();

        let mut out = MatPolyOverZ::new(self.get_num_rows(), self.get_num_columns());

        for row in 0..self.get_num_rows() {
            for col in 0..self.get_num_columns() {
                let entry: PolynomialRingZq = unsafe { self.get_entry_unchecked(row, col) };
                let res = entry.lossy_compress(&d);
                unsafe { out.set_entry_unchecked(row, col, res) };
            }
        }

        out
    }

    /// Decompresses a previously compressed value by mapping it to the closest recoverable value over the ring `Z_q`.
    /// This function modifies the value of `self` directly.
    ///
    /// The function is specified by `Decompress_d(y) := ⌈(q / 2^d) * y⌋`.
    ///
    /// Parameters:
    /// - `compressed`: specifies the compressed matrix
    /// - `d`: specifies the number of bits that was kept during compression
    /// - `modulus`: specifies the modulus of the returned value
    ///
    /// Returns a [`MatPolynomialRingZq`] with decompressed values according to the loss-factor
    /// defined by `q` and `d`.
    ///
    /// # Examples
    /// ```
    /// use qfall_math::integer_mod_q::MatPolynomialRingZq;
    /// use qfall_tools::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompressionFIPS203};
    ///
    /// let modulus = new_anticyclic(16, 257).unwrap();
    /// let mut poly = MatPolynomialRingZq::sample_uniform(2, 3, &modulus);
    ///
    /// let compressed = poly.lossy_compress(4);
    /// let decompressed = MatPolynomialRingZq::lossy_decompress(&compressed, 4, &modulus);
    /// ```
    ///
    /// # Panics ...
    /// - if `d` is smaller than `1`.
    fn lossy_decompress(
        compressed: &Self::CompressedType,
        d: impl Into<Z>,
        modulus: &Self::ModulusType,
    ) -> Self {
        let d = d.into();

        let mut out = MatPolynomialRingZq::new(
            compressed.get_num_rows(),
            compressed.get_num_columns(),
            modulus,
        );

        for row in 0..compressed.get_num_rows() {
            for col in 0..compressed.get_num_columns() {
                let entry: PolyOverZ = unsafe { compressed.get_entry_unchecked(row, col) };
                let res = PolynomialRingZq::lossy_decompress(&entry, &d, modulus);
                unsafe { out.set_entry_unchecked(row, col, res) };
            }
        }

        out
    }
}

#[cfg(test)]
mod test_compression_poly {
    use crate::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompressionFIPS203};
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
        let poly = PolynomialRingZq::sample_uniform(&modulus);

        let compressed = poly.lossy_compress(d);
        let decompressed = PolynomialRingZq::lossy_decompress(&compressed, d, &modulus);

        for i in 0..modulus.get_degree() {
            let orig_coeff: Z = poly.get_coeff(i).unwrap();
            let coeff: Z = decompressed.get_coeff(i).unwrap();

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
        let poly = PolynomialRingZq::sample_uniform(&modulus);

        let compressed = poly.lossy_compress(d);
        let decompressed = PolynomialRingZq::lossy_decompress(&compressed, d, &modulus);

        for i in 0..modulus.get_degree() {
            let orig_coeff: Z = poly.get_coeff(i).unwrap();
            let coeff: Z = decompressed.get_coeff(i).unwrap();

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
        let poly = PolynomialRingZq::sample_uniform(&modulus);

        poly.lossy_compress(d);
    }
}

#[cfg(test)]
mod test_compression_matrix {
    use crate::utils::{common_moduli::new_anticyclic, lossy_compression::LossyCompressionFIPS203};
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
        let matrix = MatPolynomialRingZq::sample_uniform(2, 2, &modulus);

        let compressed = matrix.lossy_compress(d);
        let decompressed = MatPolynomialRingZq::lossy_decompress(&compressed, d, &modulus);

        for row in 0..matrix.get_num_rows() {
            for col in 0..matrix.get_num_columns() {
                let orig_entry: PolynomialRingZq = matrix.get_entry(row, col).unwrap();
                let entry: PolynomialRingZq = decompressed.get_entry(row, col).unwrap();

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
        let poly = MatPolynomialRingZq::sample_uniform(2, 3, &modulus);

        poly.lossy_compress(d);
    }
}
