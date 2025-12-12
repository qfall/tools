// Copyright © 2025 Niklas Siemer
//
// This file is part of qFALL-tools.
//
// qFALL-tools is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Contains message encoding and decoding functions commonly used in lattice-based
//! encryption schemes.

use qfall_math::{
    error::MathError,
    integer::{PolyOverZ, Z},
    integer_mod_q::{ModulusPolynomialRingZq, PolynomialRingZq},
    traits::{GetCoefficient, SetCoefficient},
};

/// Takes any non-negative integer value, represents it with respect to `base`
/// and generates a [`PolynomialRingZq`] containing the previous representation w.r.t. `base`
/// across its coefficients multiplied by `q/base`.
/// This function is commonly used in encryption algorithms of lattice-based PKE schemes
/// and described as `⌊q/base * μ⌋`, where `μ ∈ R_{base}^n`.
/// 
/// Parameters:
/// - `value`: the non-negative integer value to encode
/// - `base`: defines the encoded representation, usually chosen as 2 for binary representation
/// - `modulus`: specifies the modulus of the returned struct and `q`
/// 
/// Returns a [`PolynomialRingZq`] containing `⌊q/base * μ⌋` as described above or a [`MathError`]
/// if `base < 2`, `value < 0`, or `value` represented w.r.t. `base` requires more than `modulus.get_degree()` coefficients.
/// 
/// # Examples
/// ```
/// use qfall_tools::utils::common_encodings::encode_value_in_polynomialringzq;
/// use qfall_tools::utils::common_moduli::new_anticyclic;
/// 
/// let modulus = new_anticyclic(16, 257).unwrap();
/// let base = 2;
/// let value = u16::MAX;
/// 
/// let encoded = encode_value_in_polynomialringzq(value, base, &modulus).unwrap();
/// ```
/// 
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
///   if `value` is negative, `value` represented w.r.t. `base` has more digits than coefficients available in `modulus.get_degree()`,
///   or `base < 2`.
pub fn encode_value_in_polynomialringzq(
    value: impl Into<Z>,
    base: impl Into<Z>,
    modulus: &ModulusPolynomialRingZq,
) -> Result<PolynomialRingZq, MathError> {
    let mut value = value.into();
    let base = base.into();
    let modulus: ModulusPolynomialRingZq = modulus.into();

    if value < Z::ZERO {
        return Err(MathError::InvalidIntegerInput(format!(
            "The given value {value} needs to be non-negative."
        )));
    }

    let min_req_degree = u64::try_from((&value + Z::ONE).log_ceil(&base)?)?;
    if min_req_degree > modulus.get_degree() as u64 {
        return Err(MathError::InvalidIntegerInput(format!(
            "The given value requires {min_req_degree} digits represented w.r.t. base {base}. Your modulus only provides space for {} digits.",
            modulus.get_degree()
        )));
    }

    // get representation of value w.r.t. base as base
    let mut base_repr = Vec::with_capacity(min_req_degree as usize);
    while value > Z::ZERO {
        let digit = &value % &base;
        base_repr.push(digit);
        value = value.div_floor(&base);
    }

    let mut res = PolyOverZ::default();
    for (i, digit) in base_repr.iter().enumerate() {
        if digit != &Z::ZERO {
            unsafe { res.set_coeff_unchecked(i as i64, digit) };
        }
    }

    // spread out each represented value by factor `q / base`
    let q_div_base = modulus.get_q().div_floor(&base);
    res *= q_div_base;

    Ok(PolynomialRingZq::from((res, modulus)))
}

/// Takes an encoded [`PolynomialRingZq`] and decodes it w.r.t. `base` and `q`, effectively performing
/// `μ = ⌈base/q * poly⌋ mod base` for `poly ∈ R_q^n`. Then, it takes any value of `μ ∈ R_{base}^n` and
/// turns it into a non-negative integer of type [`Z`]. 
/// This function is commonly used in decryption algorithms of lattice-based PKE schemes and invers to
/// [`encode_value_in_polynomialringzq`].
/// 
/// Parameters:
/// - `poly`: the [`PolynomialRingZq`] containing the encoded value
/// - `base`: defines the encoded representation, usually chosen as 2 for binary representation
/// 
/// Returns a [`Z`] containing the value of the vector `⌈base/q * poly⌋ mod base ∈ R_{base}^n` as a decimal number
/// as described above or a [`MathError`] if `base < 2`.
/// 
/// # Examples
/// ```
/// use qfall_tools::utils::common_encodings::{encode_value_in_polynomialringzq, decode_value_from_polynomialringzq};
/// use qfall_tools::utils::common_moduli::new_anticyclic;
/// 
/// let modulus = new_anticyclic(16, 257).unwrap();
/// let base = 2;
/// let value = u16::MAX;
/// 
/// let encoded = encode_value_in_polynomialringzq(value, base, &modulus).unwrap();
/// let decoded = decode_value_from_polynomialringzq(&encoded, base).unwrap();
/// 
/// assert_eq!(value, decoded);
/// ```
/// 
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
///   if `base < 2`.
pub fn decode_value_from_polynomialringzq(
    poly: &PolynomialRingZq,
    base: impl Into<Z>,
) -> Result<Z, MathError> {
    let base = base.into();
    let q = poly.get_mod().get_q();
    let q_div_2base = q.div_floor(2 * &base);

    if base <= Z::ONE {
        return Err(MathError::InvalidIntegerInput(format!(
            "The given base {base} is smaller than 2, which does not allow the encoding of any information."
        )));
    }

    let mut poly = poly.get_representative_least_nonnegative_residue();
    poly *= &base;

    let mut out = Z::default();

    for i in (0..=poly.get_degree()).rev() {
        let mut coeff = unsafe { poly.get_coeff_unchecked(i) };
        coeff += &q_div_2base;
        let res = coeff.div_floor(&q) % &base;
        out *= &base;
        out += res;
    }

    Ok(out)
}

#[cfg(test)]
mod test_encode_value_in_polynomialringzq {
    use crate::utils::{
        common_encodings::encode_value_in_polynomialringzq, common_moduli::new_anticyclic,
    };
    use qfall_math::{integer::Z, traits::GetCoefficient};

    /// Ensures that [`encode_value_in_polynomialringzq`] works properly for `base = 2`.
    #[test]
    fn binary() {
        let q = 257;
        let q_half = q / 2;
        let modulus = new_anticyclic(16, q).unwrap();

        let res0 = encode_value_in_polynomialringzq(1, 2, &modulus).unwrap();
        let res1 = encode_value_in_polynomialringzq(2, 2, &modulus).unwrap();
        let res2 = encode_value_in_polynomialringzq(3, 2, &modulus).unwrap();

        assert_eq!(GetCoefficient::<Z>::get_coeff(&res0, 0).unwrap(), q_half);
        assert_eq!(res0.get_degree(), 0);

        assert_eq!(GetCoefficient::<Z>::get_coeff(&res1, 0).unwrap(), 0);
        assert_eq!(GetCoefficient::<Z>::get_coeff(&res1, 1).unwrap(), q_half);
        assert_eq!(res1.get_degree(), 1);

        assert_eq!(GetCoefficient::<Z>::get_coeff(&res2, 0).unwrap(), q_half);
        assert_eq!(GetCoefficient::<Z>::get_coeff(&res2, 1).unwrap(), q_half);
        assert_eq!(res2.get_degree(), 1);
    }

    /// Ensures that [`encode_value_in_polynomialringzq`] works properly for `base = 3`.
    #[test]
    fn ternary() {
        let q = 257;
        let q_third = q / 3;
        let modulus = new_anticyclic(16, q).unwrap();

        let res0 = encode_value_in_polynomialringzq(1, 3, &modulus).unwrap();
        let res1 = encode_value_in_polynomialringzq(2, 3, &modulus).unwrap();
        let res2 = encode_value_in_polynomialringzq(3, 3, &modulus).unwrap();

        assert_eq!(GetCoefficient::<Z>::get_coeff(&res0, 0).unwrap(), q_third);
        assert_eq!(res0.get_degree(), 0);

        assert_eq!(
            GetCoefficient::<Z>::get_coeff(&res1, 0).unwrap(),
            2 * q_third
        );
        assert_eq!(res1.get_degree(), 0);

        assert_eq!(GetCoefficient::<Z>::get_coeff(&res2, 0).unwrap(), 0);
        assert_eq!(GetCoefficient::<Z>::get_coeff(&res2, 1).unwrap(), q_third);
        assert_eq!(res2.get_degree(), 1);
    }

    /// Ensures that [`encode_value_in_polynomialringzq`] returns an error if there is not enough space due to modulus size constraints.
    #[test]
    fn not_enough_space() {
        let modulus = new_anticyclic(16, 257).unwrap();

        let res = encode_value_in_polynomialringzq(u16::MAX as u32 + 1, 2, &modulus);

        assert!(res.is_err());
    }

    /// Ensures that [`encode_value_in_polynomialringzq`] returns an error if `base < 2`.
    #[test]
    fn too_small_base() {
        let modulus = new_anticyclic(16, 257).unwrap();

        let res = encode_value_in_polynomialringzq(4, 1, &modulus);

        assert!(res.is_err());
    }

    /// Ensures that [`encode_value_in_polynomialringzq`] returns an error if `value < 0`.
    #[test]
    fn neagive_value() {
        let modulus = new_anticyclic(16, 257).unwrap();

        let res = encode_value_in_polynomialringzq(-1, 1, &modulus);

        assert!(res.is_err());
    }
}

#[cfg(test)]
mod test_decode_value_from_polynomialringzq {
    use crate::utils::{
        common_encodings::{decode_value_from_polynomialringzq, encode_value_in_polynomialringzq},
        common_moduli::new_anticyclic,
    };
    use qfall_math::{integer::Z, integer_mod_q::PolynomialRingZq};

    /// Ensures that encoded information can be decoded without losing information for `base = 2`.
    #[test]
    fn round_trip_binary() {
        let q = 257;
        let base = 2;
        let modulus = new_anticyclic(17, q).unwrap();
        let msg = Z::sample_uniform(0, u16::MAX).unwrap();

        let encoding = encode_value_in_polynomialringzq(&msg, base, &modulus).unwrap();
        let decoding = decode_value_from_polynomialringzq(&encoding, base).unwrap();

        assert_eq!(msg, decoding);
    }

    /// Ensures that encoded information can be decoded without losing information for `base = 3`.
    #[test]
    fn round_trip_ternary() {
        let q = 257;
        let base = 3;
        let modulus = new_anticyclic(16, q).unwrap();
        let msg = Z::sample_uniform(0, u16::MAX).unwrap();

        let encoding = encode_value_in_polynomialringzq(&msg, base, &modulus).unwrap();
        let decoding = decode_value_from_polynomialringzq(&encoding, base).unwrap();

        assert_eq!(msg, decoding);
    }

    /// Ensures that [`decode_value_from_polynomialringzq`] returns an error if `base < 2`.
    #[test]
    fn too_small_base() {
        let modulus = new_anticyclic(16, 257).unwrap();
        let poly = PolynomialRingZq::sample_uniform(&modulus);

        let res = decode_value_from_polynomialringzq(&poly, 1);

        assert!(res.is_err());
    }
}
