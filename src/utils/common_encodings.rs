// Copyright © 2025 Niklas Siemer
//
// This file is part of qFALL-tools.
//
// qFALL-tools is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains functions to encode and decode data with commonly used
//! encodings.

use qfall_math::{
    integer::{PolyOverZ, Z},
    integer_mod_q::{ModulusPolynomialRingZq, PolynomialRingZq, Zq},
    traits::{Distance, GetCoefficient, SetCoefficient},
};

/// Turns a [`Z`] instance into its bit representation, converts this bit representation
/// into a [`PolynomialRingZq`] with entries q/2 for any 1-bit and 0 as coefficient for any 0-bit.
pub fn encode_z_bitwise_in_polynomialringzq(
    modulus: &ModulusPolynomialRingZq,
    mu: &Z,
) -> PolynomialRingZq {
    let modulus: ModulusPolynomialRingZq = modulus.into();
    let q_half = modulus.get_q().div_floor(2);

    let bits = mu.to_bits();
    let mut mu_q_half = PolynomialRingZq::from((PolyOverZ::default(), &modulus));
    for (i, bit) in bits.iter().enumerate() {
        if *bit {
            mu_q_half.set_coeff(i, &q_half).unwrap();
        }
    }

    mu_q_half
}

/// Checks for each coefficient of `poly` whether its value is closer to q/2 or 0
/// and sets the corresponding bit in the returned [`Z`] value to 1 or 0 respectively.
pub fn decode_z_bitwise_from_polynomialringzq(modulus: impl Into<Z>, poly: &PolynomialRingZq) -> Z {
    let q_half = modulus.into().div_floor(2);

    // check for each coefficient whether it's closer to 0 or q/2
    // if closer to q/2 -> add 2^i to result
    let mut vec = vec![];
    for i in 0..=poly.get_degree() {
        let coeff: Zq = poly.get_coeff(i).unwrap();
        let coeff: Z = coeff.get_representative_least_absolute_residue();

        if coeff.distance(&q_half) < coeff.distance(Z::ZERO) {
            vec.push(true);
        } else {
            vec.push(false);
        }
    }

    Z::from_bits(&vec)
}
