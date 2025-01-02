// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains functions to generate G-Trapdoors with default parameters.

use crate::sample::g_trapdoor::{
    gadget_classical::gen_trapdoor,
    gadget_parameters::{GadgetParameters, GadgetParametersRing},
    gadget_ring::gen_trapdoor_ring_lwe,
};
use qfall_math::{
    integer::{MatPolyOverZ, MatZ, PolyOverZ, Z},
    integer_mod_q::{MatPolynomialRingZq, MatZq, Modulus},
    rational::Q,
};

/// Computes a trapdoor with default values.
///
/// - `params` is computed using [`GadgetParameters::init_default`].
/// - `tag = I_n` is taken from [\[1\]](<index.html#:~:text=[1]>): Algorithm 1
///
/// Parameters:
/// - `n`: the security parameter
/// - `q`: the modulus for the trapdoor
///
/// Returns a matrix `a` and its gadget-trapdoor `r` as in [\[1\]](<index.html#:~:text=[1]>): Algorithm 1 for some fixed set of parameters [`GadgetParameters::init_default`].
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::gadget_default::gen_trapdoor_default;
///
/// let (a,r) = gen_trapdoor_default(42, 101);
/// ```
///
/// # Panics ...
/// - if the security parameter `n` is not in `[1, i64::MAX]`.
/// - if `q <= 1`.
pub fn gen_trapdoor_default(n: impl Into<Z>, q: impl Into<Modulus>) -> (MatZq, MatZ) {
    // panic if n < 1 (security parameter must be positive)
    let n = n.into();
    assert!(n >= Z::ONE);

    let params = GadgetParameters::init_default(n, q);

    // a_bar <-$ Z_q^{n * m_bar}
    let a_bar = MatZq::sample_uniform(&params.n, &params.m_bar, &params.q);

    // tag = I_n
    let tag = MatZq::identity(&params.n, &params.n, &params.q);

    // we can unwrap, as we compute the parameters on our own and
    // they should always work
    gen_trapdoor(&params, &a_bar, &tag).unwrap()
}

/// Computes a trapdoor with default values in a ring setting.
///
/// - `params` is computed using [`GadgetParametersRing::init_default`].
///
/// Parameters:
/// - `n`: the security parameter
/// - `q`: the modulus for the trapdoor
///
/// Returns a matrix `a` and its gadget-trapdoor `(r,e)` as in [\[2\]](<index.html#:~:text=[2]>):
/// Construction 1 for some fixed set of parameters [`GadgetParameters::init_default`].
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::gadget_default::gen_trapdoor_ring_default;
///
/// let (a,r, e) = gen_trapdoor_ring_default(100, 29, 10);;
/// ```
///
/// # Panics...
/// - if the security parameter `n` is not in `[1, i64::MAX]`.
/// - if `q <= 1`.
pub fn gen_trapdoor_ring_default(
    n: impl Into<Z>,
    q: impl Into<Modulus>,
    s: impl Into<Q>,
) -> (MatPolynomialRingZq, MatPolyOverZ, MatPolyOverZ) {
    // panic if n < 1 (security parameter must be positive)
    let n = n.into();
    assert!(n >= Z::ONE);
    let s = s.into();

    let params = GadgetParametersRing::init_default(n, q);

    // a_bar <-$ Zq[X]^n
    let a_bar = PolyOverZ::sample_uniform(&params.n, 0, params.modulus.get_q()).unwrap();

    // we can unwrap, as we compute the parameters on our own and
    // they should always work
    gen_trapdoor_ring_lwe(&params, &a_bar, &s).unwrap()
}

#[cfg(test)]
mod test_gen_trapdoor_default {
    use super::gen_trapdoor_default;
    use crate::sample::g_trapdoor::gadget_classical::gen_gadget_mat;
    use qfall_math::{
        integer::{MatZ, Z},
        traits::{Concatenate, GetNumColumns, GetNumRows, Pow},
    };

    /// Ensures that the default parameters are used correctly and the expected
    /// dimensions are returned.
    #[test]
    fn correct_default_dimensions() {
        for n in [5, 10, 50] {
            for k in [5, 10] {
                let q = 2_i64.pow(k);

                let n_log_2_pow_2 = Z::from(n).log_ceil(2).unwrap().pow(2).unwrap();
                let m_bar = n * k + n_log_2_pow_2;
                let m = &m_bar + n * k;

                let (a, r) = gen_trapdoor_default(n, q);

                assert_eq!(n as i64, a.get_num_rows());
                assert_eq!(m, Z::from(a.get_num_columns()));

                assert_eq!(m_bar, Z::from(r.get_num_rows()));
                assert_eq!((n * k) as i64, r.get_num_columns());
            }
        }
    }

    /// Ensures that for several parameter choices the generated G-Trapdoor is
    /// actually a trapdoor.
    #[test]
    fn ensure_is_trapdoor() {
        for n in [5, 10, 25] {
            for k in [5, 10] {
                let q = 2_i64.pow(k);

                let (a, r) = gen_trapdoor_default(n, q);

                let trapdoor = r.concat_vertical(&MatZ::identity(n * k, n * k)).unwrap();

                assert_eq!(
                    gen_gadget_mat(n, k, &Z::from(2)),
                    (a * trapdoor).get_representative_least_nonnegative_residue()
                )
            }
        }
    }
}
