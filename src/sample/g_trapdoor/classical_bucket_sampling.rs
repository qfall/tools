// Copyright © 2024 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation to sample from a gadget-lattice using a discrete gaussian distribution.

use super::gadget_parameters::GadgetParameters;
use qfall_math::{
    integer::{MatZ, Z},
    rational::Q,
    traits::{Concatenate, GetEntry, GetNumRows, SetEntry},
};

/// This function allows to sample from `Λ^⟂_z(G)` where `z` is a solution for `Gz = u`
/// with gaussian parameter `s` where G is the gadget matrix defined by the gadget
///  parameters `gp`. It implements the Gaussian sampling described in [\[1\]](<../index.html#:~:text=[1]>)
/// section 4.1 for an arbitrary base, but only for a perfect-base modulus.
///
/// Parameters:
/// - `gp`: The gadget parameter that define the gadget matrix
/// - `s`: The Gaussian parameter for the samples
/// - `u`: The syndrome for which we have to find the solution
///
/// Returns a preimage `x` satisfying `G*x = u` where `G` is the gadget-matrix
/// with Gaussian parameter `s`.
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
/// use qfall_math::rational::Q;
/// use qfall_crypto::primitive::psf::PSF;
/// use qfall_math::integer::MatZ;
/// use qfall_crypto::sample::g_trapdoor::classical_bucket_sampling::sampling_gaussian_gadget;
///
/// let gp = GadgetParameters::init_default(4, 64);
/// let s = Q::from(16);
/// let u = MatZ::new(4,1);
///
/// let z = sampling_gaussian_gadget(&gp, s, &u);
/// ```
///
/// # Panics...
/// - if `s <= 0`.
/// - if `base` does not fit into an [`i64`].
/// - if `u` does not have `n` rows (i.e. more or less than the size of the gadget matrix)
pub fn sampling_gaussian_gadget(gp: &GadgetParameters, s: impl Into<Q>, u: &MatZ) -> MatZ {
    assert_eq!(gp.n, u.get_num_rows().into(), "The size of the target vector and the size of the gadget matrix do not fit together. The target vector has {} rows and the gadget matrix has {} rows.", u.get_num_rows(), gp.n);

    let s = s.into();

    // a total of `base`-many buckets each containing Gaussian samples
    let mut buckets = vec![Vec::new(); i64::try_from(&gp.base).unwrap() as usize];

    // repeat the Gaussian sampling using the nearest-plane algorithm for each dimension
    let mut matz_collection = (0..u.get_num_rows()).map(|i| {
        let u_i = u.get_entry(i, 0).unwrap();
        bucket_sampling(&gp.n, &s, &gp.k, &gp.base, &u_i, &mut buckets)
    });

    // combine all samples to one single vector
    let first = matz_collection.next().unwrap();
    matz_collection.fold(first, |acc, next| acc.concat_vertical(&next).unwrap())
}

/// This is a helper function for [`sampling_gaussian_gadget`] that samples
/// a solution for `g^t * out = u` where `u` is a target value and `g` is a gadget vector.
///
/// Parameters:
/// - `n`: the security parameter
/// - `s`: the standard deviation with which the element is sampled
/// - `k`: the size of the vector
/// - `base`: the amount of buckets and base for the modulus
/// - `u`: the value which restraints the sample
/// - `buckets`: the buckets in which the samples are saved.
///
/// Returns a vector such that `g^t * out = u`
///
/// # Panics...
/// - if `s <= 0``.
/// - if `base` does not fit into an [`i64`].
fn bucket_sampling(
    n: impl Into<Z>,
    s: impl Into<Q>,
    k: impl Into<Z>,
    base: impl Into<Z>,
    u: &Z,
    buckets: &mut [Vec<Z>],
) -> MatZ {
    let n = n.into();
    let s = s.into();
    let base = base.into();

    let mut u = u.clone();

    let mut out = MatZ::new(k.into(), 1);

    for i in 0..out.get_num_rows() {
        // find out from which bucket to take a sample
        let remainder = u.modulo(&base);
        let remainder_usize = i64::try_from(remainder).unwrap() as usize;

        // check if the correct bucket has a valid sample and otherwise
        // compute new samples until the bucket contains a sample
        while buckets[remainder_usize].is_empty() {
            let sample = Z::sample_discrete_gauss(&n, &Z::ZERO, &s).unwrap();

            let sample_remainder = sample.modulo(&base);
            let sample_remainder_usize = i64::try_from(&sample_remainder).unwrap() as usize;

            buckets[sample_remainder_usize].push(sample);
        }

        // grab a sample from the respective bucket
        let sample = buckets[remainder_usize].pop().unwrap();
        out.set_entry(i, 0, &sample).unwrap();

        u = (u - &sample).div_exact(&base).unwrap();
    }
    out
}

#[cfg(test)]
mod test_sampling_gaussian_gadget {
    use super::sampling_gaussian_gadget;
    use crate::sample::g_trapdoor::{
        gadget_classical::gen_gadget_mat, gadget_parameters::GadgetParameters,
    };
    use qfall_math::{integer::MatZ, integer_mod_q::MatZq, rational::Q};
    use std::str::FromStr;

    /// This test ensures that the dimension of the returned values is correct and
    /// that they are sampled from the correct lattice-coset.
    #[test]
    fn ensure_correct_coset() {
        let targets = [
            MatZ::new(4, 1),
            MatZ::from_str("[[15],[13],[0],[1]]").unwrap(),
        ];
        for u in targets {
            let gp = GadgetParameters::init_default(4, 64);
            let s = Q::from(16);

            let z = sampling_gaussian_gadget(&gp, s, &u);
            println!("{z}");

            assert_eq!(
                MatZq::from((&(gen_gadget_mat(&gp.n, &gp.k, &gp.base) * z), &gp.q)),
                MatZq::from((&u, &gp.q))
            );
        }
    }

    /// This test ensures that the function panics if an `s` is provided that is too small.
    #[test]
    #[should_panic]
    fn too_small_gaussian_parameter() {
        let gp = GadgetParameters::init_default(4, 64);
        let s = Q::from(-1);
        let u = MatZ::new(4, 1);

        let _ = sampling_gaussian_gadget(&gp, s, &u);
    }

    /// This test ensures that the function panics if an `u` does not fit to the gadget
    /// matrix in terms of dimensions.
    #[test]
    #[should_panic]
    fn too_small_targetvector_dimension() {
        let gp = GadgetParameters::init_default(4, 64);
        let s = Q::from(16);
        let u = MatZ::new(3, 1);

        let _ = sampling_gaussian_gadget(&gp, s, &u);
    }

    /// This test ensures that the function panics if an `u` does not fit to the gadget
    /// matrix in terms of dimensions.
    #[test]
    #[should_panic]
    fn too_large_targetvector_dimension() {
        let gp = GadgetParameters::init_default(4, 64);
        let s = Q::from(16);
        let u = MatZ::new(5, 1);

        let _ = sampling_gaussian_gadget(&gp, s, &u);
    }
}
