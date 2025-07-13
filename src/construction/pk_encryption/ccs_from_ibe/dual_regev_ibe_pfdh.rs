// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! A classical implementation of the [`CCSfromIBE`] scheme using
//! the [`DualRegevIBE`] and [`PFDH`].

use super::CCSfromIBE;
use crate::construction::{identity_based_encryption::DualRegevIBE, signature::PFDHGPV};
use qfall_math::{integer::Z, integer_mod_q::Modulus, rational::Q};

impl CCSfromIBE<DualRegevIBE, PFDHGPV> {
    /// Initializes a [`CCSfromIBE`] PK encryption scheme from a [`DualRegevIBE`] and a [`PFDH`] signature.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter
    /// - `q`: specifies the modulus
    /// - `r`: specifies the Gaussian parameter used for the [`PSFGPV`] for the [`PFDH`]
    /// - `randomness_length`: specifies the number of bits added to the message before signing
    /// - `alpha`: specifies the Gaussian parameter used for encryption in
    ///   [`DualRegev`](crate::construction::pk_encryption::DualRegev) in the [`DualRegevIBE`]
    ///
    /// Returns an explicit implementation of an IND-CCA-secure public key encryption scheme.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::CCSfromIBE;
    ///
    /// let mut scheme = CCSfromIBE::init_dr_pfdh(4, 13933, 4, 10.77, 0.0021);
    /// ```
    ///
    /// # Panics ...
    /// - if `q <= 1`.
    /// - if `n < 1` or `n` does not fit into an [`i64`].
    pub fn init_dr_pfdh(
        n: impl Into<Z>, // security parameter
        q: impl Into<Modulus>,
        randomness_length: impl Into<Z>, // added to the message before signing
        r: impl Into<Q>,                 // Gaussian parameter for PSF
        alpha: impl Into<Q>,             // Gaussian parameter for Dual Regev Encryption
    ) -> Self {
        let n = n.into();
        let q = q.into();
        let r = r.into();

        let dr_ibe = DualRegevIBE::new(&n, &q, &r, alpha);
        let pfdh = PFDHGPV::setup(n, q, r, randomness_length);

        Self {
            ibe: dr_ibe,
            signature: pfdh,
        }
    }

    /// Initializes a [`CCSfromIBE`] PK encryption scheme from a [`DualRegevIBE`] and a [`PFDH`] signature
    /// from a given `n > 0`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter
    ///
    /// Returns an explicit implementation of an IND-CCA-secure public key encryption scheme
    /// chosen with appropriate parameters for given `n`..
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::CCSfromIBE;
    ///
    /// let mut scheme = CCSfromIBE::init_dr_pfdh_from_n(4);
    /// ```
    ///
    /// # Panics ...
    /// - if `n < 4` or `n` does not fit into an [`i64`].
    pub fn init_dr_pfdh_from_n(n: impl Into<Z>) -> Self {
        let n = n.into();
        assert!(
            n > 3,
            "n needs to be chosen larger than 3 for this function to work properly."
        );

        let ibe = DualRegevIBE::new_from_n(&n);
        let pfdh = PFDHGPV::setup(&n, &ibe.dual_regev.q, &ibe.psf.s, &n);

        Self {
            ibe,
            signature: pfdh,
        }
    }
}

#[cfg(test)]
mod test_ccs_from_ibe {
    use super::CCSfromIBE;
    use crate::construction::pk_encryption::PKEncryptionSchemeMut;
    use qfall_math::integer::Z;

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 0 and small n.
    #[test]
    fn cycle_zero() {
        let msg = Z::ZERO;
        let mut scheme = CCSfromIBE::init_dr_pfdh_from_n(4);

        let (pk, sk) = scheme.gen();
        let cipher = scheme.enc(&pk, &msg);
        let m = scheme.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 1 and small n.
    #[test]
    fn cycle_one() {
        let msg = Z::ONE;
        let mut scheme = CCSfromIBE::init_dr_pfdh_from_n(4);

        let (pk, sk) = scheme.gen();
        let cipher = scheme.enc(&pk, &msg);
        let m = scheme.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }
}
