// Copyright © 2023 Phil Milewski, Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! A classical implementation of the PFDH signature scheme using the [`PSFGPV`]
//! according to [\[1\]](<../index.html#:~:text=[1]>).

use crate::{
    construction::{
        hash::{sha256::HashMatZq, HashInto},
        signature::SignatureScheme,
    },
    primitive::psf::{PSF, PSFGPV},
    sample::g_trapdoor::gadget_parameters::GadgetParameters,
};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Modulus},
    rational::{MatQ, Q},
    traits::Pow,
};

/// Initializes an PFDH signature scheme from a [`PSFGPV`].
///
/// This function corresponds to an implementation of an PFDH-signature
/// scheme with the explicit PSF [`PSFGPV`] which is generated using
/// the default of [`GadgetParameters`].
///
/// Parameters:
/// - `n`: The security parameter
/// - `q`: The modulus used for the G-Trapdoors
/// - `s`: The Gaussian parameter with which is sampled
/// - `randomness_length`: the number of bits used for the randomness
///
/// Returns an explicit implementation of a PFDH-signature scheme.
///
/// # Example
/// ```
/// use qfall_crypto::construction::signature::{PFDHGPV, SignatureScheme};
///
/// let mut pfdh = PFDHGPV::setup(4, 113, 17, 128);
///
/// let m = "Hello World!";
///
/// let (pk, sk) = pfdh.gen();
/// let sigma = pfdh.sign(m.to_owned(), &sk, &pk);
///
/// assert!(pfdh.vfy(m.to_owned(), &sigma, &pk));
/// ```
///
/// # Panics ...
/// - if `q <= 1`.
pub struct PFDHGPV {
    pub psf: PSFGPV,
    pub hash: HashMatZq,
    pub randomness_length: Z,
}

impl PFDHGPV {
    pub fn setup(
        n: impl Into<Z>,
        q: impl Into<Modulus>,
        s: impl Into<Q>,
        randomness_length: impl Into<Z>,
    ) -> Self {
        let (n, q, s, randomness_length) = (n.into(), q.into(), s.into(), randomness_length.into());
        let psf = PSFGPV {
            gp: GadgetParameters::init_default(&n, &q),
            s,
        };
        let n = i64::try_from(&n).unwrap();
        Self {
            psf,
            hash: HashMatZq {
                modulus: q,
                rows: n,
                cols: 1,
            },
            randomness_length,
        }
    }
}

impl SignatureScheme for PFDHGPV {
    type SecretKey = (MatZ, MatQ);
    type PublicKey = MatZq;
    type Signature = (MatZ, Z);

    /// Generates a trapdoor by calling the `trap_gen` of the psf
    fn gen(&mut self) -> (Self::PublicKey, Self::SecretKey) {
        self.psf.trap_gen()
    }

    /// Firstly generate randomness
    /// It hashes the message and randomness into the domain and then computes a signature using
    /// `samp_p` from the psf with the trapdoor.
    fn sign(&mut self, m: String, sk: &Self::SecretKey, pk: &Self::PublicKey) -> Self::Signature {
        let randomness =
            Z::sample_uniform(0, Z::from(2).pow(&self.randomness_length).unwrap()).unwrap();
        let u = (self.hash).hash(&format!("{m} {randomness} {}", &self.randomness_length));
        let signature_part1 = self.psf.samp_p(pk, sk, &u);

        (signature_part1, randomness)
    }

    /// Checks if a signature is firstly within D_n, and then checks if
    /// the signature is actually a valid preimage under `fa` of `hash(m||r)`.
    fn vfy(&self, m: String, sigma: &Self::Signature, pk: &Self::PublicKey) -> bool {
        if !self.psf.check_domain(&sigma.0) {
            return false;
        }

        let u = (self.hash).hash(&format!("{m} {} {}", sigma.1, &self.randomness_length));

        self.psf.f_a(pk, &sigma.0) == u
    }
}

#[cfg(test)]
mod test_pfdh {
    use crate::construction::signature::{pfdh::gpv::PFDHGPV, SignatureScheme};
    use qfall_math::{integer::Z, rational::Q, traits::Pow};

    /// Ensure that the generated signature is valid.
    #[test]
    fn ensure_valid_signature_is_generated() {
        let n = Z::from(4);
        let k = Z::from(6);
        // `s >= ||\tilde short_base|| * omega(sqrt{log m})`,
        // here `log(2*n*k) = omega(sqrt{log m}))` (Theorem 4.1 - GPV08)
        let s: Q = ((&n * &k).sqrt() + 1) * Q::from(2) * (Z::from(2) * &n * &k).log(2).unwrap();
        let q = Z::from(2).pow(&k).unwrap();

        let mut pfdh = PFDHGPV::setup(n, &q, &s, 128);
        let (pk, sk) = pfdh.gen();

        for i in 0..10 {
            let m = format!("Hello World! {i}");

            let sigma = pfdh.sign(m.to_owned(), &sk, &pk);

            assert!(pfdh.vfy(m.to_owned(), &sigma, &pk))
        }
    }
}
