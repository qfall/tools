// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! A ring implementation of the FDH signature scheme using the [`PSFGPVRing`]
//! according to [\[1\]](<../index.html#:~:text=[1]>).

use crate::{
    construction::{
        hash::{sha256::HashMatPolynomialRingZq, HashInto},
        signature::SignatureScheme,
    },
    primitive::psf::{PSFGPVRing, PSF},
    sample::g_trapdoor::gadget_parameters::GadgetParametersRing,
};
use qfall_math::{
    integer::{MatPolyOverZ, Z},
    integer_mod_q::{MatPolynomialRingZq, Modulus},
    rational::Q,
};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Initializes an FDH signature scheme from a [`PSFGPVRing`].
/// The trapdoor is sampled with a Gaussian parameter of 1.005
/// as done in [\[3\]](<index.html#:~:text=[3]>) who derived it from
/// [\[5\]](<index.html#:~:text=[5]>).
///
/// This function corresponds to an implementation of an FDH-signature
/// scheme with the explicit PSF [`PSFGPVRing`] which is generated using
/// the default of [`GadgetParametersRing`].
///
/// Parameters:
/// - `n`: The security parameter
/// - `q`: The modulus used for the G-Trapdoors
/// - `s`: The Gaussian parameter with which is sampled
///
/// Returns an explicit implementation of an FDH-signature scheme.
///
/// # TODO: Example
///
/// # Panics ...
/// - if `q <= 1`.
#[derive(Serialize, Deserialize)]
pub struct FDHGPVRing {
    pub psf: PSFGPVRing,
    pub storage: HashMap<String, MatPolyOverZ>,
    pub hash: HashMatPolynomialRingZq,
}

impl FDHGPVRing {
    pub fn setup(n: impl Into<Z>, q: impl Into<Modulus>, s: impl Into<Q>) -> Self {
        let (n, q, s) = (n.into(), q.into(), s.into());
        let psf = PSFGPVRing {
            gp: GadgetParametersRing::init_default(&n, &q),
            s,
            s_td: Q::from(1.005_f64),
        };
        let modulus = psf.gp.modulus.clone();
        Self {
            psf,
            storage: HashMap::new(),
            hash: HashMatPolynomialRingZq {
                modulus,
                rows: 1,
                cols: 1,
            },
        }
    }
}

impl SignatureScheme for FDHGPVRing {
    type SecretKey = (MatPolyOverZ, MatPolyOverZ);
    type PublicKey = MatPolynomialRingZq;
    type Signature = MatPolyOverZ;

    /// Generates a trapdoor by calling the `trap_gen` of the psf
    fn gen(&mut self) -> (Self::PublicKey, Self::SecretKey) {
        self.psf.trap_gen()
    }

    /// Firstly checks if the message has been signed before, and if, return that
    /// signature, else it continues.
    /// It hashes the message into the domain and then computes a signature using
    /// `samp_p` from the psf with the trapdoor.
    fn sign(&mut self, m: String, sk: &Self::SecretKey, pk: &Self::PublicKey) -> Self::Signature {
        // check if it is in the HashMap
        if let Some(sigma) = self.storage.get(&m) {
            return sigma.clone();
        }

        let u = (self.hash).hash(&m);
        let signature = self.psf.samp_p(pk, sk, &u);

        // insert signature in HashMap
        self.storage.insert(m, signature.clone());
        signature
    }

    /// Checks if a signature is firstly within D_n, and then checks if
    /// the signature is actually a valid preimage under `fa` of `hash(m)`.
    fn vfy(&self, m: String, sigma: &Self::Signature, pk: &Self::PublicKey) -> bool {
        if !self.psf.check_domain(sigma) {
            return false;
        }

        let u = (self.hash).hash(&m);

        self.psf.f_a(pk, sigma) == u
    }
}

#[cfg(test)]
mod test_fdh {
    use crate::construction::signature::{fdh::gpv_ring::FDHGPVRing, SignatureScheme};
    use qfall_math::rational::Q;

    const MODULUS: i64 = 512;
    const N: i64 = 8;
    fn compute_s() -> Q {
        ((2 * 2 * Q::from(1.005_f64) * Q::from(N).sqrt() + 1) * 2) * 4
    }

    /// Ensure that the generated signature is valid.
    #[test]
    fn ensure_valid_signature_is_generated() {
        let mut fdh = FDHGPVRing::setup(N, MODULUS, compute_s());
        let (pk, sk) = fdh.gen();

        for i in 0..10 {
            let m = &format!("Hello World! {i}");

            let sigma = fdh.sign(m.to_owned(), &sk, &pk);

            assert!(
                fdh.vfy(m.to_owned(), &sigma, &pk),
                "This is a probabilistic test and may fail with negligible probability. \
                As n is rather small here, try to rerun the test and check whether the \
                test fails again."
            )
        }
    }

    /// Ensure that an entry is actually added to the local storage.
    #[test]
    fn storage_filled() {
        let mut fdh = FDHGPVRing::setup(N, MODULUS, compute_s());

        let m = "Hello World!";
        let (pk, sk) = fdh.gen();
        let sign_1 = fdh.sign(m.to_owned(), &sk, &pk);
        let sign_2 = fdh.sign(m.to_owned(), &sk, &pk);

        assert!(fdh.storage.contains_key(m));
        assert_eq!(sign_1, sign_2);
    }

    /// Ensure that after deserialization the HashMap still contains all entries.
    #[test]
    fn reload_hashmap() {
        let mut fdh = FDHGPVRing::setup(N, MODULUS, compute_s());

        // fill one entry in the HashMap
        let m = "Hello World!";
        let (pk, sk) = fdh.gen();
        let _ = fdh.sign(m.to_owned(), &sk, &pk);

        let fdh_string = serde_json::to_string(&fdh).expect("Unable to create a json object");

        let fdh_2: FDHGPVRing = serde_json::from_str(&fdh_string).unwrap();

        assert_eq!(fdh.storage, fdh_2.storage);
    }
}
