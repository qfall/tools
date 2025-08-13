// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! A classical implementation of the FDH signature scheme using the [`PSFGPV`]
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
};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Initializes an FDH signature scheme from a [`PSFGPV`].
///
/// This function corresponds to an implementation of an FDH-signature
/// scheme with the explicit PSF [`PSFGPV`] which is generated using
/// the default of [`GadgetParameters`].
///
/// Attributes:
/// - `psf`: Defines the PSF needed for preimage sampling.
/// - `storage`: Stores all previously constructed signatures.
/// - `hash`: Defines the hash function going from Strings to the PSFs range.
///
/// Returns an explicit implementation of a FDH-signature scheme.
///
/// # Example
/// ```
/// use qfall_crypto::construction::signature::{fdh::FDHGPV, SignatureScheme};
///
/// let m = "Hello World!";
///
/// let mut fdh = FDHGPV::setup(4, 113, 17);
/// let (pk, sk) = fdh.gen();
///
/// let sigma = fdh.sign(m.to_string(), &sk, &pk);
///
/// assert!(fdh.vfy(m.to_string(), &sigma, &pk));
/// ```
#[derive(Serialize, Deserialize)]
pub struct FDHGPV {
    pub psf: PSFGPV,
    pub storage: HashMap<String, MatZ>,
    pub hash: HashMatZq,
}

impl FDHGPV {
    /// Initializes the [`FDHGPV`] with default parameters.
    /// The setup function takes in the security parameter, the modulus, a length bound
    /// for the signatures, and the length of randomness for this construction.
    /// Then, the [`PSFGPV`] is instantiated with the default [`GadgetParameters`].
    /// This PSF with an additional storage and hash function are secured in the struct.
    ///
    /// Parameters:
    /// - `n`: The security parameter
    /// - `q`: The modulus used for the G-Trapdoors
    /// - `s`: The Gaussian parameter with which is sampled
    ///
    /// Returns an explicit instantiation of a FDH signature scheme using the default
    /// parameters.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::signature::{fdh::FDHGPV, SignatureScheme};
    ///
    /// let mut fdh = FDHGPV::setup(4, 113, 17);
    /// ```
    ///
    /// # Panics ...
    /// - if the security parameter n is not in [1, i64::MAX].    
    /// - if `q <= 1`.
    pub fn setup(n: impl Into<Z>, q: impl Into<Modulus>, s: impl Into<Q>) -> Self {
        let (n, q, s) = (n.into(), q.into(), s.into());
        let psf = PSFGPV {
            gp: GadgetParameters::init_default(&n, &q),
            s,
        };
        Self {
            psf,
            storage: HashMap::new(),
            hash: HashMatZq {
                modulus: q,
                rows: i64::try_from(n).unwrap(),
                cols: 1,
            },
        }
    }
}

impl SignatureScheme for FDHGPV {
    /// The trapdoor and a precomputed short basis that speeds up preimage sampling.
    type SecretKey = (MatZ, MatQ);
    /// The public matrix defining the PSF, for which the secret key defines a trapdoor
    type PublicKey = MatZq;
    /// Defined by the domain of the PSF.
    type Signature = MatZ;

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
    use crate::construction::signature::{fdh::gpv::FDHGPV, SignatureScheme};
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

        let mut fdh = FDHGPV::setup(n, &q, &s);
        let (pk, sk) = fdh.gen();

        for i in 0..10 {
            let m = format!("Hello World! {i}");

            let sigma = fdh.sign(m.to_owned(), &sk, &pk);

            assert_eq!(&sigma, &fdh.sign(m.to_owned(), &sk, &pk));
            assert!(fdh.vfy(m.to_owned(), &sigma, &pk))
        }
    }

    /// Ensure that an entry is actually added to the local storage.
    #[test]
    fn storage_filled() {
        let mut fdh = FDHGPV::setup(5, 1024, 10);

        let m = "Hello World!";
        let (pk, sk) = fdh.gen();
        let _ = fdh.sign(m.to_owned(), &sk, &pk);

        assert!(fdh.storage.contains_key(m))
    }

    /// Ensure that after deserialization the HashMap still contains all entries.
    #[test]
    fn reload_hashmap() {
        let mut fdh = FDHGPV::setup(5, 1024, 10);

        // fill one entry in the HashMap
        let m = "Hello World!";
        let (pk, sk) = fdh.gen();
        let _ = fdh.sign(m.to_owned(), &sk, &pk);

        let fdh_string = serde_json::to_string(&fdh).expect("Unable to create a json object");
        let fdh_2: FDHGPV = serde_json::from_str(&fdh_string).unwrap();

        assert_eq!(fdh.storage, fdh_2.storage);
    }
}
