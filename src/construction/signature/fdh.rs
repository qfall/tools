// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This Module contains a implementations of the full domain hash signature scheme,
//! which only has to be instantiated with a corresponding PSF, a storage and
//! a corresponding hash function.
//!
//! The constructions follow the general definition of a hash-then-sign signature scheme
//! that uses a hash function as in [\[1\]](<index.html#:~:text=[1]>) and a PSF.
//! This signature scheme uses a storage, so it is stateful.
//!
//! Requirements
//! - `psf`: The PSF which has to implement the [`PSF`](crate::primitive::psf::PSF) trait
//!   and must also be (de-)serializable.
//! - `storage`: A Hashmap that safes all previously signed messages and their signature
//! - `hash`: The hash-function which has to map a string into the correct domain
//!
//! # Example
//! ## Signature Scheme from [`PSFGPV`](crate::primitive::psf::PSFGPV)
//! ```
//! use qfall_crypto::construction::signature::{fdh::FDHGPV, SignatureScheme};
//!
//! let mut fdh = FDHGPV::setup(4, 113, 17);
//!
//! let m = "Hello World!";
//!
//! let (pk, sk) = fdh.gen();
//! let sigma = fdh.sign(m.to_owned(), &sk, &pk);
//!
//! assert!(fdh.vfy(m.to_owned(), &sigma, &pk));
//! ```

mod gpv;
mod gpv_ring;

pub use gpv::FDHGPV;
pub use gpv_ring::FDHGPVRing;
