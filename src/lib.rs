// Copyright © 2023 Niklas Siemer, Marvin Beckmann
//
// This file is part of qFALL-tools.
//
// qFALL-tools is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! `qFALL` is a prototyping library for lattice-based cryptography.
//! `qFALL-tools` collects common sub-modules and features used by lattice-based constructions
//! to simplify and accelerate the development of such.
//! Among these are:
//! - [Compression techniques](crate::compression),
//! - [Primitives such as Preimage Samplable Functions (PSF)](crate::primitive),
//! - [Sampling algorithm using trapdoors](crate::sample), and
//! - [common functions for efficient prototyping](crate::utils) such as
//!   - [common message encodings for encryption](crate::utils::common_encodings),
//!   - [quick instantiations of common moduli for rings](crate::utils::common_moduli), as well as
//!   - [rotation matrices](crate::utils::rotation_matrix).
//! 
//! The `qFALL` project contains two more crates called [`qFALL-math`](https://crates.io/crates/qfall-math)
//! and [`qFALL-schemes`](https://crates.io/crates/qfall-schemes) to support prototyping.
//! - Find further information on [our website](https://qfall.github.io/).
//! - We recommend [our tutorial](https://qfall.github.io/book) to start working with qFALL.
//! 
//! ## Quick Example
//! ```
//! use qfall_tools::utils::{common_moduli::new_anticyclic, common_encodings::encode_value_in_polynomialringzq};
//! use qfall_math::integer::Z;
//! 
//! // Create X^256 + 1 mod 3329
//! let poly_mod = new_anticyclic(256, 3329).unwrap();
//! // Generate integer from string
//! let message = Z::from_utf8("Hello!");
//! // Turn string into encoding q/2 and 0 for each 1 and 0 bit respectively
//! let mu_q_half = encode_value_in_polynomialringzq(message, 2, &poly_mod).unwrap();
//! ```

pub mod compression;
pub mod primitive;
pub mod sample;
pub mod utils;
