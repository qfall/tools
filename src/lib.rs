// Copyright © 2023 Niklas Siemer, Marvin Beckmann
//
// This file is part of qFALL-tools.
//
// qFALL-tools is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! # What is qFALL-tools?
//! qFall-tools provides cryptographic basics such as commonly used tools in lattice-based cryptography,
//! mathematical primitives, and samplable distributions to prototype
//! lattice-based cryptographic constructions and more.
//! Actual constructions can be found in [qfall-schemes](https://github.com/qfall/schemes)
//!
//! Our library provides further primitives useful for prototyping such as
//! [`PSFs`](primitive::psf::PSF) that can be used to implement constructions.
//!
//! qFALL-tools is free software: you can redistribute it and/or modify it under
//! the terms of the Mozilla Public License Version 2.0 as published by the
//! Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.
//!
//! ## Tutorial + Website
//! You can find a dedicated [tutorial](https://qfall.github.io/book/index.html) to qFALL-tools on our [website](https://qfall.github.io/).
//! The tutorial explains the basic steps starting from installation and
//! continues with basic usage.
//! qFALL-tools is co-developed together with qFALL-math and qFALL-schemes.

pub mod primitive;
pub mod sample;
pub mod utils;
