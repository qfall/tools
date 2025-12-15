// Copyright © 2025 Niklas Siemer
//
// This file is part of qFALL-tools.
//
// qFALL-tools is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Contains commonly used compression techniques in lattice-based cryptography.
//!
//! References:
//! - \[1\] National Institute of Standards and Technology (2024).
//!   Module-Lattice-Based Key-Encapsulation Mechanism Standard.
//!   Federal Information Processing Standards Publication (FIPS 203).
//!   <https://doi.org/10.6028/NIST.FIPS.203>

mod lossy_compression_fips203;

pub use lossy_compression_fips203::LossyCompressionFIPS203;
