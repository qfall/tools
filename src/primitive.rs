// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-tools.
//
// qFALL-tools is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Contains primitives that are useful for cryptographic
//! constructions, but are solely targeted to be used in other constructions.
//! 
//! More specifically, this module contains primitives that do not provide security
//! guarantees targeted at end-users such as confidentiality, integrity, etc.

pub mod psf;
