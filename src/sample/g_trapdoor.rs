// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-tools.
//
// qFALL-tools is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! A G-Trapdoor is a form of a trapdoor for lattices
//! that allows for very efficient sampling.
//! This module contains implementations for G-Trapdoors in the classical and
//! in the ring setting.
//!
//! The main references are listed in the following
//! and will be further referenced in submodules by these numbers:
//! - \[1\] Micciancio, D., Peikert, C. (2012).
//!   Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller.
//!   In: Pointcheval, D., Johansson, T. (eds) Advances in Cryptology – EUROCRYPT 2012.
//!   EUROCRYPT 2012. Lecture Notes in Computer Science, vol 7237.
//!   Springer, Berlin, Heidelberg. <https://doi.org/10.1007/978-3-642-29011-4_41>
//! - \[2\] El Bansarkhani, R., Buchmann, J. (2014). Improvement and Efficient
//!   Implementation of a Lattice-Based Signature Scheme. In: Lange, T., Lauter, K.,
//!   Lisoněk, P. (eds) Selected Areas in Cryptography -- SAC 2013. SAC 2013. Lecture Notes
//!   in Computer Science(), vol 8282. Springer, Berlin, Heidelberg.
//!   <https://doi.org/10.1007/978-3-662-43414-7_3>
//! - \[3\] Gür, K.D., Polyakov, Y., Rohloff, K., Ryan, G.W. and Savas, E., 2018,
//!   January. Implementation and evaluation of improved Gaussian sampling for lattice
//!   trapdoors. In Proceedings of the 6th Workshop on Encrypted Computing & Applied
//!   Homomorphic Cryptography (pp. 61-71). <https://dl.acm.org/doi/pdf/10.1145/3267973.3267975>
//! - \[4\] Cash, D., Hofheinz, D., Kiltz, E., & Peikert, C. (2012).
//!   Bonsai trees, or how to delegate a lattice basis. Journal of cryptology, 25, 601-639.
//!   <https://doi.org/10.1007/s00145-011-9105-2>
//! - \[5\] Chen, Yuanmi, and Phong Q. Nguyen. "BKZ 2.0: Better lattice security
//!   estimates." International Conference on the Theory and Application of Cryptology and
//!   Information Security. Berlin, Heidelberg: Springer Berlin Heidelberg, 2011.

pub mod gadget_classical;
pub mod gadget_default;
pub mod gadget_parameters;
pub mod gadget_ring;
pub mod short_basis_classical;
pub mod short_basis_ring;
pub mod trapdoor_distribution;
