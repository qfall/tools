// Copyright © 2023 Marcel Luca Schmidt, Marvin Beckmann
//
// This file is part of qFALL-tools.
//
// qFALL-tools is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Contains implementations of Preimage Samplable Functions, short [`PSF`].
//!
//! The main references are listed in the following
//! and will be further referenced in submodules by these numbers:
//! - \[1\] Micciancio, D., Peikert, C. (2012).
//!   Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller.
//!   In: Pointcheval, D., Johansson, T. (eds) Advances in Cryptology – EUROCRYPT 2012.
//!   EUROCRYPT 2012. Lecture Notes in Computer Science, vol 7237.
//!   Springer, Berlin, Heidelberg. <https://doi.org/10.1007/978-3-642-29011-4_41>
//! - \[2\] Gür, K.D., Polyakov, Y., Rohloff, K., Ryan, G.W. and Savas, E., 2018,
//!   January. Implementation and evaluation of improved Gaussian sampling for lattice
//!   trapdoors. In Proceedings of the 6th Workshop on Encrypted Computing & Applied
//!   Homomorphic Cryptography (pp. 61-71). <https://dl.acm.org/doi/pdf/10.1145/3267973.3267975>
//! - \[3\] Peikert, Chris.
//!   An efficient and parallel Gaussian sampler for lattices.
//!   In: Annual Cryptology Conference - CRYPTO 2010.
//!   Springer, Berlin, Heidelberg. <https://doi.org/10.1007/978-3-642-14623-7_5>

mod gpv;
mod gpv_ring;
mod mp_perturbation;

pub use gpv::PSFGPV;
pub use gpv_ring::PSFGPVRing;
pub use mp_perturbation::PSFPerturbation;

/// This trait should be implemented by all constructions that are
/// actual implementations of a preimage sampleable function.
/// A formal definition for these PSFs can be found in
/// [\[1\]](<index.html#:~:text=[1]>)
pub trait PSF {
    type A;
    type Trapdoor;
    type Domain;
    type Range;

    /// Samples a parity-check matrix and a trapdoor for that matrix.
    ///
    /// Returns the parity-check matrix and the trapdoor.
    fn trap_gen(&self) -> (Self::A, Self::Trapdoor);

    /// Samples an element in the domain according to a specified distribution.
    ///
    /// Returns the sampled element.
    fn samp_d(&self) -> Self::Domain;

    /// Samples an element `e` in the domain according to a specified distribution
    /// conditioned on `f_a(a, e) = u`.
    ///
    /// Parameters:
    /// - `a`: The parity-check matrix
    /// - `r`: The G-Trapdoor for `a`
    /// - `u`: The syndrome from the range
    ///
    /// Returns a sample `e` from the domain on the conditioned discrete
    /// Gaussian distribution `f_a(a,e) = u`.
    fn samp_p(&self, a: &Self::A, r: &Self::Trapdoor, u: &Self::Range) -> Self::Domain;

    /// Implements the efficiently computable function `f_a`,
    /// which is uniquely classified by `a`.
    ///
    /// Parameters:
    /// - `a`: The parity-check matrix of dimensions `n x m`
    /// - `sigma`: A column vector of length `m`
    ///
    /// Returns the result of `f_a`.
    fn f_a(&self, a: &Self::A, sigma: &Self::Domain) -> Self::Range;

    /// Checks whether an element is in the correct domain (and not just the correct type).
    ///
    /// Returns the result of the check as a boolean.
    fn check_domain(&self, sigma: &Self::Domain) -> bool;
}
