// Copyright © 2025 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implements a Perturbation MP12 PSF according to [\[1\]](<../index.html#:~:text=[1]>)
//! using G-Trapdoors and corresponding trapdoor.

use super::PSF;
use crate::sample::g_trapdoor::{
    gadget_classical::short_basis_gadget,
    gadget_classical::{find_solution_gadget_mat, gen_trapdoor},
    gadget_parameters::GadgetParameters,
};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    rational::{MatQ, Q},
    traits::{Concatenate, MatrixDimensions, Pow},
};
use serde::{Deserialize, Serialize};

/// A lattice-based implementation of a [`PSF`] according to
/// [\[1\]](<index.html#:~:text=[1]>) using
/// G-Trapdoors where D_n = {e ∈ Z^m | |e| <= s sqrt(m)}
/// and R_n = Z_q^n.
///
/// Attributes
/// - `gp`: Describes the gadget parameters with which the G-Trapdoor is generated
/// - `r`: The rounding parameter
/// - `s`: The Gaussian parameter with which is sampled
///
/// # Examples
/// ```
/// use qfall_crypto::primitive::psf::PSFPerturbation;
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
/// use qfall_math::rational::Q;
/// use qfall_crypto::primitive::psf::PSF;
///
/// let psf = PSFPerturbation {
///     gp: GadgetParameters::init_default(8, 64),
///     r: Q::from(3),
///     s: Q::from(25),
/// };
///
/// let (a, td) = psf.trap_gen();
/// let domain_sample = psf.samp_d();
/// let range_fa = psf.f_a(&a, &domain_sample);
/// let preimage = psf.samp_p(&a, &td, &range_fa);
///
/// assert!(psf.check_domain(&preimage));
/// assert_eq!(a * preimage, range_fa);
/// ```
#[derive(Serialize, Deserialize)]
pub struct PSFPerturbation {
    pub gp: GadgetParameters,
    pub r: Q,
    pub s: Q,
}

impl PSFPerturbation {
    /// Computes √Σ_2 = √(r^2 * (b^2 + 1) * [R^t | I]^t * [R^t | I] - r^2 * I)
    /// to perform non-spherical Gaussian sampling according to Algorithm 1 in
    /// [\[3\]](<index.html#:~:text=[3]>).
    /// This matrix is the second part of the secret key and needs to be precomputed
    /// to execute [PSFPerturbation::samp_p].
    /// [PSFPerturbation::trap_gen] outputs this matrix for `s^2 * I`, i.e. for discrete
    /// Gaussian preimages centered around `0`. This function enables changing the
    /// covariance matrix to any covariance matrix s.t. Σ_2 is positive definite.
    ///
    /// Parameters:
    /// - `mat_r`: The trapdoor matrix `R`
    /// - `mat_sigma`: The covariance matrix `Σ` to sample [`Self::samp_p`] with
    ///
    /// Returns a [`MatQ`] containing √Σ_2 = √(r^2 * (b^2 + 1) * [R^t | I]^t * [R^t | I] - r^2 * I)
    /// if Σ_2 was positive definite.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSFPerturbation;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::{Q, MatQ};
    /// use qfall_crypto::primitive::psf::PSF;
    /// use qfall_math::traits::*;
    ///
    /// let psf = PSFPerturbation {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     r: Q::from(3),
    ///     s: Q::from(25),
    /// };
    ///
    /// let (a, td) = psf.trap_gen();
    ///
    /// let cov_mat = psf.s.pow(2).unwrap() * &psf.r * MatQ::identity(a.get_num_columns(), a.get_num_columns());
    /// let mat_sqrt_sigma_2 = psf.compute_sqrt_sigma_2(&td.0, &cov_mat);
    /// let new_td = (td.0, mat_sqrt_sigma_2, td.2);
    ///
    /// let domain_sample = psf.samp_d();
    /// let range_fa = psf.f_a(&a, &domain_sample);
    /// let preimage = psf.samp_p(&a, &new_td, &range_fa);
    ///
    /// assert!(psf.check_domain(&preimage));
    /// ```
    ///
    /// # Panics ...
    /// - if Σ_2 is not positive definite.
    pub fn compute_sqrt_sigma_2(&self, mat_r: &MatZ, mat_sigma: &MatQ) -> MatQ {
        // Normalization factor according to MP12, Section 2.3
        let normalization_factor = 1.0 / (2.0 * Q::PI);

        // full_td = [R^t | I]^t
        let full_td = mat_r
            .concat_vertical(&MatZ::identity(
                mat_sigma.get_num_rows() - mat_r.get_num_rows(),
                mat_r.get_num_columns(),
            ))
            .unwrap();

        // Assemble Σ_p = Σ - [R^t|I]^t * Σ_G * [R^t|I] with √Σ_G ≥ η (Λ^⊥(G))
        // and assumption Σ_G = (base^2 + 1) * I as ||S|| = √(b^2 + 1), Theorem 4.1, MP12
        let mat_sigma_p: MatQ =
            mat_sigma - (self.gp.base.pow(2).unwrap() + 1) * &full_td * full_td.transpose();

        // r^2 * Σ_p <=> r * √Σ_p
        // Compute Σ_2 according to the requirements of Algorithm 1 in [3]
        // Assume Σ_1 = r^2 * B_1 * B_1^t with B_1 = I as basis for ZZ^n
        // Then, Σ_2 = Σ - Σ_1, where Σ = r^2 * Σ_p
        let sigma_2: MatQ = normalization_factor
            * self.r.pow(2).unwrap()
            * (&mat_sigma_p
                - MatQ::identity(mat_sigma_p.get_num_rows(), mat_sigma_p.get_num_columns()));

        // Compute √Σ_2
        sigma_2.cholesky_decomposition_flint()
    }
}

/// Samples a preimage with respect to the gadget matrix `G` of `vec_u`.
///
/// Parameters:
/// - `psf`: The [`PSFPerturbation`] that sets the parameters for this function
/// - `vec_u`: The target vector
/// - `short_basis_gadget`: The short basis of the corresponding gadget-matrix - just required to speed up the algorithm
/// - `short_basis_gadget_gso`: The GSO of the short basis of the gadget-matrix - just required to speed up the algorithm
///
/// Returns a discrete Gaussian sampled preimage of `u` under `G` with Gaussian parameter `r * √(b^2 + 1)`.
///
/// # Examples
/// ```compile_fail
/// use qfall_crypto::primitive::psf::PSFPerturbation;
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
/// use qfall_crypto::sample::g_trapdoor::gadget_classical::short_basis_gadget;
/// use qfall_math::rational::Q;
/// use qfall_crypto::primitive::psf::PSF;
///
/// let psf = PSFPerturbation {
///     gp: GadgetParameters::init_default(8, 64),
///     r: Q::from(3),
///     s: Q::from(25),
/// };
/// 
/// let short_basis_g = short_basis_gadget(&psf.gp);
/// let short_basis_g_gso = MatQ::from(&short_basis_g).gso();
///
/// let target = MatZq::sample_uniform(8, 1, 64);
/// let preimage = randomized_nearest_plane_gadget(&osf, &target, &short_basis_g, &short_basis_g_gso);
/// ```
pub(crate) fn randomized_nearest_plane_gadget(
    psf: &PSFPerturbation,
    vec_u: &MatZq,
    short_basis_gadget: &MatZ,
    short_basis_gadget_gso: &MatQ,
) -> MatZ {
    // s = r * √(b^2 + 1) according to Algorithm 3 in [1]
    let s = &psf.r * (psf.gp.base.pow(2).unwrap() + Z::ONE).sqrt();

    // find solution s.t. G * long_solution = vec_u
    let long_solution = find_solution_gadget_mat(vec_u, &psf.gp.k, &psf.gp.base);

    let center = MatQ::from(&(-1 * &long_solution));

    // just as PSFGPV::samp_p
    long_solution
        + MatZ::sample_d_precomputed_gso(
            short_basis_gadget,
            short_basis_gadget_gso,
            &psf.gp.n,
            &center,
            &s,
        )
        .unwrap()
}

impl PSF for PSFPerturbation {
    type A = MatZq;
    type Trapdoor = (MatZ, MatQ, (MatZ, MatQ));
    type Domain = MatZ;
    type Range = MatZq;

    /// Computes a G-Trapdoor according to the [`GadgetParameters`] defined in the struct.
    /// It returns a matrix `A` together with the short trapdoor matrix `R` and a precomputed √Σ_2
    /// for covariance matrix `s^2 * I`, i.e. for preimage sampling centered around `0`.
    /// The last part of the trapdoor tuple contains a short basis for the gadget matrix `G` and its
    /// GSO, which removes the need to recompute the GSO for each iteration and speeds up preimage sampling
    /// drastically.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSFPerturbation;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFPerturbation {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     r: Q::from(3),
    ///     s: Q::from(25),
    /// };
    ///
    /// let (mat_a, (mat_r, mat_sqrt_sigma_2, (sh_b_g, sh_b_g_gso))) = psf.trap_gen();
    /// ```
    fn trap_gen(&self) -> (MatZq, (MatZ, MatQ, (MatZ, MatQ))) {
        let mat_a_bar = MatZq::sample_uniform(&self.gp.n, &self.gp.m_bar, &self.gp.q);
        let tag = MatZq::identity(&self.gp.n, &self.gp.n, &self.gp.q);

        let (mat_a, mat_r) = gen_trapdoor(&self.gp, &mat_a_bar, &tag).unwrap();

        let mat_sqrt_sigma_2 = self.compute_sqrt_sigma_2(
            &mat_r,
            &(&self.s.pow(2).unwrap()
                * MatQ::identity(mat_a.get_num_columns(), mat_a.get_num_columns())),
        );

        let short_basis_gadget = short_basis_gadget(&self.gp);
        let short_basis_gadget_gso = MatQ::from(&short_basis_gadget).gso();

        (
            mat_a,
            (
                mat_r,
                mat_sqrt_sigma_2,
                (short_basis_gadget, short_basis_gadget_gso),
            ),
        )
    }

    /// Samples in the domain using SampleD with the standard basis and center `0`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSFPerturbation;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFPerturbation {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     r: Q::from(3),
    ///     s: Q::from(25),
    /// };
    /// let (mat_a, td) = psf.trap_gen();
    ///
    /// let domain_sample = psf.samp_d();
    /// ```
    fn samp_d(&self) -> MatZ {
        let m = &self.gp.n * &self.gp.k + &self.gp.m_bar;
        MatZ::sample_d_common(&m, &self.gp.n, &self.s).unwrap()
    }

    /// Samples an `e` in the domain using SampleD that is generated
    /// from the G-Trapdoor from the conditioned
    /// discrete Gaussian with `f_a(a,e) = u` for a provided syndrome `u`.
    ///
    /// *Note*: the provided parameters `mat_a, mat_r, vec_u` must fit together,
    /// otherwise unexpected behavior such as panics may occur.
    ///
    /// Parameters:
    /// - `mat_a`: The parity-check matrix
    /// - `mat_r`: The short trapdoor matrix `R`
    /// - `mat_sqrt_sigma_2`: The precomputed √Σ_2
    /// - `vec_u`: The syndrome from the range
    ///
    /// Returns a sample `e` from the domain on the conditioned discrete
    /// Gaussian distribution `f_a(a,e) = u` with covariance matrix Σ depending on `mat_sqrt_sigma_2`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSFPerturbation;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFPerturbation {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     r: Q::from(3),
    ///     s: Q::from(25),
    /// };
    /// let (mat_a, td) = psf.trap_gen();
    /// let domain_sample = psf.samp_d();
    /// let range_fa = psf.f_a(&mat_a, &domain_sample);
    ///
    /// let preimage = psf.samp_p(&mat_a, &td, &range_fa);
    /// assert_eq!(range_fa, psf.f_a(&mat_a, &preimage))
    /// ```
    fn samp_p(
        &self,
        mat_a: &MatZq,
        (mat_r, mat_sqrt_sigma_2, (short_basis_gadget, short_basis_gadget_gso)): &(
            MatZ,
            MatQ,
            (MatZ, MatQ),
        ),
        vec_u: &MatZq,
    ) -> MatZ {
        // Sample perturbation p <- D_{ZZ^m, r * √Σ_p} - not correct for now. √Σ_p := as √Σ_2
        let vec_p =
            MatZ::sample_d_common_non_spherical(&self.gp.n, mat_sqrt_sigma_2, &self.r).unwrap();

        // v = u - A * p
        let vec_v = vec_u - mat_a * &vec_p;

        // z <- D_{Λ_v^⊥(G), r * √Σ_G}
        let vec_z = randomized_nearest_plane_gadget(
            self,
            &vec_v,
            short_basis_gadget,
            short_basis_gadget_gso,
        );

        let full_td = mat_r
            .concat_vertical(&MatZ::identity(
                mat_r.get_num_columns(),
                mat_r.get_num_columns(),
            ))
            .unwrap();

        vec_p + full_td * vec_z
    }

    /// Implements the efficiently computable function `f_a` which here corresponds to
    /// `mat_a * sigma`. The sigma must be from the domain, i.e. D_n.
    ///
    /// Parameters:
    /// - `mat_a`: The parity-check matrix of dimensions `n x m`
    /// - `sigma`: A column vector of length `m`
    ///
    /// Returns `mat_a * sigma`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSFPerturbation;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFPerturbation {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     r: Q::from(3),
    ///     s: Q::from(25),
    /// };
    /// let (mat_a, td) = psf.trap_gen();
    /// let domain_sample = psf.samp_d();
    /// let range_fa = psf.f_a(&mat_a, &domain_sample);
    /// ```
    ///
    /// # Panics ...
    /// - if `sigma` is not in the domain.
    fn f_a(&self, mat_a: &MatZq, sigma: &MatZ) -> MatZq {
        assert!(self.check_domain(sigma));
        mat_a * sigma
    }

    /// Checks whether a value `sigma` is in D_n = {e ∈ Z^m | |e| <= s * r * sqrt(m)}.
    ///
    /// Parameters:
    /// - `sigma`: The value for which is checked, if it is in the domain
    ///
    /// Returns true, if `sigma` is in D_n.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSF;
    /// use qfall_crypto::primitive::psf::PSFPerturbation;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    ///
    /// let psf = PSFPerturbation {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     r: Q::from(3),
    ///     s: Q::from(25),
    /// };
    /// let (mat_a, td) = psf.trap_gen();
    ///
    /// let vector = psf.samp_d();
    ///
    /// assert!(psf.check_domain(&vector));
    /// ```
    fn check_domain(&self, sigma: &MatZ) -> bool {
        let m = &self.gp.n * &self.gp.k + &self.gp.m_bar;
        sigma.is_column_vector()
            && m == sigma.get_num_rows()
            && sigma.norm_eucl_sqrd().unwrap()
                <= self.s.pow(2).unwrap() * &m * &self.r.pow(2).unwrap()
    }
}

#[cfg(test)]
mod test_psf_perturbation {
    use super::PSFPerturbation;
    use super::PSF;
    use crate::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    use qfall_math::integer::MatZ;
    use qfall_math::rational::Q;
    use qfall_math::traits::*;

    /// Ensures that `samp_d` actually computes values that are in D_n.
    #[test]
    fn samp_d_samples_from_dn() {
        for (n, q) in [(5, 256), (10, 128), (15, 157)] {
            let psf = PSFPerturbation {
                gp: GadgetParameters::init_default(n, q),
                r: Q::from(n).log(2).unwrap(),
                s: Q::from(25),
            };

            for _ in 0..5 {
                assert!(psf.check_domain(&psf.samp_d()));
            }
        }
    }

    /// Ensures that `samp_p` actually computes preimages that are also in the correct
    /// domain.
    #[test]
    fn samp_p_preimage_and_domain() {
        for (n, q) in [(5, 256), (6, 128)] {
            let psf = PSFPerturbation {
                gp: GadgetParameters::init_default(n, q),
                r: Q::from(n).log(2).unwrap(),
                s: Q::from(25),
            };
            let (a, r) = psf.trap_gen();
            let domain_sample = psf.samp_d();
            let range_fa = psf.f_a(&a, &domain_sample);

            let preimage = psf.samp_p(&a, &r, &range_fa);
            assert_eq!(range_fa, psf.f_a(&a, &preimage));
            assert!(psf.check_domain(&preimage));
        }
    }

    /// Ensures that `f_a` returns `a * sigma`.
    #[test]
    fn f_a_works_as_expected() {
        for (n, q) in [(5, 256), (6, 128)] {
            let psf = PSFPerturbation {
                gp: GadgetParameters::init_default(n, q),
                r: Q::from(n).log(2).unwrap(),
                s: Q::from(25),
            };
            let (a, _) = psf.trap_gen();
            let domain_sample = psf.samp_d();

            assert_eq!(&a * &domain_sample, psf.f_a(&a, &domain_sample));
        }
    }

    /// Ensures that `f_a` panics if a value is provided, that is not within the domain.
    /// Sigma is not a vector.
    #[test]
    #[should_panic]
    fn f_a_sigma_not_in_domain_matrix() {
        let psf = PSFPerturbation {
            gp: GadgetParameters::init_default(8, 128),
            r: Q::from(8).log(2).unwrap(),
            s: Q::from(25),
        };
        let (a, _) = psf.trap_gen();
        let not_in_domain = MatZ::new(a.get_num_columns(), 2);

        let _ = psf.f_a(&a, &not_in_domain);
    }

    /// Ensures that `f_a` panics if a value is provided, that is not within the domain.
    /// Sigma is not of the correct length.
    #[test]
    #[should_panic]
    fn f_a_sigma_not_in_domain_incorrect_length() {
        let psf = PSFPerturbation {
            gp: GadgetParameters::init_default(8, 128),
            r: Q::from(8).log(2).unwrap(),
            s: Q::from(25),
        };
        let (a, _) = psf.trap_gen();
        let not_in_domain = MatZ::new(a.get_num_columns() - 1, 1);

        let _ = psf.f_a(&a, &not_in_domain);
    }

    /// Ensures that `f_a` panics if a value is provided, that is not within the domain.
    /// Sigma is too long.
    #[test]
    #[should_panic]
    fn f_a_sigma_not_in_domain_too_long() {
        let psf = PSFPerturbation {
            gp: GadgetParameters::init_default(8, 128),
            r: Q::from(8).log(2).unwrap(),
            s: Q::from(25),
        };
        let (a, _) = psf.trap_gen();
        let not_in_domain =
            psf.s.round() * a.get_num_columns() * MatZ::identity(a.get_num_columns(), 1);

        let _ = psf.f_a(&a, &not_in_domain);
    }

    /// Ensures that `check_domain` works for vectors with the correct length.
    #[test]
    fn check_domain_as_expected() {
        let psf = PSFPerturbation {
            gp: GadgetParameters::init_default(8, 128),
            r: Q::from(8).log(2).unwrap(),
            s: Q::from(25),
        };
        let (a, _) = psf.trap_gen();
        let value = psf.s.round();
        let mut in_domain = MatZ::new(a.get_num_columns(), 1);
        for i in 0..in_domain.get_num_rows() {
            in_domain.set_entry(i, 0, &value).unwrap();
        }

        assert!(psf.check_domain(&MatZ::new(a.get_num_columns(), 1)));
        assert!(psf.check_domain(&in_domain));
    }

    /// Ensures that `check_domain` returns false for values that are not in the domain.
    #[test]
    fn check_domain_not_in_dn() {
        let psf = PSFPerturbation {
            gp: GadgetParameters::init_default(8, 128),
            r: Q::from(8).log(2).unwrap(),
            s: Q::from(25),
        };
        let (a, _) = psf.trap_gen();

        let matrix = MatZ::new(a.get_num_columns(), 2);
        let too_short = MatZ::new(a.get_num_columns() - 1, 1);
        let too_long = MatZ::new(a.get_num_columns() + 1, 1);
        let entry_too_large =
            psf.s.round() * a.get_num_columns() * MatZ::identity(a.get_num_columns(), 1);

        assert!(!psf.check_domain(&matrix));
        assert!(!psf.check_domain(&too_long));
        assert!(!psf.check_domain(&too_short));
        assert!(!psf.check_domain(&entry_too_large));
    }
}
