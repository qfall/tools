// Copyright © 2024 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implements a GPV PSF according to [\[1\]](<../index.html#:~:text=[1]>)
//! using G-Trapdoors and efficient Gaussian sampling by utilizing perturbation sampling.

use super::PSF;
use crate::sample::g_trapdoor::{
    classical_bucket_sampling::sampling_gaussian_gadget, gadget_classical::gen_trapdoor,
    gadget_parameters::GadgetParameters,
};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    rational::{MatQ, Q},
    traits::{Concatenate, GetNumColumns, GetNumRows, Pow},
};
use serde::{Deserialize, Serialize};

/// A lattice-based implementation of a [`PSF`] according to
/// [\[1\]](<index.html#:~:text=[1]>) using
/// G-Trapdoors where D_n = {e ∈ Z^m | |e| <= s * rounding_parameter * sqrt(m)}
/// and R_n = Z_q^n.
///
/// *Note*: This implementation only works reliably for perfect base-moduli.
/// It is very similar to [`PSFGPV`](super::PSFGPV) with the exception of `samp_p`,
/// because `samp_p` uses perturbation sampling instead of computing a short basis for
/// `Λ^⟂(A)` and that the samples have Gaussian parameter `s * rounding_parameter`
/// instead of only `s`.
///
/// Attributes
/// - `gp`: Describes the gadget parameters with which the G-Trapdoor is generated
/// - `s`: The Gaussian parameter with which is sampled
/// - `rounding_parameter`: The rounding parameter used for the discrete Gaussian
/// - `sigma_`: TODO
///
/// # Examples
/// ```
/// use qfall_crypto::primitive::psf::PSFGPVPerturbation;
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
/// use qfall_math::rational::Q;
/// use qfall_crypto::primitive::psf::PSF;
///
/// let psf = PSFGPVPerturbation {
///     gp: GadgetParameters::init_default(4, 64),
///     s: Q::from(18),
///     rounding_parameter: Q::from(2)
/// };
///
/// let (a, td) = psf.trap_gen();
/// let domain_sample = psf.samp_d();
/// let range_fa = psf.f_a(&a, &domain_sample);
/// let preimage = psf.samp_p(&a, &td, &range_fa);
///
/// assert!(psf.check_domain(&preimage));
/// ```
#[derive(Serialize, Deserialize)]
pub struct PSFGPVPerturbation {
    pub gp: GadgetParameters,
    pub s: Q,
    pub rounding_parameter: Q,
}

impl PSF<MatZq, (MatZ, MatQ), MatZ, MatZq> for PSFGPVPerturbation {
    /// Computes a G-Trapdoor according to the [`GadgetParameters`]
    /// defined in the struct.
    ///
    /// Returns a parity check matrix and a corresponding G-trapdoor.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSFGPVPerturbation;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFGPVPerturbation {
    ///     gp: GadgetParameters::init_default(4, 64),
    ///     s: Q::from(18),
    ///     rounding_parameter: Q::from(2)
    /// };
    ///
    /// let (a, td) = psf.trap_gen();
    /// ```
    fn trap_gen(&self) -> (MatZq, (MatZ, MatQ)) {
        let a_bar = MatZq::sample_uniform(&self.gp.n, &self.gp.m_bar, &self.gp.q);

        let tag = MatZq::identity(&self.gp.n, &self.gp.n, &self.gp.q);

        let (a, g_trapdoor) = gen_trapdoor(&self.gp, &a_bar, &tag).unwrap();

        // use `g_trapdoor` to compute actual trapdoor: `[[td],[I]]`
        let td = g_trapdoor
            .concat_vertical(&MatZ::identity(
                a.get_num_columns() - g_trapdoor.get_num_rows(),
                g_trapdoor.get_num_columns(),
            ))
            .unwrap();

        // `s^2*I - trapdoor^t * sigma_g * trapdoor`
        // we have to subtract the trapdoor because we will sample a perturbation to a given solution
        // which is sampled using `sigma_g` and then multiplied with the trapdoor
        let sigma_p = self.s.pow(2).unwrap() * MatQ::identity(td.get_num_rows(), td.get_num_rows())
            - MatQ::from(&(&td * &self.gp.base.clone().pow(2).unwrap() * td.transpose()));

        // compute actual convolution matrix for the perturbation sampling `r * sqrt(Sigma - I) = sqrt(r^2 Sigma_p - r^2 I)`
        // we have to subtract the `r^2 I` because the randomized rounding adds additional length for which we have to account
        let convolution_matrix = &self.rounding_parameter
            * (sigma_p - MatQ::identity(a.get_num_columns(), a.get_num_columns()))
                .cholesky_decomposition_flint();

        (a, (g_trapdoor, convolution_matrix))
    }

    /// Samples in the domain using SampleD with the standard basis and center `0`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSFGPVPerturbation;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFGPVPerturbation {
    ///     gp: GadgetParameters::init_default(4, 64),
    ///     s: Q::from(18),
    ///     rounding_parameter: Q::from(2)
    /// };
    /// let (a, td) = psf.trap_gen();
    ///
    /// let domain_sample = psf.samp_d();
    /// ```
    fn samp_d(&self) -> MatZ {
        let m = &self.gp.n * &self.gp.k + &self.gp.m_bar;
        MatZ::sample_d_common(&m, &self.gp.n, &self.s * &self.rounding_parameter).unwrap()
    }

    /// Samples an `e` in the domain using SampleD with a G-trapdoor `td` for
    /// the parity-check matrix `a`.
    /// This function samples a discrete Gaussian under the condition `f_a(a,e) = u`
    /// for a provided syndrome `u` and fixed `a` with gaussian parameter
    /// `s * rounding_parameter`.
    ///
    /// This algorithm is implemented according algorithm 3 from
    /// [\[1\]](<../index.html#:~:text=[1]>).
    ///
    /// *Note*: The provided parameters `a,td,u` must fit together,
    /// otherwise unexpected behavior such as panics may occur.
    ///
    /// Parameters:
    /// - `a`: The parity-check matrix
    /// - `td`: The G-trapdoor for the parity-check matrix
    /// - `u`: The syndrome from the range
    ///
    /// Returns a sample `e` from the domain on the conditioned discrete
    /// Gaussian distribution `f_a(a,e) = u`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSFGPVPerturbation;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFGPVPerturbation {
    ///     gp: GadgetParameters::init_default(4, 64),
    ///     s: Q::from(18),
    ///     rounding_parameter: Q::from(2)
    /// };
    /// let (a, td) = psf.trap_gen();
    /// let domain_sample = psf.samp_d();
    /// let range_fa = psf.f_a(&a, &domain_sample);
    ///
    /// let preimage = psf.samp_p(&a, &td, &range_fa);
    /// assert_eq!(range_fa, psf.f_a(&a, &preimage))
    /// ```
    fn samp_p(&self, a: &MatZq, (td, convolution_matrix): &(MatZ, MatQ), u: &MatZq) -> MatZ {
        // use `td` to compute actual G-Trapdoor: `[[td],[I]]`
        let trapdoor = td
            .concat_vertical(&MatZ::identity(
                a.get_num_columns() - td.get_num_rows(),
                td.get_num_columns(),
            ))
            .unwrap();

        // we sample a perturbation value with convolution `r * sqrt(Sigma_p)`
        // the convolution matrix that is provided has to account for the additional
        // randomized rounding that has to occure, hence it is generated in the `TrapGen`
        let perturbation = MatZ::sample_d_common_non_spherical(
            &self.gp.n,
            &convolution_matrix,
            &self.rounding_parameter,
        )
        .unwrap();

        // compute the new target vector that acocunts that includes the perturbation
        let v = u - a * &perturbation;

        let z: MatZ = sampling_gaussian_gadget(
            &self.gp,
            &self.rounding_parameter * &self.gp.base,
            &v.get_mat(),
        );

        perturbation + trapdoor * z
    }

    /// Implements the efficiently computable function `f_a` which here corresponds to
    /// `a*sigma`. The sigma must be from the domain, i.e. D_n.
    ///
    /// Parameters:
    /// - `a`: The parity-check matrix of dimensions `n x m`
    /// - `sigma`: A column vector of length `m`
    ///
    /// Returns `a*sigma`
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSFGPVPerturbation;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFGPVPerturbation {
    ///     gp: GadgetParameters::init_default(4, 64),
    ///     s: Q::from(18),
    ///     rounding_parameter: Q::from(2)
    /// };
    /// let (a, td) = psf.trap_gen();
    /// let domain_sample = psf.samp_d();
    /// let range_fa = psf.f_a(&a, &domain_sample);
    /// ```
    ///
    /// # Panics ...
    /// - if `sigma` is not in the domain.
    fn f_a(&self, a: &MatZq, sigma: &MatZ) -> MatZq {
        assert!(self.check_domain(sigma));
        a * sigma
    }

    /// Checks whether a value `sigma` is in
    /// D_n = {e ∈ Z^m | |e| <= s * rounding_parameter *  sqrt(m)}.
    ///
    /// Parameters:
    /// - `sigma`: The value for which is checked, if it is in the domain
    ///
    /// Returns true, if `sigma` is in D_n.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSF;
    /// use qfall_crypto::primitive::psf::PSFGPVPerturbation;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    ///
    /// let psf = PSFGPVPerturbation {
    ///     gp: GadgetParameters::init_default(4, 64),
    ///     s: Q::from(18),
    ///     rounding_parameter: Q::from(2)
    /// };
    /// let (a, td) = psf.trap_gen();
    ///
    /// let vector = psf.samp_d();
    ///
    /// assert!(psf.check_domain(&vector));
    /// ```
    fn check_domain(&self, sigma: &MatZ) -> bool {
        let m = &self.gp.n * &self.gp.k + &self.gp.m_bar;
        sigma.is_column_vector()
            && m == Z::from(sigma.get_num_rows())
            && Q::from(&sigma.norm_eucl_sqrd().unwrap())
                <= self.s.pow(2).unwrap() * &m * &self.rounding_parameter.pow(2).unwrap()
    }
}

#[cfg(test)]
mod test_gpv_psf_perturbation {
    use super::{PSFGPVPerturbation, PSF};
    use crate::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    use qfall_math::integer::MatZ;
    use qfall_math::rational::Q;
    use qfall_math::traits::{GetNumColumns, GetNumRows, SetEntry};

    /// Ensures that `samp_d` actually computes values that are in D_n.
    #[test]
    fn samp_d_samples_from_dn() {
        for (n, q) in [(5, 256), (10, 128), (15, 157)] {
            let psf = PSFGPVPerturbation {
                gp: GadgetParameters::init_default(n, q),
                s: Q::from(20),
                rounding_parameter: Q::from(n).log(2).unwrap(),
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
        {
            let (n, q) = (4, 128);
            let psf = PSFGPVPerturbation {
                gp: GadgetParameters::init_default(n, q),
                s: Q::from(20),
                rounding_parameter: Q::from(n).log(2).unwrap(),
            };
            let (a, r) = psf.trap_gen();
            let domain_sample = psf.samp_d();
            let range_fa = psf.f_a(&a, &domain_sample);

            let preimage = psf.samp_p(&a, &r, &range_fa);
            assert_eq!(range_fa, psf.f_a(&a, &preimage));
            assert!(psf.check_domain(&preimage));
        }
    }

    /// Ensures that `f_a` returns `a*sigma`.
    #[test]
    fn f_a_works_as_expected() {
        for (n, q) in [(5, 256), (6, 128)] {
            let psf = PSFGPVPerturbation {
                gp: GadgetParameters::init_default(n, q),
                s: Q::from(20),
                rounding_parameter: Q::from(n).log(2).unwrap(),
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
        let psf = PSFGPVPerturbation {
            gp: GadgetParameters::init_default(8, 128),
            s: Q::from(20),
            rounding_parameter: Q::from(3),
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
        let psf = PSFGPVPerturbation {
            gp: GadgetParameters::init_default(8, 128),
            s: Q::from(20),
            rounding_parameter: Q::from(3),
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
        let psf = PSFGPVPerturbation {
            gp: GadgetParameters::init_default(8, 128),
            s: Q::from(20),
            rounding_parameter: Q::from(3),
        };
        let (a, _) = psf.trap_gen();
        let not_in_domain = (&psf.s * &psf.rounding_parameter).round()
            * a.get_num_columns()
            * MatZ::identity(a.get_num_columns(), 1);

        let _ = psf.f_a(&a, &not_in_domain);
    }

    /// Ensures that `check_domain` works for vectors with the correct length.
    #[test]
    fn check_domain_as_expected() {
        let psf = PSFGPVPerturbation {
            gp: GadgetParameters::init_default(4, 128),
            s: Q::from(20),
            rounding_parameter: Q::from(3),
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
        let psf = PSFGPVPerturbation {
            gp: GadgetParameters::init_default(4, 128),
            s: Q::from(20),
            rounding_parameter: Q::from(3),
        };
        let (a, _) = psf.trap_gen();

        let matrix = MatZ::new(a.get_num_columns(), 2);
        let too_short = MatZ::new(a.get_num_columns() - 1, 1);
        let too_long = MatZ::new(a.get_num_columns() + 1, 1);
        let entry_too_large = psf.rounding_parameter.round()
            * psf.s.round()
            * a.get_num_columns()
            * MatZ::identity(a.get_num_columns(), 1);

        assert!(!psf.check_domain(&matrix));
        assert!(!psf.check_domain(&too_long));
        assert!(!psf.check_domain(&too_short));
        assert!(!psf.check_domain(&entry_too_large));
    }
}
