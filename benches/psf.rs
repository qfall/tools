// Copyright © 2024 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

use criterion::{criterion_group, Criterion};
use qfall_crypto::{
    primitive::psf::{PSFGPVPerturbation, PSF, PSFGPV},
    sample::g_trapdoor::gadget_parameters::GadgetParameters,
};
use qfall_math::{integer_mod_q::MatZq, rational::Q};
use std::str::FromStr;

/// Benchmark [bench_psf] with `n = 4`.
///
/// This benchmark can be run with for example:
/// - `cargo criterion PSF\ without\ Perturbation\ n=4`
/// - `cargo bench --bench benchmarks PSF\ without\ Perturbation\ n=4`
/// - `cargo flamegraph --bench benchmarks -- --bench PSF\ without\ Perturbation\ n=4`
///
/// Shorter variants or regex expressions can also be used to specify the
/// benchmark name. The `\ ` is used to escape the space, alternatively,
/// quotation marks can be used.
fn bench_psf(c: &mut Criterion) {
    let (n, q) = (4, 128);

    let psf = PSFGPV {
        gp: GadgetParameters::init_default(n, q),
        // multiply with the rounding parameter from next test to have the same samplign parameter
        s: Q::from(20) * Q::from(n).log(2).unwrap(),
    };

    let target = MatZq::from_str("[[1],[2],[3],[4]] mod 128").unwrap();
    let (a, r) = psf.trap_gen();

    c.bench_function("PSF without Perturbation n=4", |b| {
        b.iter(|| psf.samp_p(&a, &r, &target))
    });
}

/// Benchmark [bench_psf_perturbation] with `n = 4`.
///
/// This benchmark can be run with for example:
/// - `cargo criterion PSF\ Perturbation\ n=4`
/// - `cargo bench --bench benchmarks PSF\ Perturbation\ n=4`
/// - `cargo flamegraph --bench benchmarks -- --bench PSF\ Perturbation\ n=4`
///
/// Shorter variants or regex expressions can also be used to specify the
/// benchmark name. The `\ ` is used to escape the space, alternatively,
/// quotation marks can be used.
fn bench_psf_perturbation(c: &mut Criterion) {
    let (n, q) = (4, 128);

    let psf = PSFGPVPerturbation {
        gp: GadgetParameters::init_default(n, q),
        s: Q::from(20),
        rounding_parameter: Q::from(n).log(2).unwrap(),
    };

    let target = MatZq::from_str("[[1],[2],[3],[4]] mod 128").unwrap();
    let (a, r) = psf.trap_gen();

    c.bench_function("PSF Perturbation n=4", |b| {
        b.iter(|| psf.samp_p(&a, &r, &target))
    });
}

criterion_group!(benches, bench_psf, bench_psf_perturbation);
