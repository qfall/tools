// Copyright © 2025 Niklas Siemer
//
// This file is part of qFALL-tools.
//
// qFALL-tools is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

use criterion::{criterion_group, Criterion};
use qfall_math::{integer_mod_q::MatZq, rational::Q};
use qfall_tools::{
    primitive::psf::{PSFPerturbation, PSF, PSFGPV},
    sample::g_trapdoor::gadget_parameters::GadgetParameters,
};

/// Benchmark [bench_psf] with `n = 8`.
///
/// This benchmark can be run with for example:
/// - `cargo criterion PSF\ GPV\ n=8`
/// - `cargo bench --bench benchmarks PSF\ GPV\ n=8`
/// - `cargo flamegraph --bench benchmarks -- --bench PSF\ GPV\ n=8`
///
/// Shorter variants or regex expressions can also be used to specify the
/// benchmark name. The `\ ` is used to escape the space, alternatively,
/// quotation marks can be used.
fn bench_psf(c: &mut Criterion) {
    let (n, q) = (8, 128);

    let psf = PSFGPV {
        gp: GadgetParameters::init_default(n, q),
        // multiply with the rounding parameter from next test to have the same samplign parameter
        s: Q::from(30) * Q::from(n).log(2).unwrap(),
    };

    let target = MatZq::sample_uniform(n, 1, q);
    let (a, r) = psf.trap_gen();

    c.bench_function("PSF GPV n=8", |b| b.iter(|| psf.samp_p(&a, &r, &target)));
}

/// Benchmark [bench_psf_perturbation] with `n = 8`.
///
/// This benchmark can be run with for example:
/// - `cargo criterion PSF\ Perturbation\ n=8`
/// - `cargo bench --bench benchmarks PSF\ Perturbation\ n=8`
/// - `cargo flamegraph --bench benchmarks -- --bench PSF\ Perturbation\ n=8`
///
/// Shorter variants or regex expressions can also be used to specify the
/// benchmark name. The `\ ` is used to escape the space, alternatively,
/// quotation marks can be used.
fn bench_psf_perturbation(c: &mut Criterion) {
    let (n, q) = (8, 128);

    let psf = PSFPerturbation {
        gp: GadgetParameters::init_default(n, q),
        s: Q::from(30),
        r: Q::from(n).log(2).unwrap(),
    };

    let target = MatZq::sample_uniform(n, 1, q);
    let (a, r) = psf.trap_gen();

    c.bench_function("PSF Perturbation n=8", |b| {
        b.iter(|| psf.samp_p(&a, &r, &target))
    });
}

/// Benchmark [bench_psf_perturbation] with `n = 64`.
///
/// This benchmark can be run with for example:
/// - `cargo criterion PSF\ Perturbation\ n=64`
/// - `cargo bench --bench benchmarks PSF\ Perturbation\ n=64`
/// - `cargo flamegraph --bench benchmarks -- --bench PSF\ Perturbation\ n=64`
///
/// Shorter variants or regex expressions can also be used to specify the
/// benchmark name. The `\ ` is used to escape the space, alternatively,
/// quotation marks can be used.
fn bench_psf_perturbation_larger(c: &mut Criterion) {
    let (n, q) = (64, 128);

    let psf = PSFPerturbation {
        gp: GadgetParameters::init_default(n, q),
        s: Q::from(100),
        r: Q::from(n).log(2).unwrap(),
    };

    let target = MatZq::sample_uniform(n, 1, q);
    let (a, r) = psf.trap_gen();

    c.bench_function("PSF Perturbation n=64", |b| {
        b.iter(|| psf.samp_p(&a, &r, &target))
    });
}

criterion_group!(
    benches,
    bench_psf,
    bench_psf_perturbation,
    bench_psf_perturbation_larger,
);
