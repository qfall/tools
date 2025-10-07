// Copyright © 2025 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

use criterion::*;
use qfall_crypto::construction::pk_encryption::PKEncryptionScheme;
use qfall_crypto::construction::pk_encryption::KPKE;

/// Performs a full-cycle of gen, enc, dec with [`KPKE`].
fn kpke_cycle(k_pke: &KPKE) {
    let (pk, sk) = k_pke.gen();
    let cipher = k_pke.enc(&pk, 1);
    let _ = k_pke.dec(&sk, &cipher);
}

/// Benchmark [kpke_cycle] with [KPKE::ml_kem_512].
///
/// This benchmark can be run with for example:
/// - `cargo criterion K-PKE\ cycle\ 512`
/// - `cargo bench --bench benchmarks K-PKE\ cycle\ 512`
/// - `cargo flamegraph --bench benchmarks -- --bench K-PKE\ cycle\ 512`
fn bench_kpke_cycle_512(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_512();

    c.bench_function("K-PKE cycle 512", |b| b.iter(|| kpke_cycle(&k_pke)));
}

/// Benchmark [KPKE::gen] with [KPKE::ml_kem_512].
fn bench_kpke_gen_512(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_512();

    c.bench_function("K-PKE gen 512", |b| b.iter(|| k_pke.gen()));
}

/// Benchmark [KPKE::enc] with [KPKE::ml_kem_512].
fn bench_kpke_enc_512(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_512();
    let (pk, _) = k_pke.gen();
    let msg = i64::MAX;

    c.bench_function("K-PKE enc 512", |b| b.iter(|| k_pke.enc(&pk, msg)));
}

/// Benchmark [KPKE::dec] with [KPKE::ml_kem_512].
fn bench_kpke_dec_512(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_512();
    let (pk, sk) = k_pke.gen();
    let cipher = k_pke.enc(&pk, i64::MAX);

    c.bench_function("K-PKE dec 512", |b| b.iter(|| k_pke.dec(&sk, &cipher)));
}

/// Benchmark [kpke_cycle] with [KPKE::ml_kem_768].
///
/// This benchmark can be run with for example:
/// - `cargo criterion K-PKE\ cycle\ 768`
/// - `cargo bench --bench benchmarks K-PKE\ cycle\ 768`
/// - `cargo flamegraph --bench benchmarks -- --bench K-PKE\ cycle\ 768`
fn bench_kpke_cycle_768(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_768();

    c.bench_function("K-PKE cycle 768", |b| b.iter(|| kpke_cycle(&k_pke)));
}

/// Benchmark [KPKE::gen] with [KPKE::ml_kem_768].
fn bench_kpke_gen_768(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_768();

    c.bench_function("K-PKE gen 768", |b| b.iter(|| k_pke.gen()));
}

/// Benchmark [KPKE::enc] with [KPKE::ml_kem_768].
fn bench_kpke_enc_768(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_768();
    let (pk, _) = k_pke.gen();
    let msg = i64::MAX;

    c.bench_function("K-PKE enc 768", |b| b.iter(|| k_pke.enc(&pk, msg)));
}

/// Benchmark [KPKE::dec] with [KPKE::ml_kem_768].
fn bench_kpke_dec_768(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_768();
    let (pk, sk) = k_pke.gen();
    let cipher = k_pke.enc(&pk, i64::MAX);

    c.bench_function("K-PKE dec 768", |b| b.iter(|| k_pke.dec(&sk, &cipher)));
}

/// Benchmark [kpke_cycle] with [KPKE::ml_kem_1024].
///
/// This benchmark can be run with for example:
/// - `cargo criterion K-PKE\ cycle\ 1024`
/// - `cargo bench --bench benchmarks K-PKE\ cycle\ 1024`
/// - `cargo flamegraph --bench benchmarks -- --bench K-PKE\ cycle\ 1024`
fn bench_kpke_cycle_1024(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_1024();

    c.bench_function("K-PKE cycle 1024", |b| b.iter(|| kpke_cycle(&k_pke)));
}

/// Benchmark [KPKE::gen] with [KPKE::ml_kem_1024].
fn bench_kpke_gen_1024(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_1024();

    c.bench_function("K-PKE gen 1024", |b| b.iter(|| k_pke.gen()));
}

/// Benchmark [KPKE::enc] with [KPKE::ml_kem_1024].
fn bench_kpke_enc_1024(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_1024();
    let (pk, _) = k_pke.gen();
    let msg = i64::MAX;

    c.bench_function("K-PKE enc 1024", |b| b.iter(|| k_pke.enc(&pk, msg)));
}

/// Benchmark [KPKE::dec] with [KPKE::ml_kem_1024].
fn bench_kpke_dec_1024(c: &mut Criterion) {
    let k_pke = KPKE::ml_kem_1024();
    let (pk, sk) = k_pke.gen();
    let cipher = k_pke.enc(&pk, i64::MAX);

    c.bench_function("K-PKE dec 1024", |b| b.iter(|| k_pke.dec(&sk, &cipher)));
}

criterion_group!(
    benches,
    bench_kpke_cycle_512,
    bench_kpke_gen_512,
    bench_kpke_enc_512,
    bench_kpke_dec_512,
    bench_kpke_cycle_768,
    bench_kpke_gen_768,
    bench_kpke_enc_768,
    bench_kpke_dec_768,
    bench_kpke_cycle_1024,
    bench_kpke_gen_1024,
    bench_kpke_enc_1024,
    bench_kpke_dec_1024,
);
