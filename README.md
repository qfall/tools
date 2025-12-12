# qFALL-tools
[<img alt="github" src="https://img.shields.io/badge/qfall--tools-github?style=for-the-badge&logo=github&label=github&color=8da0cb" height="20">](https://github.com/qfall/tools)
[<img alt="crates.io" src="https://img.shields.io/badge/qfall--tools-cratesio?style=for-the-badge&logo=rust&label=crates&color=fc8d62" height="20">](https://crates.io/crates/qfall-tools)
[<img alt="docs.rs" src="https://img.shields.io/badge/qfall--tools-docs?style=for-the-badge&logo=docs.rs&label=docs.rs&color=66c2a5" height="20">](https://docs.rs/qfall-tools)
[<img alt="tutorial" src="https://img.shields.io/badge/book-tutorial?style=for-the-badge&logo=mdBook&label=Tutorial&color=ffd92f" height="20">](https://qfall.github.io/book)
[<img alt="build" src="https://img.shields.io/github/actions/workflow/status/qfall/tools/push.yml?branch=main&style=for-the-badge" height="20">](https://github.com/qfall/tools/actions/workflows/push.yml)
[<img alt="license" src="https://img.shields.io/badge/License-MPL_2.0-blue.svg?style=for-the-badge" height="20">](https://github.com/qfall/tools/blob/dev/LICENSE)

`qFALL` is a prototyping library for lattice-based cryptography.
This `tools`-crate collects common sub-modules and features used by lattice-based constructions to simplify and accelerate the development of such.

## Quick-Start
First, ensure that you use a Unix-like distribution (Linux or MacOS). Setup [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) if you're using Windows. This is required due to this crate's dependency on FLINT.
Then, make sure your `rustc --version` is `1.85` or newer. 

Furthermore, it's required that `m4`, a C-compiler such as `gcc`, and `make` are installed.
```bash
sudo apt-get install m4 gcc make
```
Then, add you can add this crate to your project by executing the following command.
```bash
cargo add qfall-tools
```
- Find further information on [our website](https://qfall.github.io/). Also check out [`qfall-math`](https://crates.io/crates/qfall-math) and [`qfall-schemes`](https://crates.io/crates/qfall-schemes).
- Read the [documentation of this crate](https://docs.rs/qfall-tools).
- We recommend [our tutorial](https://qfall.github.io/book) to start working with qFALL.

## What does qFALL-tools offer?
qFALL-tools offers several commonly used sub-modules specific to lattice-based cryptography.
- [Preimage Samplable Functions (PSF)](https://docs.rs/qfall-tools/latest/qfall_tools/primitive/psf/index.html)
  - [GPV-based PSF over Z_q](https://docs.rs/qfall-tools/latest/qfall_tools/primitive/psf/struct.PSFGPV.html)
  - [GPV-based PSF over R_q](https://docs.rs/qfall-tools/latest/qfall_tools/primitive/psf/struct.PSFGPVRing.html)
  - [MP12 / Perturbation-based PSF over Z_q](https://docs.rs/qfall-tools/latest/qfall_tools/primitive/psf/struct.PSFPerturbation.html)
- [Trapdoors](https://docs.rs/qfall-tools/latest/qfall_tools/sample/g_trapdoor/index.html)
  - [G-trapdoor incl. short basis](https://docs.rs/qfall-tools/latest/qfall_tools/sample/g_trapdoor/gadget_classical/index.html)
  - [Ring-based G-trapdoor incl. short basis](https://docs.rs/qfall-tools/latest/qfall_tools/sample/g_trapdoor/gadget_ring/index.html)

Furthermore, this crate simplifies the implementation of your prototype by supporting a range of utility functions to quickly instantiate commonly used moduli, matrices, or encodings.
- [Utility functions](https://docs.rs/qfall-tools/latest/qfall_tools/utils/index.html)
  - [Common encodings](https://docs.rs/qfall-tools/latest/qfall_tools/utils/common_encodings/index.html)
  - [Common moduli](https://docs.rs/qfall-tools/latest/qfall_tools/utils/common_moduli/index.html)
  - [Lossy Compression](https://docs.rs/qfall-tools/latest/qfall_tools/utils/lossy_compression/index.html)
  - [Rotation matrices](https://docs.rs/qfall-tools/latest/qfall_tools/utils/rotation_matrix/index.html)

## Quick Examples
From String to Encoding for Encryption
```rust
use qfall_tools::utils::{common_moduli::new_anticyclic, common_encodings::encode_value_in_polynomialringzq};
use qfall_math::integer::Z;

// Create X^256 + 1 mod 3329
let poly_mod = new_anticyclic(256, 3329).unwrap();
// Generate integer from string
let message = Z::from_utf8("Hello!");
// Turn string into encoding q/2 and 0 for each 1 and 0 bit respectively
let mu_q_half = encode_value_in_polynomialringzq(message, 2, &poly_mod).unwrap();
```

Preimage Sampling using a PSF
```rust
use qfall_tools::primitive::psf::{PSF, PSFPerturbation};
use qfall_tools::sample::g_trapdoor::gadget_parameters::GadgetParameters;
use qfall_math::rational::Q;

let psf = PSFPerturbation {
    gp: GadgetParameters::init_default(8, 64),
    r: Q::from(3),
    s: Q::from(25),
};

// Generate matrix A with a trapdoor
let (a, td) = psf.trap_gen();
// Choose a random target
let domain_sample = psf.samp_d();
let target = psf.f_a(&a, &domain_sample);
// Sample a preimage for the given target
let preimage = psf.samp_p(&a, &td, &target);

assert!(psf.check_domain(&preimage));
assert_eq!(a * preimage, target);
```

## Bugs
Please report bugs through the [GitHub issue tracker](https://github.com/qfall/tools/issues).

## Contributions
Contributors are:
- Marvin Beckmann
- Jan Niklas Siemer

See [Contributing](https://github.com/qfall/tools/blob/dev/CONTRIBUTING.md) for details on how to contribute.

## Citing

Please use the following bibtex entry to cite [qFALL](https://qfall.github.io).

```text
TODO: Update to eprint
```

## Dependencies
This project is based on [qfall-math](https://crates.io/crates/qfall-math), which builds on top of the C-based, optimised math-library [FLINT](https://flintlib.org/). We utilise [serde](https://crates.io/crates/serde) and [serde_json](https://crates.io/crates/serde_json) to (de-)serialize objects to and from JSON. This crate relies on [criterion](https://crates.io/crates/criterion) for benchmarking purposes. An extensive list can be found in our `Cargo.toml` file.

## License

This library is distributed under the [Mozilla Public License Version 2.0](https://github.com/qfall/tools/blob/dev/LICENSE).
Permissions of this weak copyleft license are conditioned on making the source code of licensed files and modifications of those files available under the same license (or in certain cases, under one of the GNU licenses). Copyright and license notices must be preserved. Contributors provide an express grant of patent rights. However, a larger work using the licensed work may be distributed under different terms and without source code for files added to the larger work.
