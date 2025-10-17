# qFALL-tools

[![made-with-rust](https://img.shields.io/badge/Made%20with-Rust-1f425f.svg)](https://www.rust-lang.org/)
[![CI](https://github.com/qfall/tools/actions/workflows/push.yml/badge.svg?branch=dev)](https://github.com/qfall/tools/actions/workflows/pull_request.yml)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL_2.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)

This repository is currently being developed by the project group [qFALL - quantum resistant fast lattice library](https://cs.uni-paderborn.de/cuk/lehre/veranstaltungen/ws-2022-23/project-group-qfall) in the winter term 2022 and summer term 2023 by the Codes and Cryptography research group in Paderborn.

The main objective of this project is to provide researchers and students with the possibility to easily and quickly prototype (lattice-based) cryptography.

## Disclaimer

Currently, we are in the development phase and interfaces might change.
Feel free to check out the current progress, but be aware, that the content will
change in the upcoming weeks and months. An official release will most likely be published in the second half of 2023.

## Quick-Start

Please refer to [our website](https://qfall.github.io/) as central information point.

To install and add our library to your project, please refer to [our tutorial](https://qfall.github.io/book/index.html).
It provides a step-by-step guide to install the required libraries and gives further insights in the usage of our crates.

## What does qFALL-tools offer?

qFALL-tools offers a variety of implementations of commonly used tools in lattice-based cryptography.
We provide a brief overview in the following list.
For a more detailed description, please refer to [our tutorial section](https://qfall.github.io/book/crypto/features.html).

- [Preimage Samplable Functions (PSF)](https://github.com/qfall/tools/blob/dev/src/primitive/psf.rs)
- [Trapdoors](https://github.com/qfall/tools/blob/dev/src/sample/g_trapdoor.rs)
  - [G-trapdoor incl. short basis](https://github.com/qfall/tools/blob/dev/src/sample/g_trapdoor/gadget_classical.rs)
  - [Ring-based G-trapdoor incl. short basis](https://github.com/qfall/tools/blob/dev/src/sample/g_trapdoor/gadget_ring.rs)
- [Utility functions for quick instantiations](https://github.com/qfall/tools/blob/dev/src/utils/)
  - [Common moduli](https://github.com/qfall/tools/blob/dev/src/utils/common_moduli.rs)
  - [Rotation matrices](https://github.com/qfall/tools/blob/dev/src/utils/rotation_matrix.rs)
  - [Common encodings](https://github.com/qfall/tools/blob/dev/src/utils/common_encodings.rs)

## License

This library is distributed under the **Mozilla Public License Version 2.0** which can be found here [License](https://github.com/qfall/tools/blob/dev/LICENSE).
Permissions of this weak copyleft license are conditioned on making available source code of licensed files and modifications of those files under the same license (or in certain cases, one of the GNU licenses). Copyright and license notices must be preserved. Contributors provide an express grant of patent rights. However, a larger work using the licensed work may be distributed under different terms and without source code for files added in the larger work.

## Citing

Please use the following bibtex entry to cite [qFALL-tools](https://github.com/qfall/tools):

```text
@misc{qFALL-tools,
    author = {Porzenheim, Laurens and Beckmann, Marvin and Kramer, Paul and Milewski, Phil and Moog, Sven and Schmidt, Marcel and Siemer, Niklas},
    title = {qFALL-tools v0.0},
    howpublished = {Online: \url{https://github.com/qfall/tools}},
    month = Mar,
    year = 2023,
    note = {University Paderborn,  Codes and Cryptography}
}
```

## Get in Touch

One can contact the members of the project group with our mailing list `pg-qfall(at)lists.upb.de`.
