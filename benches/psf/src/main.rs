use std::{fs::File, io::Write};

use qfall_crypto::{
    primitive::psf::{PSFGPVPerturbation, PSF},
    sample::g_trapdoor::gadget_parameters::GadgetParameters,
};
use qfall_math::{integer::Z, rational::Q, traits::Pow};

/// Run this code to generate an instance for the psf that can be used for benchmarking.
///
/// Running this code with the current parameter set takes around 20 minutes.
fn main() {
    let (n, q) = (64, Z::from(2).pow(24).unwrap());

    let psf = PSFGPVPerturbation {
        gp: GadgetParameters::init_default(n, q),
        s: Q::from(418),
        rounding_parameter: Q::from(4.5),
    };

    let (a, (r, convolution_matrix)) = psf.trap_gen();

    let a_string = serde_json::to_string(&a).unwrap();
    let r_string = serde_json::to_string(&r).unwrap();
    let conv_matrix_string = serde_json::to_string(&convolution_matrix).unwrap();

    let mut file_a = File::create("parity_check_matrix.txt").unwrap();
    let mut file_r = File::create("trapdoor.txt").unwrap();
    let mut file_conv = File::create("convolution_matrix.txt").unwrap();
    file_a.write_all(a_string.as_bytes()).unwrap();
    file_r.write_all(r_string.as_bytes()).unwrap();
    file_conv.write_all(conv_matrix_string.as_bytes()).unwrap();
}
