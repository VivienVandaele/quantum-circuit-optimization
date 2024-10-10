mod circuit;
mod pauli_product;
mod phase_polynomial;
mod tableau;
mod bit_vector;
mod h_opt;
mod t_merge;
mod t_opt;
use crate::h_opt::internal_h_opt;
use crate::circuit::Circuit;
use crate::t_merge::*;
use std::path::Path;

fn help() {
    println!("cargo run -r [OPTIONS] file.qc\n\nOptional arguments (case-insensitive, no order):");
    println!("'BBMerge': runs the BBMerge algorithm");
    println!("'FastTMerge': runs the FastTMerge algorithm");
    println!("'InternalHOpt': runs the InternalHOpt algorithm");
    println!("'TOHPE': runs the TOHPE algorithm");
    println!("'FastTODD': runs the FastTODD algorithm");
    std::process::exit(1);
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args.iter().any(|s| s.to_lowercase().ends_with("help")) { help(); }
    let file_index = args.iter().position(|s| s.ends_with(".qc"));
    if file_index == None {
        println!("No .qc file provided");
        help();
    }

    let do_bb_merge = args.iter().any(|s| s.to_lowercase().ends_with("bbmerge"));
    let mut do_fast_t_merge = args.iter().any(|s| s.to_lowercase().ends_with("fasttmerge"));
    let mut do_internal_h_opt = args.iter().any(|s| s.to_lowercase().ends_with("internalhopt"));
    let do_tohpe = args.iter().any(|s| s.to_lowercase().ends_with("tohpe"));
    let mut do_fast_todd = args.iter().any(|s| s.to_lowercase().ends_with("fasttodd"));

    if !(do_bb_merge || do_fast_t_merge || do_internal_h_opt || do_tohpe || do_fast_todd) {
        do_fast_t_merge = true;
        do_internal_h_opt = true;
        do_fast_todd = true;
    }

    let filename = Path::new(&args[file_index.unwrap()]).file_name().unwrap().to_str().unwrap();
    let output_filename = &("circuits/outputs/".to_string() + &filename);
    let (mut c, header, qubits_mapping) = Circuit::from_qc(&args[file_index.unwrap()]);
    println!("File {} processed\n", filename);
     {
        if do_bb_merge { println!("Running BBMerge algorithm"); c = bb_merge(c); }
        if do_fast_t_merge { println!("Running FastTMerge algorithm"); c = fast_t_merge(c); }
        if do_internal_h_opt { println!("Running InternalHOpt algorithm"); c = internal_h_opt(&c); }
        if do_tohpe || do_fast_todd { println!("Internal Hadamard gates gadgetization"); c = c.hadamard_gadgetization(); }
        if do_tohpe { println!("Running TOHPE algorithm"); c = c.t_opt("TOHPE".to_string()); }
        if do_fast_todd { println!("Running FastTODD algorithm"); c = c.t_opt("FastTODD".to_string()); }

        let (h_count, internal_h_count, t_count) = c.get_statistics();
        println!("\nOptimized circuit:\nH-count: {}\nInternal H-count: {}\nT-count: {}", h_count, internal_h_count, t_count);
        c.to_qc(output_filename, header, qubits_mapping);
    }
}
