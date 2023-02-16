#![feature(stdsimd)]
mod circuit;
mod pauli_product;
mod tableau;
mod bit_vector;
mod h_opt;
use crate::h_opt::internal_h_opt;
use crate::circuit::Circuit;
use std::path::Path;

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args.len() < 2 {
        println!("No input file provided");
        std::process::exit(1);
    }
    let path = Path::new(&args[1]);
    let filename = path.file_name().unwrap().to_str().unwrap();
    let output_filename = &("circuits/outputs/".to_string() + &filename);
    let (circ, header, qubits_mapping) = Circuit::from_qc(&args[1]);
    println!("File {} processed\nStarting optimization", filename);
    unsafe {
        let start = std::time::Instant::now();
        let c = internal_h_opt(&circ);
        println!("Time: {} seconds", start.elapsed().as_millis() as f32/1000.0);
        let (h_count, internal_h_count, t_count) = c.get_statistics();
        println!("H-count: {}\nInternal H-count: {}\nT-count: {}", h_count, internal_h_count, t_count);
        c.to_qc(output_filename, header, qubits_mapping);
    }
}
