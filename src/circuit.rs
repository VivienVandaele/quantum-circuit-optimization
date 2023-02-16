use regex::Regex;
use std::fs::{File};
use std::io::{BufRead, BufReader, Write};
use std::collections::HashMap;

#[derive(Clone)]
pub struct Circuit {
    pub circ: Vec<(String, Vec<usize>)>,
    pub nb_qubits: usize,
}

impl Circuit {
    pub fn new(nb_qubits: usize) -> Self {
        Circuit {
            circ: Vec::new(),
            nb_qubits: nb_qubits,
        }
    }

    pub fn from_qc(filename: &str) -> (Circuit, String, HashMap<usize, String>) {
        let mut c = Circuit::new(0);
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let re = Regex::new(r"\s([[:alnum:]]*)").unwrap();
        let re_gate = Regex::new(r"(\.*[[:alpha:]]+)\s").unwrap();
        let mut header: String = "".to_string();
        let mut qubits_mapping = HashMap::new();
        let mut rev_qubits_mapping = HashMap::new();
        for (_, line) in reader.lines().enumerate() {
            let line = line.unwrap(); 
            if line.len() == 0 || line.chars().next().unwrap() == '#' { continue }
            let mut gate: Vec<String> = re_gate.captures_iter(&line).map(|x| x.get(1).unwrap().as_str().parse().unwrap()).collect();
            if gate.len() == 0 { continue }
            if gate[0] == ".v" {
                for qubit in re.captures_iter(&line).map(|x| x.get(1).unwrap().as_str()) {
                    qubits_mapping.insert(qubit.to_string(), c.nb_qubits);
                    rev_qubits_mapping.insert(c.nb_qubits, qubit.to_string());
                    c.nb_qubits += 1;
                }
            }
            if gate[0].chars().next().unwrap() == '.' {
                header.push_str(&line.to_string()); 
                header.push_str("\n"); 
                continue
            }
            let qubits: Vec<usize> = re.captures_iter(&line).map(|x| *qubits_mapping.get(x.get(1).unwrap().as_str()).unwrap()).collect();
            if gate[0] == "tof" && qubits.len() == 3 { gate[0] = "tof".to_string() }
            else if (gate[0] == "Zd" || gate[0] == "Z") && qubits.len() == 3 { gate[0] = "ccz".to_string() }
            else if gate[0] == "cnot" || (gate[0] == "tof" && qubits.len() == 2) { gate[0] = "cx".to_string() }
            else if gate[0] == "H" && qubits.len() == 1 { gate[0] = "h".to_string() }
            else if gate[0] == "X" && qubits.len() == 1 { gate[0] = "x".to_string() }
            else if gate[0] == "Z" && qubits.len() == 1 { gate[0] = "z".to_string() }
            else if gate[0] == "S" && qubits.len() == 1 { gate[0] = "s".to_string() }
            else if gate[0] == "T" && qubits.len() == 1 { gate[0] = "t".to_string() }
            else { println!("Operator not implemented: {}", gate[0]); std::process::exit(1) }
            c.circ.push((gate[0].to_string(), qubits));
        }
        (c, header, rev_qubits_mapping)
    }

    pub fn to_qc(&self, filename: &str, header: String, map: HashMap<usize, String>) {
        let mut file = File::create(filename).unwrap();
        write!(file, "{}", header).unwrap();
        write!(file, "\nBEGIN\n").unwrap();
        for (gate, q) in &self.circ {
            match &gate[..] {
                "h" => write!(file, "H {}\n", map.get(&q[0]).unwrap()).unwrap(),
                "x" => write!(file, "X {}\n", map.get(&q[0]).unwrap()).unwrap(),
                "z" => write!(file, "Z {}\n", map.get(&q[0]).unwrap()).unwrap(),
                "s" => write!(file, "S {}\n", map.get(&q[0]).unwrap()).unwrap(),
                "t" => write!(file, "T {}\n", map.get(&q[0]).unwrap()).unwrap(),
                "cx" => write!(file, "cnot {} {}\n", map.get(&q[0]).unwrap(), map.get(&q[1]).unwrap()).unwrap(),
                "ccx" => write!(file, "tof {} {} {}\n", map.get(&q[0]).unwrap(), map.get(&q[1]).unwrap(), map.get(&q[2]).unwrap()).unwrap(),
                "ccz" => write!(file, "Z {} {} {}\n", map.get(&q[0]).unwrap(), map.get(&q[1]).unwrap(), map.get(&q[2]).unwrap()).unwrap(),
                _ => {println!("Operator not implemented: {}", gate); std::process::exit(1)},
            }
        }
        write!(file, "END").unwrap();
    }

    pub fn get_statistics(&self) -> (usize, usize, usize) {
        let mut h_count = 0;
        let mut internal_h_count = 0;
        let mut t_count = 0;
        let mut flag = false;
        for (gate, _) in &self.circ {
            if gate == "h" {
                h_count += 1; 
                if flag { internal_h_count += 1; }
            }
            if gate == "t" { t_count += 1; flag = true; }
        }
        if flag {
            for (gate, _) in self.circ.iter().rev() {
                if gate == "h" { internal_h_count -= 1; }
                if gate == "t" { break; }
            }
        }
        (h_count, internal_h_count, t_count)
    }

    pub fn append(&mut self, mut circ: Vec<(String, Vec<usize>)>) {
        self.circ.append(&mut circ);
    }
}
