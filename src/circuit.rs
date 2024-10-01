use regex::Regex;
use std::fs::{File};
use std::io::{BufRead, BufReader, Write};
use std::collections::HashMap;
use crate::phase_polynomial::PhasePolynomial;
use crate::tableau::TableauColumnMajor;
use crate::t_opt::{tohpe, fast_todd};

#[derive(Debug, Clone)]
pub struct Circuit {
    pub circ: Vec<(String, Vec<usize>)>,
    pub nb_qubits: usize,
    pub ancillas: HashMap::<usize, usize>,
}

impl Circuit {
    pub fn new(nb_qubits: usize) -> Self {
        Circuit {
            circ: Vec::new(),
            nb_qubits: nb_qubits,
            ancillas: HashMap::new(),
        }
    }

    pub fn from_qc(filename: &str) -> (Circuit, String, HashMap<usize, String>) {
        let mut c = Circuit::new(0);
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let re = Regex::new(r"\s([[:alnum:]]*)").unwrap();
        let re_gate = Regex::new(r"(\.*[[:alpha:]]+\*?)\s").unwrap();
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
            else if (gate[0] == "S" || gate[0] == "P") && qubits.len() == 1 { gate[0] = "s".to_string() }
            else if (gate[0] == "S*" || gate[0] == "P*") && qubits.len() == 1 { 
                c.circ.push(("z".to_string(), qubits.clone()));
                gate[0] = "s".to_string();
            }
            else if gate[0] == "T" && qubits.len() == 1 { gate[0] = "t".to_string() }
            else if gate[0] == "T*" && qubits.len() == 1 { 
                c.circ.push(("z".to_string(), qubits.clone()));
                c.circ.push(("s".to_string(), qubits.clone()));
                gate[0] = "t".to_string();
            }
            else { println!("Operator not implemented: {}", gate[0]); std::process::exit(1) }
            c.circ.push((gate[0].to_string(), qubits));
        }
        (c, header, rev_qubits_mapping)
    }

    pub fn decompose_tof(&self) -> Circuit {
        let mut c = Circuit::new(self.nb_qubits);
        for (gate, qubits) in &self.circ {
            if (gate == "ccz" || gate == "tof") && qubits.len() == 3 { 
                if gate == "tof" {
                    c.circ.push(("h".to_string(), vec![qubits[2]]));
                }
                for i in 0..3 {
                    c.circ.push(("t".to_string(), vec![qubits[i]]));
                }
                c.circ.push(("cx".to_string(), vec![qubits[1], qubits[0]]));
                c.circ.push(("x".to_string(), vec![qubits[0]]));
                c.circ.push(("t".to_string(), vec![qubits[0]]));
                c.circ.push(("x".to_string(), vec![qubits[0]]));
                c.circ.push(("cx".to_string(), vec![qubits[2], qubits[0]]));
                c.circ.push(("t".to_string(), vec![qubits[0]]));
                c.circ.push(("cx".to_string(), vec![qubits[1], qubits[0]]));
                c.circ.push(("x".to_string(), vec![qubits[0]]));
                c.circ.push(("t".to_string(), vec![qubits[0]]));
                c.circ.push(("x".to_string(), vec![qubits[0]]));
                c.circ.push(("cx".to_string(), vec![qubits[2], qubits[0]]));
                c.circ.push(("cx".to_string(), vec![qubits[2], qubits[1]]));
                c.circ.push(("x".to_string(), vec![qubits[1]]));
                c.circ.push(("t".to_string(), vec![qubits[1]]));
                c.circ.push(("x".to_string(), vec![qubits[1]]));
                c.circ.push(("cx".to_string(), vec![qubits[2], qubits[1]]));
                if gate == "tof" {
                    c.circ.push(("h".to_string(), vec![qubits[2]]));
                }
                continue
            }
            else {
                c.circ.push((gate.to_string(), qubits.to_vec()));
            }
        }
        c
    }

    pub fn to_qc(&self, filename: &str, header: String, mut map: HashMap<usize, String>) {
        let mut file = File::create(filename).unwrap();
        let mut index = map.len();
        let mut val = map.len();
        for s in header.split("\n") {
            write!(file, "{}", s).unwrap();
            let s2 = s.split(" ").collect::<Vec<_>>();
            if s2[0] == ".v"  {
                for _ in self.ancillas.keys() {
                    while map.values().any(|x| *x == val.to_string()) {
                        val += 1;
                    }
                    write!(file, " {}", val).unwrap();
                    map.insert(index, val.to_string());
                    index += 1;
                }
            }
            write!(file, "\n").unwrap();
        }
        write!(file, "BEGIN\n").unwrap();
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

    pub fn hadamard_gadgetization(&self) -> Circuit {
        let mut c = Circuit::new(self.nb_qubits);
        let mut anc = Circuit::new(self.nb_qubits);
        let mut flag = false;
        let mut last = 0;
        let mut parent_ancilla = Vec::new();
        for i in 0..self.nb_qubits {
            parent_ancilla.push(i);
        }
        for (i, (gate, _)) in self.circ.clone().into_iter().enumerate() {
            if gate == "t" { last = i; }
        }
        for (i, (gate, qubits)) in self.circ.clone().into_iter().enumerate() {
            if gate == "t" { flag = true; }
            if gate == "h" && i < last && flag {
                anc.circ.push(("h".to_string(), vec![anc.nb_qubits]));
                c.circ.push(("s".to_string(), vec![anc.nb_qubits]));
                c.circ.push(("s".to_string(), qubits.to_vec()));
                c.circ.push(("cx".to_string(), vec![qubits[0], anc.nb_qubits]));
                c.circ.push(("s".to_string(), vec![anc.nb_qubits]));
                c.circ.push(("z".to_string(), vec![anc.nb_qubits]));
                c.circ.push(("cx".to_string(), vec![anc.nb_qubits, qubits[0]]));
                c.circ.push(("cx".to_string(), vec![qubits[0], anc.nb_qubits]));
                anc.ancillas.insert(anc.nb_qubits, parent_ancilla[qubits[0]]);
                parent_ancilla[qubits[0]] = anc.nb_qubits;
                anc.nb_qubits += 1;
            }
            else {
                c.circ.push((gate.to_string(), qubits.to_vec()));
            }
        }
        let mut c_out = anc.clone();
        c_out.append(c.circ);
        c_out.append(anc.circ);
        c_out
    }

    pub fn t_opt(&self, optimizer: String) -> Circuit {
        SlicedCircuit::from_circ(self).t_opt(optimizer)
    }
}


#[derive(Debug, Clone)]
pub struct SlicedCircuit {
    pub nb_qubits: usize,
    pub init_circuit: Circuit,
    pub tableau_vec: Vec<TableauColumnMajor>,
    pub phase_polynomials: Vec<PhasePolynomial>,
}

impl SlicedCircuit {
    pub fn new(nb_qubits: usize) -> Self {
        SlicedCircuit {
            nb_qubits: nb_qubits,
            init_circuit: Circuit::new(nb_qubits),
            tableau_vec: Vec::new(),
            phase_polynomials: Vec::new(),
        }
    }

    pub fn from_circ(c: &Circuit) -> SlicedCircuit {
        let mut sliced_c = SlicedCircuit::new(c.nb_qubits);
        sliced_c.init_circuit.ancillas = c.ancillas.clone();
        let mut first_t = 0;
        for (i, (gate, q)) in c.circ.clone().into_iter().enumerate() {
            if gate == "t" {
                first_t = i;
                break;
            }
            sliced_c.init_circuit.circ.push((gate.to_string(), q.to_vec()));
        }
        let mut tab = TableauColumnMajor::new(c.nb_qubits);
        let mut p = PhasePolynomial::new(c.nb_qubits);
        for (i, (gate, q)) in c.circ.clone().into_iter().enumerate() {
            if i < first_t { continue; }
            match &gate[..] {
                "h" => { 
                    if p.table.len() > 0 {
                        sliced_c.phase_polynomials.push(p);
                        p = PhasePolynomial::new(c.nb_qubits);
                    }
                    tab.prepend_h(q[0]);
                },
                "x" => { tab.prepend_x(q[0]); },
                "z" => { tab.prepend_z(q[0]); },
                "s" => { tab.prepend_s(q[0]); tab.prepend_z(q[0]); },
                "cx" => { tab.prepend_cx(q.to_vec()); },
                "t" => { 
                    if p.table.len() == 0 && sliced_c.phase_polynomials.len() > 0 {
                        sliced_c.tableau_vec.push(tab);
                        tab = TableauColumnMajor::new(c.nb_qubits);
                    }
                    p.table.push(tab.stabs[q[0]].z.clone());
                    if tab.stabs[q[0]].sign {
                        tab.prepend_s(q[0]);
                        tab.prepend_z(q[0]);
                    }
                },
                _ => {println!("Operator not implemented: {}", gate); std::process::exit(1)},
            }
        }
        if p.table.len() > 0 {
            sliced_c.phase_polynomials.push(p);
        }
        sliced_c.tableau_vec.push(tab);
        sliced_c
    }

    pub fn t_opt(&mut self, optimizer: String) -> Circuit {
        let mut c = self.init_circuit.clone();
        for i in 0..self.phase_polynomials.len() {
            let table = self.phase_polynomials[i].table.clone();
            if optimizer == "FastTODD" {
                self.phase_polynomials[i].table = fast_todd(table.clone(), self.nb_qubits);
            }
            else if optimizer == "TOHPE" {
                self.phase_polynomials[i].table = tohpe(table.clone(), self.nb_qubits);
            }
            else {
                println!("Optimizer not implemented: {}", optimizer); std::process::exit(1);
            }
            c.append(self.phase_polynomials[i].clifford_correction(&table, self.nb_qubits).to_circ(false).circ);
            c.append(self.phase_polynomials[i].to_circ().circ);
            if self.tableau_vec.len() > i {
                c.append(self.tableau_vec[i].to_circ(true).circ);
            }
        }
        c
    }
}
