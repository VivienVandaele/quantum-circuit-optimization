use crate::tableau::{Tableau, TableauColumnMajor};
use crate::circuit::Circuit;
use std::collections::HashMap;

pub fn bb_merge(c_in: Circuit) -> Circuit {
    let nb_qubits = c_in.nb_qubits;
    let v = rank_vector(&c_in);
    let mut r = vec![1; v.len()];
    let mut tab = TableauColumnMajor::new(nb_qubits);
    let mut pauli_products = Vec::new();
    let mut map: HashMap::<_, Vec<(usize, bool)>> = HashMap::new();
    let mut t = 0;
    let c_in = c_in.decompose_tof();
    for (gate, q) in &c_in.circ {
        match &gate[..] {
            "h" => { tab.prepend_h(q[0]); },
            "x" => { tab.prepend_x(q[0]); },
            "z" => { tab.prepend_z(q[0]); },
            "s" => { tab.prepend_s(q[0]); tab.prepend_z(q[0]); },
            "cx" => { tab.prepend_cx(q.to_vec()); },
            "t" => { 
                let p = tab.stabs[q[0]].clone();
                let vec = p.get_boolean_vec(nb_qubits);
                let mut merge = map.contains_key(&vec);
                let mut value = Vec::new();
                if merge {
                    value = map.remove(&vec).unwrap();
                    let (index, sign) = value.pop().unwrap();
                    for i in (index+1)..t {
                        if v[i] && !p.is_commuting(&pauli_products[i]) {
                            merge = false;
                            break;
                        }
                    }
                    if merge {
                        r[index] = 0;
                        r[t] = 0;
                        if sign == p.sign { r[t] = 2; tab.prepend_s(q[0]); }
                    }
                }
                if !merge {
                    value.push((t, p.sign));
                    map.insert(vec, value);
                }
                pauli_products.push(p.clone());
                t += 1;
            },
            _ => {println!("Operator not implemented: {}", gate); std::process::exit(1)},
        }
    }
    let mut c = Circuit::new(nb_qubits);
    r.reverse();
    for (gate, q) in &c_in.circ {
        if gate == "t" {
            let val = r.pop().unwrap();
            if val == 1 { c.circ.push((gate.to_string(), q.to_vec())); }
            else if val == 2 { c.circ.push(("s".to_string(), q.to_vec())); }
        }
        else {
            c.circ.push((gate.to_string(), q.to_vec()));
        }
    }
    c
}

pub fn fast_t_merge(c_in: Circuit) -> Circuit {
    let nb_qubits = c_in.nb_qubits;
    let v = rank_vector(&c_in);
    let mut w = v.clone();
    let mut r = vec![1; v.len()];
    let mut tab = TableauColumnMajor::new(nb_qubits);
    let mut pauli_products = Vec::new();
    let mut map: HashMap::<_, Vec<(usize, bool)>> = HashMap::new();
    let mut t = 0;
    let c_in = c_in.decompose_tof();
    for (gate, q) in &c_in.circ {
        match &gate[..] {
            "h" => { tab.prepend_h(q[0]); },
            "x" => { tab.prepend_x(q[0]); },
            "z" => { tab.prepend_z(q[0]); },
            "s" => { tab.prepend_s(q[0]); tab.prepend_z(q[0]); },
            "cx" => { tab.prepend_cx(q.to_vec()); },
            "t" => { 
                let p = tab.stabs[q[0]].clone();
                let vec = p.get_boolean_vec(nb_qubits);
                let mut merge = map.contains_key(&vec);
                let mut value = Vec::new();
                if merge {
                    value = map.remove(&vec).unwrap();
                    let (index, sign) = value.pop().unwrap();
                    for i in (index+1)..t {
                        if v[i] && !p.is_commuting(&pauli_products[i]) {
                            if r[i] == 1 {
                                merge = false;
                            }
                            else {
                                for j in (i+1)..t {
                                    if w[j] && r[j] == 1 && !p.is_commuting(&pauli_products[j]) {
                                        merge = false;
                                        break;
                                    }
                                }
                            }
                            break;
                        }
                    }
                    if merge {
                        if v[index] {
                            for i in (index+1)..t {
                                w[i] = true;
                            }
                        }
                        w[index] = false;
                        r[index] = 0;
                        r[t] = 0;
                        if sign == p.sign { r[t] = 2; tab.prepend_s(q[0]); }
                    }
                }
                if !merge {
                    value.push((t, p.sign));
                    map.insert(vec, value);
                }
                pauli_products.push(p.clone());
                t += 1;
            },
            _ => {println!("Operator not implemented: {}", gate); std::process::exit(1)},
        }
    }
    let mut c = Circuit::new(c_in.nb_qubits);
    r.reverse();
    for (gate, q) in &c_in.circ {
        if gate == "t" {
            let val = r.pop().unwrap();
            if val == 1 { c.circ.push((gate.to_string(), q.to_vec())); }
            else if val == 2 { c.circ.push(("s".to_string(), q.to_vec())); }
        }
        else {
            c.circ.push((gate.to_string(), q.to_vec()));
        }
    }
    c
}

 fn diagonalize_pauli_rotation(tab: &mut Tableau, col: usize) -> bool {
    if let Some(pivot) = tab.x.iter().position(|x| x.get(col)) {
        for j in 0..tab.nb_qubits {
            if tab.x[j].get(col) && j != pivot {
                tab.append_cx(vec![pivot, j]);
            }
        }
        if tab.z[pivot].get(col) {
            tab.append_s(pivot);
        }
        tab.append_h(pivot);
        return true;
    }
    false
}

 fn diagonalize_tof(tab: &mut Tableau, cols: Vec::<usize>, h_gate: bool) -> Vec::<bool> {
    let mut vec = Vec::new();
    vec.push(diagonalize_pauli_rotation(tab, cols[0]));
    vec.push(diagonalize_pauli_rotation(tab, cols[1]));
    vec.push(diagonalize_pauli_rotation(tab, cols[2] + tab.nb_qubits * (h_gate as usize)));
    for _ in 0..4 {
        vec.push(false);
    }
    vec
}

 fn reverse_diagonalization(c_in: &Circuit) -> Tableau {
    let mut tab = Tableau::new(c_in.nb_qubits);
    for (gate, q) in &c_in.circ {
        match &gate[..] {
            "h" => { tab.prepend_h(q[0]); },
            "x" => { tab.prepend_x(q[0]); },
            "z" => { tab.prepend_z(q[0]); },
            "s" => { tab.prepend_s(q[0]); tab.prepend_z(q[0]); },
            "cx" => { tab.prepend_cx(q.to_vec()); },
            "t" | "ccz" | "tof" => continue,
            _ => {println!("Operator not implemented: {}", gate); std::process::exit(1)},
        }
    }
    for (gate, q) in c_in.circ.iter().rev() {
        match &gate[..] {
            "h" => { tab.prepend_h(q[0]); },
            "x" => { tab.prepend_x(q[0]); },
            "z" => { tab.prepend_z(q[0]); },
            "s" => { tab.prepend_s(q[0]); },
            "cx" => { tab.prepend_cx(q.to_vec()); },
            "t" => { diagonalize_pauli_rotation(&mut tab, q[0]); },
            "tof" => { diagonalize_tof(&mut tab, q.to_vec(), true); },
            "ccz" => { diagonalize_tof(&mut tab, q.to_vec(), false); },
            _ => {println!("Operator not implemented: {}", gate); std::process::exit(1)},
        }
    }
    tab
}

pub fn rank_vector(c_in: &Circuit) -> Vec::<bool> {
    let mut tab = reverse_diagonalization(c_in);
    let mut vec = Vec::new();
    for (gate, q) in &c_in.circ {
        match &gate[..] {
            "h" => { tab.prepend_h(q[0]); },
            "x" => { tab.prepend_x(q[0]); },
            "z" => { tab.prepend_z(q[0]); },
            "s" => { tab.prepend_s(q[0]); tab.prepend_z(q[0]); },
            "cx" => { tab.prepend_cx(q.to_vec()); },
            "t" => { vec.push(diagonalize_pauli_rotation(&mut tab, q[0])); },
            "tof" => { vec.append(&mut diagonalize_tof(&mut tab, q.to_vec(), true)); },
            "ccz" => { vec.append(&mut diagonalize_tof(&mut tab, q.to_vec(), false)); },
            _ => {println!("Operator not implemented: {}", gate); std::process::exit(1)},
        }
    }
    vec
}
