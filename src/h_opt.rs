use crate::pauli_product::PauliProduct;
use crate::tableau::Tableau;
use crate::circuit::Circuit;

 fn implement_pauli_z_rotation_from_pauli_product(tab: &mut Tableau, p: &PauliProduct) -> Circuit {
    let mut c = Circuit::new(tab.nb_qubits);
    let mut cnot_circ = Circuit::new(tab.nb_qubits);
    let pivot = p.z.get_first_one();
    let mut indices = p.z.get_all_ones(tab.nb_qubits);
    indices.swap_remove(0);
    for j in indices {
        cnot_circ.circ.push(("cx".into(), vec![j, pivot]));
    }
    c.append(cnot_circ.clone().circ);
    c.circ.push(("t".into(), vec![pivot]));
    if p.sign {
        c.circ.push(("s".into(), vec![pivot]));
        c.circ.push(("z".into(), vec![pivot]));
    }
    c.append(cnot_circ.circ);
    c
}

 fn implement_pauli_z_rotation(tab: &mut Tableau, col: usize) -> Circuit {
    let pivot = tab.z.iter().position(|z| z.get(col)).unwrap();
    let mut c = Circuit::new(tab.nb_qubits);
    let mut cnot_circ = Circuit::new(tab.nb_qubits);
    for j in 0..tab.nb_qubits {
        if tab.z[j].get(col) && j != pivot {
            cnot_circ.circ.push(("cx".into(), vec![j, pivot]));
        }
    }
    c.append(cnot_circ.clone().circ);
    c.circ.push(("t".into(), vec![pivot]));
    if tab.signs.get(col) {
        c.circ.push(("s".into(), vec![pivot]));
        c.circ.push(("z".into(), vec![pivot]));
    }
    c.append(cnot_circ.circ);
    c
}

 fn implement_pauli_rotation(tab: &mut Tableau, col: usize) -> Circuit {
    let mut c = Circuit::new(tab.nb_qubits);
    if let Some(pivot) = tab.x.iter().position(|x| x.get(col)) {
        for j in 0..tab.nb_qubits {
            if tab.x[j].get(col) && j != pivot {
                tab.append_cx(vec![pivot, j]);
                c.circ.push(("cx".into(), vec![pivot, j]));
            }
        }
        if tab.z[pivot].get(col) {
            tab.append_s(pivot);
            c.circ.push(("s".into(), vec![pivot]));
        }
        tab.append_h(pivot);
        c.circ.push(("h".into(), vec![pivot]));
    }
    c.append(implement_pauli_z_rotation(tab, col).circ);
    c
}


 fn implement_tof(tab: &mut Tableau, cols: Vec::<usize>, h_gate: bool) -> Circuit {
    let mut c = Circuit::new(tab.nb_qubits);
    c.append(implement_pauli_rotation(tab, cols[0]).circ);
    c.append(implement_pauli_rotation(tab, cols[1]).circ);
    c.append(implement_pauli_rotation(tab, cols[2] + tab.nb_qubits * (h_gate as usize)).circ);
    let mut p0 = tab.extract_pauli_product(cols[0]);
    let mut p1 = tab.extract_pauli_product(cols[1]);
    let p2 = tab.extract_pauli_product(cols[2] + tab.nb_qubits * (h_gate as usize));
    p0.z.xor(&p1.z);
    p0.sign ^= p1.sign ^ true;
    c.append(implement_pauli_z_rotation_from_pauli_product(tab, &p0).circ);
    p0.z.xor(&p2.z);
    p0.sign ^= p2.sign ^ true;
    c.append(implement_pauli_z_rotation_from_pauli_product(tab, &p0).circ);
    p0.z.xor(&p1.z);
    p0.sign ^= p1.sign ^ true;
    c.append(implement_pauli_z_rotation_from_pauli_product(tab, &p0).circ);
    p1.z.xor(&p2.z);
    p1.sign ^= p2.sign ^ true;
    c.append(implement_pauli_z_rotation_from_pauli_product(tab, &p1).circ);
    c
}

 fn h_opt_reverse(c_in: &Circuit) -> Tableau{
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
            "t" => { implement_pauli_rotation(&mut tab, q[0]); },
            "tof" => { implement_tof(&mut tab, q.to_vec(), true); },
            "ccz" => { implement_tof(&mut tab, q.to_vec(), false); },
            _ => {println!("Operator not implemented: {}", gate); std::process::exit(1)},
        }
    }
    tab
}

pub fn internal_h_opt(c_in: &Circuit) -> Circuit {
    let mut tab = h_opt_reverse(c_in);
    let mut c = tab.to_circ(false);
    for (gate, q) in &c_in.circ {
        match &gate[..] {
            "h" => { tab.prepend_h(q[0]); },
            "x" => { tab.prepend_x(q[0]); },
            "z" => { tab.prepend_z(q[0]); },
            "s" => { tab.prepend_s(q[0]); tab.prepend_z(q[0]); },
            "cx" => { tab.prepend_cx(q.to_vec()); },
            "t" => { c.append(implement_pauli_rotation(&mut tab, q[0]).circ); },
            "tof" => { c.append(implement_tof(&mut tab, q.to_vec(), true).circ); },
            "ccz" => { c.append(implement_tof(&mut tab, q.to_vec(), false).circ); },
            _ => {println!("Operator not implemented: {}", gate); std::process::exit(1)},
        }
    }
    c.append(tab.to_circ(true).circ);
    c
}
