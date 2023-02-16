use crate::bit_vector::BitVector;
use crate::pauli_product::PauliProduct;
use crate::circuit::Circuit;

#[derive(Debug, Clone)]
pub struct Tableau {
    pub nb_qubits: usize,
    pub z: Vec<BitVector>,
    pub x: Vec<BitVector>,
    pub signs: BitVector,
}

impl Tableau {
    pub unsafe fn new(nb_qubits: usize) -> Self {
        Tableau {
            nb_qubits: nb_qubits,
            z: Tableau::init_z(nb_qubits),
            x: Tableau::init_x(nb_qubits),
            signs: BitVector::new(nb_qubits << 1),
        }
    }

    unsafe fn init_z(nb_qubits: usize) -> Vec<BitVector> {
        let mut vec = Vec::new();
        for i in 0..nb_qubits {
            let mut bv = BitVector::new(nb_qubits << 1);
            bv.xor_bit(i);
            vec.push(bv);
        }
        vec
    }

    unsafe fn init_x(nb_qubits: usize) -> Vec<BitVector> {
        let mut vec = Vec::new();
        for i in 0..nb_qubits {
            let mut bv = BitVector::new(nb_qubits << 1);
            bv.xor_bit(i + nb_qubits);
            vec.push(bv);
        }
        vec
    }

    pub unsafe fn append_x(&mut self, qubit: usize) {
        self.signs.xor(&self.z[qubit]);
    }

    pub unsafe fn append_z(&mut self, qubit: usize) {
        self.signs.xor(&self.x[qubit]);
    }

    pub unsafe fn append_v(&mut self, qubit: usize) {
        let mut a = self.x[qubit].clone();
        a.negate();
        a.and(&self.z[qubit]);
        self.signs.xor(&a);
        self.x[qubit].xor(&self.z[qubit]);
    }

    pub unsafe fn append_s(&mut self, qubit: usize) {
        let mut a = self.z[qubit].clone();
        a.and(&self.x[qubit]);
        self.signs.xor(&a);
        self.z[qubit].xor(&self.x[qubit]);
    }

    pub unsafe fn append_h(&mut self, qubit: usize) {
        self.append_s(qubit);
        self.append_v(qubit);
        self.append_s(qubit);
    }

    pub unsafe fn append_cx(&mut self, qubits: Vec<usize>) {
        let mut a =  self.z[qubits[0]].clone();
        a.negate();
        a.xor(&self.x[qubits[1]]);
        a.and(&self.z[qubits[1]]);
        a.and(&self.x[qubits[0]]);
        self.signs.xor(&a);
        let a = self.z[qubits[1]].clone();
        self.z[qubits[0]].xor(&a);
        let a = self.x[qubits[0]].clone();
        self.x[qubits[1]].xor(&a);
    }

    pub unsafe fn extract_pauli_product(&self, col: usize) -> PauliProduct {
        let mut z = BitVector::new(self.nb_qubits);
        let mut x = BitVector::new(self.nb_qubits);
        for i in 0..self.nb_qubits {
            if self.z[i].get(col) { z.xor_bit(i); }
            if self.x[i].get(col) { x.xor_bit(i); }
        }
        PauliProduct::new(z, x, self.signs.get(col))
    }

    pub unsafe fn insert_pauli_product(&mut self, p: PauliProduct, col: usize) {
        let p_x = p.x.get_boolean_vec();
        let p_z = p.z.get_boolean_vec();
        for i in 0..self.nb_qubits {
            if p_z[i] ^ self.z[i].get(col) {
                self.z[i].xor_bit(col);
            }
            if p_x[i] ^ self.x[i].get(col) {
                self.x[i].xor_bit(col);
            }
        }
        if p.sign ^ self.signs.get(col) {
            self.signs.xor_bit(col);
        }
    }

    pub unsafe fn prepend_x(&mut self, qubit: usize) {
        self.signs.xor_bit(qubit);
    }

    pub unsafe fn prepend_z(&mut self, qubit: usize) {
        self.signs.xor_bit(qubit + self.nb_qubits);
    }

    pub unsafe fn prepend_s(&mut self, qubit: usize) {
        let stab = self.extract_pauli_product(qubit);
        let mut destab = self.extract_pauli_product(qubit + self.nb_qubits);
        destab.pauli_product_mult(&stab);
        self.insert_pauli_product(destab, qubit + self.nb_qubits);
    }

    pub unsafe fn prepend_h(&mut self, qubit: usize) {
        let stab = self.extract_pauli_product(qubit);
        let destab = self.extract_pauli_product(qubit + self.nb_qubits);
        self.insert_pauli_product(destab, qubit);
        self.insert_pauli_product(stab, qubit + self.nb_qubits);
    }

    pub unsafe fn prepend_cx(&mut self, qubits: Vec<usize>) {
        let stab_ctrl = self.extract_pauli_product(qubits[0]);
        let mut stab_targ = self.extract_pauli_product(qubits[1]);
        let mut destab_ctrl = self.extract_pauli_product(qubits[0] + self.nb_qubits);
        let destab_targ = self.extract_pauli_product(qubits[1] + self.nb_qubits);
        stab_targ.pauli_product_mult(&stab_ctrl);
        destab_ctrl.pauli_product_mult(&destab_targ);
        self.insert_pauli_product(stab_targ, qubits[1]);
        self.insert_pauli_product(destab_ctrl, qubits[0] + self.nb_qubits);
    }

    pub unsafe fn to_circ(&self, inverse: bool) -> Circuit {
        let mut tab = self.clone();
        let mut c = Circuit::new(self.nb_qubits);
        for i in 0..self.nb_qubits {
            if let Some(index) = tab.x.iter().position(|x| x.get(i)) {
                for j in (i+1)..self.nb_qubits {
                    if tab.x[j].get(i) && j != index {
                        tab.append_cx(vec![index, j]);
                        c.circ.push(("cx".into(), vec![index, j]));
                    }
                }
                if tab.z[index].get(i) {
                    tab.append_s(index);
                    c.circ.push(("s".into(), vec![index]));
                }
                tab.append_h(index);
                c.circ.push(("h".into(), vec![index]));
            }
            if !tab.z[i].get(i) {
                let index = tab.z.iter().position(|z| z.get(i)).unwrap();
                tab.append_cx(vec![i, index]);
                c.circ.push(("cx".into(), vec![i, index]));
            }
            for j in 0..self.nb_qubits {
                if tab.z[j].get(i) && j != i {
                    tab.append_cx(vec![j, i]);
                    c.circ.push(("cx".into(), vec![j, i]));
                }
            }
            for j in 0..self.nb_qubits {
                if tab.x[j].get(i + self.nb_qubits) && j != i {
                    tab.append_cx(vec![i, j]);
                    c.circ.push(("cx".into(), vec![i, j]));
                }
            }
            for j in 0..self.nb_qubits {
                if tab.z[j].get(i + self.nb_qubits) && j != i {
                    tab.append_cx(vec![i, j]);
                    c.circ.push(("cx".into(), vec![i, j]));
                    tab.append_s(j);
                    c.circ.push(("s".into(), vec![j]));
                    tab.append_cx(vec![i, j]);
                    c.circ.push(("cx".into(), vec![i, j]));
                }
            }
            if tab.z[i].get(i + self.nb_qubits) {
                tab.append_s(i);
                c.circ.push(("s".into(), vec![i]));
            }
            if tab.signs.get(i) {
                tab.append_x(i);
                c.circ.push(("x".into(), vec![i]));
            }
            if tab.signs.get(i + self.nb_qubits) {
                tab.append_z(i);
                c.circ.push(("z".into(), vec![i]));
            }
        }
        if !inverse {
            let mut c2 = Circuit::new(self.nb_qubits);
            for (gate, qubits) in c.circ.into_iter().rev() {
                c2.circ.push((gate.to_string(), qubits.to_vec()));
                if gate == "s" { c2.circ.push(("z".into(), qubits.to_vec())); }
            }
            return c2;
        }
        c
    }
}
