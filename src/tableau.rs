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
    pub fn new(nb_qubits: usize) -> Self {
        Tableau {
            nb_qubits: nb_qubits,
            z: Tableau::init_z(nb_qubits),
            x: Tableau::init_x(nb_qubits),
            signs: BitVector::new(nb_qubits << 1),
        }
    }

     fn init_z(nb_qubits: usize) -> Vec<BitVector> {
        let mut vec = Vec::new();
        for i in 0..nb_qubits {
            let mut bv = BitVector::new(nb_qubits << 1);
            bv.xor_bit(i);
            vec.push(bv);
        }
        vec
    }

     fn init_x(nb_qubits: usize) -> Vec<BitVector> {
        let mut vec = Vec::new();
        for i in 0..nb_qubits {
            let mut bv = BitVector::new(nb_qubits << 1);
            bv.xor_bit(i + nb_qubits);
            vec.push(bv);
        }
        vec
    }

    pub fn append_x(&mut self, qubit: usize) {
        self.signs.xor(&self.z[qubit]);
    }

    pub fn append_z(&mut self, qubit: usize) {
        self.signs.xor(&self.x[qubit]);
    }

    pub fn append_v(&mut self, qubit: usize) {
        let mut a = self.x[qubit].clone();
        a.negate();
        a.and(&self.z[qubit]);
        self.signs.xor(&a);
        self.x[qubit].xor(&self.z[qubit]);
    }

    pub fn append_s(&mut self, qubit: usize) {
        let mut a = self.z[qubit].clone();
        a.and(&self.x[qubit]);
        self.signs.xor(&a);
        self.z[qubit].xor(&self.x[qubit]);
    }

    pub fn append_h(&mut self, qubit: usize) {
        self.append_s(qubit);
        self.append_v(qubit);
        self.append_s(qubit);
    }

    pub fn append_cx(&mut self, qubits: Vec<usize>) {
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

    pub fn append_cz(&mut self, qubits: Vec<usize>) {
        self.append_s(qubits[0]);
        self.append_s(qubits[1]);
        self.append_cx(qubits.to_vec());
        self.append_s(qubits[1]);
        self.append_z(qubits[1]);
        self.append_cx(qubits);
    }

    pub fn extract_pauli_product(&self, col: usize) -> PauliProduct {
        let mut z = BitVector::new(self.nb_qubits);
        let mut x = BitVector::new(self.nb_qubits);
        for i in 0..self.nb_qubits {
            if self.z[i].get(col) { z.xor_bit(i); }
            if self.x[i].get(col) { x.xor_bit(i); }
        }
        PauliProduct::new(z, x, self.signs.get(col))
    }

    pub fn insert_pauli_product(&mut self, p: PauliProduct, col: usize) {
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

    pub fn prepend_x(&mut self, qubit: usize) {
        self.signs.xor_bit(qubit);
    }

    pub fn prepend_z(&mut self, qubit: usize) {
        self.signs.xor_bit(qubit + self.nb_qubits);
    }

    pub fn prepend_s(&mut self, qubit: usize) {
        let stab = self.extract_pauli_product(qubit);
        let mut destab = self.extract_pauli_product(qubit + self.nb_qubits);
        destab.pauli_product_mult(&stab);
        self.insert_pauli_product(destab, qubit + self.nb_qubits);
    }

    pub fn prepend_h(&mut self, qubit: usize) {
        let stab = self.extract_pauli_product(qubit);
        let destab = self.extract_pauli_product(qubit + self.nb_qubits);
        self.insert_pauli_product(destab, qubit);
        self.insert_pauli_product(stab, qubit + self.nb_qubits);
    }

    pub fn prepend_cx(&mut self, qubits: Vec<usize>) {
        let stab_ctrl = self.extract_pauli_product(qubits[0]);
        let mut stab_targ = self.extract_pauli_product(qubits[1]);
        let mut destab_ctrl = self.extract_pauli_product(qubits[0] + self.nb_qubits);
        let destab_targ = self.extract_pauli_product(qubits[1] + self.nb_qubits);
        stab_targ.pauli_product_mult(&stab_ctrl);
        destab_ctrl.pauli_product_mult(&destab_targ);
        self.insert_pauli_product(stab_targ, qubits[1]);
        self.insert_pauli_product(destab_ctrl, qubits[0] + self.nb_qubits);
    }

    pub fn to_circ(&self, inverse: bool) -> Circuit {
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

#[derive(Debug, Clone)]
pub struct TableauColumnMajor {
    pub nb_qubits: usize,
    pub stabs: Vec<PauliProduct>,
    pub destabs: Vec<PauliProduct>,
}

impl TableauColumnMajor {
    pub fn new(nb_qubits: usize) -> Self {
        TableauColumnMajor {
            nb_qubits: nb_qubits,
            stabs: TableauColumnMajor::init_stabs(nb_qubits),
            destabs: TableauColumnMajor::init_destabs(nb_qubits),
        }
    }

     fn init_stabs(nb_qubits: usize) -> Vec<PauliProduct> {
        let mut vec = Vec::new();
        for i in 0..nb_qubits {
            let mut bv = BitVector::new(nb_qubits);
            bv.xor_bit(i);
            vec.push(PauliProduct::new(bv, BitVector::new(nb_qubits), false));
        }
        vec
    }

     fn init_destabs(nb_qubits: usize) -> Vec<PauliProduct> {
        let mut vec = Vec::new();
        for i in 0..nb_qubits {
            let mut bv = BitVector::new(nb_qubits);
            bv.xor_bit(i);
            vec.push(PauliProduct::new(BitVector::new(nb_qubits), bv, false));
        }
        vec
    }

    pub fn prepend_x(&mut self, qubit: usize) {
        self.stabs[qubit].sign ^= true;
    }

    pub fn prepend_z(&mut self, qubit: usize) {
        self.destabs[qubit].sign ^= true;
    }

    pub fn prepend_v(&mut self, qubit: usize) {
        self.stabs[qubit].pauli_product_mult(&self.destabs[qubit]);
    }

    pub fn prepend_s(&mut self, qubit: usize) {
        self.destabs[qubit].pauli_product_mult(&self.stabs[qubit]);
    }

    pub fn prepend_h(&mut self, qubit: usize) {
        self.prepend_s(qubit);
        self.prepend_v(qubit);
        self.prepend_s(qubit);
    }

    pub fn prepend_cx(&mut self, qubits: Vec<usize>) {
        let p = self.stabs[qubits[0]].clone();
        self.stabs[qubits[1]].pauli_product_mult(&p);
        let p = self.destabs[qubits[1]].clone();
        self.destabs[qubits[0]].pauli_product_mult(&p);
    }

    pub fn to_circ(&self, inverse: bool) -> Circuit {
        let mut tab = self.clone();
        let mut c = Circuit::new(tab.nb_qubits);
        for i in 0..tab.nb_qubits {
            if let Some(index) = tab.stabs.iter().position(|p| p.x.get(i) ) {
                for j in (i+1)..tab.nb_qubits {
                    if tab.stabs[j].x.get(i) && j != index {
                        tab.prepend_cx(vec![index, j]);
                        c.circ.push(("cx".into(), vec![index, j]));
                    }
                }
                if tab.destabs[index].x.get(i) {
                    tab.prepend_s(index);
                    c.circ.push(("s".into(), vec![index]));
                }
                tab.prepend_h(index);
                c.circ.push(("h".into(), vec![index]));
            }
            if !tab.destabs[i].x.get(i) {
                let index = tab.destabs.iter().position(|p| p.x.get(i)).unwrap();
                tab.prepend_cx(vec![i, index]);
                c.circ.push(("cx".into(), vec![i, index]));
            }
            for j in 0..tab.nb_qubits {
                if tab.destabs[j].x.get(i) && j != i {
                    tab.prepend_cx(vec![j, i]);
                    c.circ.push(("cx".into(), vec![j, i]));
                }
            }
            for j in 0..tab.nb_qubits {
                if tab.stabs[j].z.get(i) && j != i {
                    tab.prepend_cx(vec![i, j]);
                    c.circ.push(("cx".into(), vec![i, j]));
                }
            }
            for j in 0..tab.nb_qubits {
                if tab.destabs[j].z.get(i) && j != i {
                    tab.prepend_cx(vec![i, j]);
                    c.circ.push(("cx".into(), vec![i, j]));
                    tab.prepend_s(j);
                    c.circ.push(("s".into(), vec![j]));
                    tab.prepend_cx(vec![i, j]);
                    c.circ.push(("cx".into(), vec![i, j]));
                }
            }
            if tab.destabs[i].z.get(i) {
                tab.prepend_s(i);
                c.circ.push(("s".into(), vec![i]));
            }
            if tab.stabs[i].sign {
                tab.prepend_x(i);
                c.circ.push(("x".into(), vec![i]));
            }
            if tab.destabs[i].sign {
                tab.prepend_z(i);
                c.circ.push(("z".into(), vec![i]));
            }
        }
        c.circ.reverse();
        if !inverse {
            let mut c2 = Circuit::new(tab.nb_qubits);
            for (gate, qubits) in c.circ.into_iter().rev() {
                c2.circ.push((gate.to_string(), qubits.to_vec()));
                if gate == "s" { c2.circ.push(("z".into(), qubits.to_vec())); }
            }
            return c2;
        }
        c
    }
}
