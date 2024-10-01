use crate::bit_vector::BitVector;
use crate::circuit::Circuit;
use crate::tableau::Tableau;

#[derive(Debug, Clone)]
pub struct PhasePolynomial {
    pub nb_qubits: usize,
    pub table: Vec<BitVector>,
}

impl PhasePolynomial {
    pub fn new(nb_qubits: usize) -> Self {
        PhasePolynomial {
            nb_qubits: nb_qubits,
            table: Vec::new(),
        }
    }

    pub fn clifford_correction(&self, table: &Vec<BitVector>, nb_qubits: usize) -> Tableau {
        let mut tab = Tableau::new(nb_qubits);
        for i in 0..nb_qubits {
            for j in (i+1)..nb_qubits {
                let z1 = (0..table.len()).filter(|&k| table[k].get(i) & table[k].get(j)).count();
                let z2 = (0..self.table.len()).filter(|&k| self.table[k].get(i) & self.table[k].get(j)).count();
                for _ in 0..((((z1 - z2) % 8 + 8) % 8) / 2) {
                    tab.append_cz(vec![i, j]);
                }
            }
            let z1 = (0..table.len()).filter(|&k| table[k].get(i)).count();
            let z2 = (0..self.table.len()).filter(|&k| self.table[k].get(i)).count();
            for _ in 0..((((z1 - z2) % 8 + 8) % 8) / 2) {
                tab.append_s(i);
            }
        }
        tab
    }

    pub fn to_circ(&self) -> Circuit {
        let mut c = Circuit::new(self.nb_qubits);
        for z in &self.table {
            let mut cnot_circ = Circuit::new(self.nb_qubits);
            let pivot = z.get_first_one();
            let mut indices = z.get_all_ones(self.nb_qubits);
            indices.swap_remove(0);
            for j in indices {
                cnot_circ.circ.push(("cx".into(), vec![j, pivot]));
            }
            c.append(cnot_circ.clone().circ);
            c.circ.push(("t".into(), vec![pivot]));
            c.append(cnot_circ.circ);
        }
        c
    }
}
