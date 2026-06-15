use crate::bit_vector::BitVector;

#[derive(Debug, Clone)]
pub struct PauliProduct {
    pub z: BitVector,
    pub x: BitVector,
    pub sign: bool,
}

impl PauliProduct {
    pub fn new(z: BitVector, x: BitVector, sign: bool) -> Self {
        PauliProduct {
            z: z,
            x: x,
            sign: sign,
        }
    }

    pub fn is_commuting(&self, p: &PauliProduct) -> bool {
        let mut x1z2 = self.z.clone();
        x1z2.and(&p.x);
        let mut ac = self.x.clone();
        ac.and(&p.z);
        ac.xor(&x1z2);
        ac.popcount() % 2 == 0
    }

    pub fn get_integer_vec(&self) -> Vec<i128> {
        let mut vec = self.z.get_integer_vec();
        vec.append(&mut self.x.get_integer_vec());
        vec
    }

    pub fn pauli_product_mult(&mut self, p: &PauliProduct) {
        let mut x1z2 = self.z.clone();
        x1z2.and(&p.x);
        let mut ac = self.x.clone();
        ac.and(&p.z);
        ac.xor(&x1z2);
        self.x.xor(&p.x);
        self.z.xor(&p.z);
        x1z2.xor(&self.x);
        x1z2.xor(&self.z);
        x1z2.and(&ac);
        self.sign ^= p.sign ^ (((ac.popcount() + 2*x1z2.popcount()) % 4) > 1);
    }
}
