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

    pub unsafe fn pauli_product_mult(&mut self, p: &PauliProduct) {
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
