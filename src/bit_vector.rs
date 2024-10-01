// A wrapper type to enforce 32-byte alignment for SIMD loads
#[repr(align(32))]
struct BitLanes([i32; 8]);

#[derive(Debug, Clone, Copy)]
pub struct BitBlock {
    #[cfg(target_feature = "avx2")]
    inner: std::arch::x86_64::__m256i,

    #[cfg(not(target_feature = "avx2"))]
    inner: [i32; 8]
}

impl BitBlock {
    #[cfg(target_feature = "avx2")]
    fn constant(a: i32) -> Self {
        BitBlock {
            inner: unsafe { std::arch::x86_64::_mm256_set1_epi32(a) }
        }
    }

    #[cfg(not(target_feature = "avx2"))]
    fn constant(a: i32) -> Self {
        BitBlock {
            inner: [a; 8]
        }
    }

    #[cfg(target_feature = "avx2")]
    fn load(arr: &BitLanes) -> Self {
        BitBlock {
            // SAFETY: BitLanes is 32-byte aligned
            inner: unsafe { std::arch::x86_64::_mm256_load_si256(arr.0.as_ptr() as *const _) }
        }
    }

    #[cfg(not(target_feature = "avx2"))]
    fn load(arr: &BitLanes) -> Self {
        BitBlock {
            inner: arr.0.clone()
        }
    }
    
    #[cfg(target_feature = "avx2")]
    fn zero() -> Self {
        BitBlock {
            inner: unsafe { std::arch::x86_64::_mm256_setzero_si256() }
        }
    }

    #[cfg(not(target_feature = "avx2"))]
    fn zero() -> Self {
        Self::constant(0)
    }

    #[cfg(target_feature = "avx2")]
    fn extract(&self) -> [i32; 8] {
        let mut arr = BitLanes([0; 8]);
        // SAFETY: BitLanes is 32-byte aligned
        unsafe { std::arch::x86_64::_mm256_store_si256(arr.0.as_mut_ptr() as *mut _, self.inner); }
        arr.0
    }

    #[cfg(not(target_feature = "avx2"))]
    fn extract(&self) -> [i32; 8] {
        self.inner
    }
}

impl std::ops::BitXorAssign for BitBlock {
    #[cfg(target_feature = "avx2")]
    fn bitxor_assign(&mut self, rhs: Self) {
        self.inner = unsafe { std::arch::x86_64::_mm256_xor_si256(self.inner, rhs.inner) };
    }

    #[cfg(not(target_feature = "avx2"))]
    fn bitxor_assign(&mut self, rhs: Self) {
        for i in 0..8 {
            self.inner[i] ^= rhs.inner[i];
        }
    }
}

impl std::ops::BitAndAssign for BitBlock {
    #[cfg(target_feature = "avx2")]
    fn bitand_assign(&mut self, rhs: Self) {
        self.inner = unsafe { std::arch::x86_64::_mm256_and_si256(self.inner, rhs.inner) };
    }

    #[cfg(not(target_feature = "avx2"))]
    fn bitand_assign(&mut self, rhs: Self) {
        for i in 0..8 {
            self.inner[i] &= rhs.inner[i];
        }
    }
}

#[derive(Debug, Clone)]
pub struct BitVector {
    pub blocks: Vec<BitBlock>,
}

impl BitVector {
    const LANES: usize = 8;
    const LANE_SIZE: usize = 32;
    const BLOCK_SIZE: usize = 256;

    pub fn new(nb_bits: usize) -> Self {
        BitVector {
            blocks: BitVector::init_blocks(nb_bits),
        }
    }

    pub fn new_block_size(nb_blocks: usize) -> Self {
        BitVector {
            blocks: BitVector::init_blocks(nb_blocks * BitVector::BLOCK_SIZE - 1),
        }
    }

    pub fn from_integer_vec(vec: Vec<i128>) -> Self {
        let mut bv = BitVector::new(vec.len() * 128 - 1);
        let mut arr = BitLanes([0; BitVector::LANES]);
        let mut block_index = 0;
        let mut index = 0;
        for v in vec {
            let mut val = v.clone();
            for i in 0..4 {
                arr.0[index] = val as u32 as i32;
                index += 1;
                if i < 3 {
                    val = val >> 32;
                }
            }
            if index == 8 {
                index = 0;
                bv.blocks[block_index] = BitBlock::load(&arr);
                block_index += 1;
                arr = BitLanes([0; BitVector::LANES]);
            }
        }
        if index > 0 {
            bv.blocks[block_index] = BitBlock::load(&arr);
        }
        bv
    }

    fn init_blocks(nb_bits: usize) -> Vec<BitBlock> {
        let mut vec: Vec<BitBlock> = Vec::with_capacity(nb_bits / BitVector::BLOCK_SIZE + 1);
        for _ in 0..vec.capacity() {
            vec.push(BitBlock::zero());
        }
        vec
    }

    pub fn size(&self) -> usize {
        self.blocks.len() * BitVector::BLOCK_SIZE
    }

    pub fn xor_bit(&mut self, mut bit: usize) {
        let block_index = bit / BitVector::BLOCK_SIZE;
        bit = bit % BitVector::BLOCK_SIZE;
        let lane_index = bit / BitVector::LANE_SIZE;
        bit = bit % BitVector::LANE_SIZE;
        let mut arr = BitLanes([0; BitVector::LANES]);
        arr.0[lane_index] ^= 1 << bit;
        self.blocks[block_index] ^= BitBlock::load(&arr);
    }

    pub fn get(&self, mut bit: usize) -> bool {
        let block_index = bit / BitVector::BLOCK_SIZE;
        bit = bit % BitVector::BLOCK_SIZE;
        let lane_index = bit / 32;
        bit = bit % 32;
        self.extract_block(block_index)[lane_index] & (1 << bit) != 0
    }

    pub fn get_first_one(&self) -> usize {
        for i in 0..self.blocks.len() {
            let block = self.extract_block(i);
            for j in 0..8 {
                for k in 0..32 {
                    if block[j] & (1 << k) != 0 { return i * BitVector::BLOCK_SIZE + j * BitVector::LANE_SIZE + k; }
                }
            }
        }
        0
    }

    pub fn get_all_ones(&self, nb_bits: usize) -> Vec<usize> {
        let mut index = 0;
        let mut vec = Vec::new();
        for i in 0..self.blocks.len() {
            let block = self.extract_block(i);
            for j in 0..8 {
                for k in 0..32 {
                    if block[j] & (1 << k) != 0 { vec.push(index); }
                    index += 1;
                    if index >= nb_bits { return vec; }
                }
            }
        }
        vec
    }

    pub fn xor(&mut self, bv: &BitVector) {
        for i in 0..self.blocks.len() {
            self.blocks[i] ^= bv.blocks[i];
        }
    }

    pub fn and(&mut self, bv: &BitVector) {
        for i in 0..self.blocks.len() {
            self.blocks[i] &= bv.blocks[i];
        }
    }

    pub fn negate(&mut self) {
        let a: i32 = !0;
        for i in 0..self.blocks.len() {
            self.blocks[i] ^= BitBlock::constant(a);
        }
    }
    
    pub fn extend_vec(&mut self, vec: Vec<bool>, nb_bits: usize) {
        let mut bit = nb_bits;
        let mut block_index = bit / BitVector::BLOCK_SIZE;
        bit = bit % BitVector::BLOCK_SIZE;
        let mut lane_index = bit / BitVector::LANE_SIZE;
        bit = bit % BitVector::LANE_SIZE;
        let mut arr = BitLanes([0; BitVector::LANES]);

        for val in vec {
            if bit == BitVector::LANE_SIZE {
                bit = 0;
                lane_index += 1;
                if lane_index == BitVector::LANES {
                    lane_index = 0;
                    self.blocks[block_index] ^= BitBlock::load(&arr);
                    block_index += 1;
                    self.blocks.push(BitBlock::zero());
                    arr = BitLanes([0; BitVector::LANES]);
                }
            }
            if val {
                arr.0[lane_index] ^= 1 << bit;
            }
            bit += 1;
        }
        self.blocks[block_index] ^= BitBlock::load(&arr);
    }

    pub fn get_boolean_vec(&self) -> Vec<bool> {
        let mut vec: Vec<bool> = Vec::with_capacity(self.blocks.len() * BitVector::BLOCK_SIZE);
        for block_index in 0..self.blocks.len() {
            let arr = self.extract_block(block_index);
            for j in 0..8 {
                for i in 0..32 {
                    vec.push(arr[j] & (1 << i) != 0);
                }
            }
        }
        vec
    }

    pub fn get_integer_vec(&self) -> Vec<i128> {
        let mut vec: Vec<i128> = Vec::with_capacity(self.blocks.len() * 2);
        for block_index in 0..self.blocks.len() {
            let arr = self.extract_block(block_index);
            for k in 0..2 {
                let mut integer: i128 = 0;
                for j in 0..4 {
                    integer ^= (arr[k*4 + j] as u32 as i128) << (32 * j);
                }
                vec.push(integer);
            }
        }
        vec
    }

    pub fn popcount(&self) -> i32 {
        let mut sum: i32 = 0;
        for block_index in 0..self.blocks.len() {
            let arr = self.extract_block(block_index);
            for j in 0..8 {
                sum += arr[j].count_ones() as i32;
            }
        }
        sum
    }

    fn extract_block(&self, block: usize) -> [i32; 8] {
        self.blocks[block].extract()
    }
}
