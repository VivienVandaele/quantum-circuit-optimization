use std::arch::x86_64::__m256i;
use std::arch::x86_64::{_mm256_setzero_si256, _mm256_load_epi32, _mm256_extract_epi32, _mm256_xor_si256 as xor_256, _popcnt32, _mm256_and_si256 as and_256, _mm256_set1_epi32};

#[derive(Debug, Clone)]
pub struct BitVector {
    pub blocks: Vec<__m256i>,
}

impl BitVector {
    const LANES: usize = 8;
    const LANE_SIZE: usize = 32;
    const BLOCK_SIZE: usize = 256;

    pub unsafe fn new(nb_bits: usize) -> Self {
        BitVector {
            blocks: BitVector::init_blocks(nb_bits),
        }
    }

    unsafe fn init_blocks(nb_bits: usize) -> Vec<__m256i> {
        let mut vec: Vec<__m256i> = Vec::with_capacity(nb_bits / BitVector::BLOCK_SIZE + 1);
        for _ in 0..vec.capacity() {
            vec.push(_mm256_setzero_si256());
        }
        vec
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn xor_bit(&mut self, mut bit: usize) {
        let block_index = bit / BitVector::BLOCK_SIZE;
        bit = bit % BitVector::BLOCK_SIZE;
        let lane_index = bit / BitVector::LANE_SIZE;
        bit = bit % BitVector::LANE_SIZE;
        let mut arr: [i32; BitVector::LANES] = [0; BitVector::LANES];
        arr[lane_index] ^= 1 << bit;
        self.blocks[block_index] = xor_256(_mm256_load_epi32(arr.as_ptr()), self.blocks[block_index]);
    }

    pub unsafe fn get(&self, mut bit: usize) -> bool {
        let block_index = bit / BitVector::BLOCK_SIZE;
        bit = bit % BitVector::BLOCK_SIZE;
        let lane_index = bit / 32;
        bit = bit % 32;
        self.extract_block(block_index)[lane_index] & (1 << bit) != 0
    }

    pub unsafe fn get_first_one(&self) -> usize {
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

    pub unsafe fn get_all_ones(&self, nb_bits: usize) -> Vec<usize> {
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

    #[target_feature(enable = "avx2")]
    pub unsafe fn xor(&mut self, bv: &BitVector) {
        for i in 0..self.blocks.len() {
            self.blocks[i] = xor_256(self.blocks[i], bv.blocks[i]);
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn and(&mut self, bv: &BitVector) {
        for i in 0..self.blocks.len() {
            self.blocks[i] = and_256(self.blocks[i], bv.blocks[i]);
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn negate(&mut self) {
        let a: i32 = !0;
        for i in 0..self.blocks.len() {
            self.blocks[i] = xor_256(self.blocks[i], _mm256_set1_epi32(a));
        }
    }

    pub unsafe fn get_boolean_vec(&self) -> Vec<bool> {
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

    #[target_feature(enable = "avx2")]
    pub unsafe fn popcount(&self) -> i32 {
        let mut sum: i32 = 0;
        for block_index in 0..self.blocks.len() {
            let arr = self.extract_block(block_index);
            for j in 0..8 {
                sum += _popcnt32(arr[j]);
            }
        }
        sum
    }

    #[target_feature(enable = "avx2")]
    unsafe fn extract_block(&self, index: usize) -> [i32; 8] {
        [_mm256_extract_epi32::<0>(self.blocks[index]), _mm256_extract_epi32::<1>(self.blocks[index]),
        _mm256_extract_epi32::<2>(self.blocks[index]), _mm256_extract_epi32::<3>(self.blocks[index]),
        _mm256_extract_epi32::<4>(self.blocks[index]), _mm256_extract_epi32::<5>(self.blocks[index]),
        _mm256_extract_epi32::<6>(self.blocks[index]), _mm256_extract_epi32::<7>(self.blocks[index])]
    }
}
