use crate::bit_vector::BitVector;
use hashbrown::HashMap;

pub fn proper(mut table: Vec<BitVector>) -> Vec<BitVector> {
    let mut map = HashMap::new();
    let mut to_remove: Vec<usize> = Vec::new();
    for i in 0..table.len() {
        let col = table[i].get_boolean_vec();
        if !table[i].get(table[i].get_first_one()) {
            to_remove.push(i);
        }
        else if map.contains_key(&col) {
            to_remove.push(*map.get(&col).unwrap());
            to_remove.push(i);
            map.remove(&col);
        }
        else {
            map.insert(col, i);
        }
    }
    to_remove.sort_by(|a, b| b.cmp(a));
    for i in to_remove {
        table.swap_remove(i);
    }
    table
}

pub fn to_remove(table: &Vec<BitVector>) -> Vec<usize> {
    let mut map = HashMap::new();
    let mut to_remove: Vec<usize> = Vec::new();
    for i in 0..table.len() {
        let col = table[i].get_integer_vec();
        if !table[i].get(table[i].get_first_one()) {
            to_remove.push(i);
        }
        else if map.contains_key(&col) {
            to_remove.push(*map.get(&col).unwrap());
            to_remove.push(i);
            map.remove(&col);
        }
        else {
            map.insert(col, i);
        }
    }
    to_remove
}

pub fn kernel(matrix: &mut Vec<BitVector>,augmented_matrix: &mut Vec<BitVector>,
                               pivots: &mut HashMap<usize, usize>) -> Option<BitVector> {
    for i in 0..matrix.len() {
        if pivots.contains_key(&i) { continue; }
        for (key, value) in pivots.clone() {
            if matrix[i].get(value) {
                let pivot = matrix[key].clone();
                let augmented_pivot = augmented_matrix[key].clone();
                matrix[i].xor(&pivot);
                augmented_matrix[i].xor(&augmented_pivot);
            }
        }
        let index = matrix[i].get_first_one();
        if matrix[i].get(index) {
            let pivot = matrix[i].clone();
            let augmented_pivot = augmented_matrix[i].clone();
            for j in pivots.keys() {
                if matrix[*j].get(index) {
                    matrix[*j].xor(&pivot);
                    augmented_matrix[*j].xor(&augmented_pivot);
                }
            }
            pivots.insert(i, index);
        }
        else {
            return Some(augmented_matrix[i].clone());
        }
    }
    None
}

pub fn tohpe(mut table: Vec<BitVector>, nb_qubits: usize) -> Vec<BitVector> {
     fn clear_column(i: usize, matrix: &mut Vec<BitVector>, augmented_matrix: &mut Vec<BitVector>,
                           pivots: &mut HashMap<usize, usize>) {
        if !pivots.contains_key(&i) { return; }
        let val = pivots.remove(&i).unwrap();
        if !augmented_matrix[i].get(i) {
            for j in 0..matrix.len() {
                if !augmented_matrix[j].get(i) { continue; }
                pivots.insert(j, val);
                let col = matrix[j].clone();
                let augmented_col = augmented_matrix[j].clone();
                matrix[j] = matrix[i].clone();
                augmented_matrix[j] = augmented_matrix[i].clone();
                matrix[i] = col;
                augmented_matrix[i] = augmented_col;
                break;
            }
        }
        let col = matrix[i].clone();
        let augmented_col = augmented_matrix[i].clone();
        for j in 0..matrix.len() {
            if augmented_matrix[j].get(i) && i != j {
                matrix[j].xor(&col);
                augmented_matrix[j].xor(&augmented_col);
            }
        }
    }

    let mut matrix = table.clone();
    for i in 0..table.len() {
        let mut t_vec = table[i].get_boolean_vec();
        t_vec.truncate(nb_qubits);
        let mut vec = Vec::<bool>::new();
        for _ in 0..nb_qubits {
            if t_vec.pop().unwrap() {
                vec.append(&mut t_vec.clone());
            }
            else {
                vec.append(&mut vec![false; t_vec.len()]);
            }
        }
        matrix[i].extend_vec(vec, nb_qubits);
    }
    let mut pivots: HashMap::<usize, usize> = HashMap::new();
    let mut augmented_matrix: Vec<BitVector> = Vec::with_capacity(table.len());
    for i in 0..table.len() {
        let mut bv = BitVector::new(table.len());
        bv.xor_bit(i);
        augmented_matrix.push(bv);
    }
    loop {
        if let Some(y) = kernel(&mut matrix, &mut augmented_matrix, &mut pivots) {
            let mut map = HashMap::new();
            let parity = y.popcount() & 1 == 1;
            for i in 0..table.len() {
                if parity && !y.get(i) { 
                    map.insert(table[i].get_integer_vec(), 1);
                }
                else if !parity && y.get(i) { 
                    map.insert(table[i].get_integer_vec(), 1);
                }
            }
            for i in 0..table.len() {
                if !y.get(i) { continue; }
                for j in 0..table.len() {
                    if y.get(j) { continue; }
                    let mut z = table[i].clone();
                    z.xor(&table[j]);
                    let score = map.entry(z.get_integer_vec()).or_insert(0);
                    *score += 2;
                }
            }
            let mut max_cost = 0;
            let mut max_key: Option<&Vec<i128>> = None;
            for (key, val) in map.iter() {
                if *val > max_cost || (*val == max_cost && max_key.is_some() && key < max_key.unwrap()) {
                    max_cost = *val;
                    max_key = Some(key);
                }
            }
            if max_cost <= 0 { break; }
            let z = BitVector::from_integer_vec(max_key.unwrap().to_vec());
            let mut to_update = y.get_boolean_vec();
            to_update.truncate(table.len());
            if (y.popcount() & 1) == 1 {
                table.push(BitVector::new(table[0].size() - 1));
                matrix.push(BitVector::new(table[0].size() - 1));
                let mut bv = BitVector::new(table.len());
                bv.xor_bit(augmented_matrix.len());
                augmented_matrix.push(bv);
                to_update.push(true);
            }
            let indices = to_update.iter().enumerate().filter(|(_, &r)| r).map(|(index, _)| index).collect::<Vec<_>>();
            for i in indices {
                table[i].xor(&z);
            }
            let mut to_remove = to_remove(&table);
            to_remove.sort_by(|a, b| b.cmp(a));
            for i in to_remove {
                clear_column(i, &mut matrix, &mut augmented_matrix, &mut pivots);
                table.swap_remove(i);
                matrix.swap_remove(i);
                augmented_matrix.swap_remove(i);
                to_update.swap_remove(i);
                if pivots.contains_key(&table.len()) {
                    let tmp = pivots.remove(&table.len()).unwrap();
                    pivots.insert(i, tmp);
                }
                for j in 0..augmented_matrix.len() {
                    if augmented_matrix[j].get(i) != augmented_matrix[j].get(table.len()) {
                        augmented_matrix[j].xor_bit(i);
                    }
                    if augmented_matrix[j].get(table.len()) {
                        augmented_matrix[j].xor_bit(table.len());
                    }
                }
            }
            let size = BitVector::new(table.len()).size();
            for i in 0..table.len() {
                while augmented_matrix[i].size() > size {
                    augmented_matrix[i].blocks.pop();
                }
            }
            let indices = to_update.iter().enumerate().filter(|(_, &r)| r).map(|(index, _)| index).collect::<Vec<_>>();
            for i in &indices {
                clear_column(*i, &mut matrix, &mut augmented_matrix, &mut pivots);
                matrix[*i] = table[*i].clone();
                let mut bv = BitVector::new(table.len());
                bv.xor_bit(*i);
                augmented_matrix[*i] = bv;
                let mut t_vec = table[*i].get_boolean_vec();
                t_vec.truncate(nb_qubits);
                let mut vec = Vec::<bool>::new();
                for _ in 0..nb_qubits {
                    if t_vec.pop().unwrap() {
                        vec.append(&mut t_vec.clone());
                    }
                    else {
                        vec.append(&mut vec![false; t_vec.len()]);
                    }
                }
                matrix[*i].extend_vec(vec, nb_qubits);
            }
        }
        else { break; }
    }
    table
}

pub fn fast_todd(mut table: Vec<BitVector>, nb_qubits: usize) -> Vec<BitVector> {
    loop {
    table = tohpe(table.clone(), nb_qubits);
    let mut matrix = table.clone();
    for i in 0..table.len() {
        let mut t_vec = table[i].get_boolean_vec();
        t_vec.truncate(nb_qubits);
        let mut vec = Vec::<bool>::new();
        for _ in 0..nb_qubits {
            if t_vec.pop().unwrap() {
                vec.append(&mut t_vec.clone());
            }
            else {
                vec.append(&mut vec![false; t_vec.len()]);
            }
        }
        matrix[i].extend_vec(vec, nb_qubits);
    }
    let mut pivots: HashMap::<usize, usize> = HashMap::new();
    let mut augmented_matrix: Vec<BitVector> = Vec::with_capacity(table.len());
    for i in 0..table.len() {
        let mut bv = BitVector::new(table.len());
        bv.xor_bit(i);
        augmented_matrix.push(bv);
    }
    kernel(&mut matrix, &mut augmented_matrix, &mut pivots);
    pivots = pivots.iter().map(|(k, v)| (*v, *k)).collect();

    let mut map = HashMap::new();
    for i in 0..table.len() {
        map.insert(table[i].get_integer_vec(), i);
    }
    let mut max_score = 0;
    let mut max_z = None;
    let mut max_y = None;
    for i in 0..table.len() {
        for j in i+1..table.len() {
            let mut z = table[i].clone();
            z.xor(&table[j]);
            let z_vec = z.get_boolean_vec();
            let mut r_mat: Vec<BitVector> = Vec::new();
            let mut augmented_r_mat: Vec<BitVector> = Vec::new();
            for k in 0..nb_qubits {
                let mut col = BitVector::new_block_size(matrix[0].blocks.len());
                let mut augmented_col = BitVector::new_block_size(augmented_matrix[0].blocks.len());
                let mut l = 0;
                for a in (0..nb_qubits).rev() {
                    for b in 0..a {
                        if (a == k && z_vec[b]) || (b == k && z_vec[a]) {
                            col.xor_bit(nb_qubits + l);
                            if pivots.contains_key(&(nb_qubits + l)) {
                                let val = pivots.get(&(nb_qubits + l)).unwrap();
                                col.xor(&matrix[*val]);
                                augmented_col.xor(&augmented_matrix[*val]);
                            }
                        }
                        l += 1;
                    }
                }
                r_mat.push(col);
                augmented_r_mat.push(augmented_col);
            }

            let mut col = BitVector::new_block_size(matrix[0].blocks.len());
            let mut augmented_col = BitVector::new_block_size(augmented_matrix[0].blocks.len());
            let mut l = 0;
            for a in (0..nb_qubits).rev() {
                for b in 0..a {
                    if z_vec[a] && z_vec[b] {
                        col.xor_bit(nb_qubits + l);
                        if pivots.contains_key(&(nb_qubits + l)) {
                            let val = pivots.get(&(nb_qubits + l)).unwrap();
                            col.xor(&matrix[*val]);
                            augmented_col.xor(&augmented_matrix[*val]);
                        }
                    }
                    l += 1;
                }
                if z_vec[a] {
                    col.xor_bit(a);
                    if pivots.contains_key(&a) {
                        let val = pivots.get(&a).unwrap();
                        col.xor(&matrix[*val]);
                        augmented_col.xor(&augmented_matrix[*val]);
                    }
                }
            }
            r_mat.push(col);
            augmented_r_mat.push(augmented_col);

            for k in 0..r_mat.len() {
                let index = r_mat[k].get_first_one();
                if r_mat[k].get(index) {
                    let pivot = r_mat[k].clone();
                    let augmented_pivot = augmented_r_mat[k].clone();
                    for l in (k+1)..r_mat.len() {
                        if r_mat[l].get(index) {
                            r_mat[l].xor(&pivot);
                            augmented_r_mat[l].xor(&augmented_pivot);
                        }
                    }
                }
                else if augmented_r_mat[k].get(i) ^ augmented_r_mat[k].get(j) {
                    let mut score = 0;
                    let y = augmented_r_mat[k].clone();
                    for l in 0..table.len() {
                        if y.get(l) {
                            table[l].xor(&z);
                            if map.contains_key(&table[l].get_integer_vec()) && !y.get(*map.get(&table[l].get_integer_vec()).unwrap()) {
                                score += 2;
                            }
                            table[l].xor(&z);
                        }
                    }
                    if y.popcount() & 1 == 1 {
                        if map.contains_key(&z.get_integer_vec()) {
                            score += 1;
                        }
                        else {
                            score -= 1;
                        }
                    }
                    if score > max_score {
                        max_score = score;
                        max_z = Some(z.clone());
                        max_y = Some(y);
                    }
                }
            }
        }
    }
    if max_score == 0 { break; }
    let y = max_y.unwrap();
    let z = max_z.unwrap();
    for l in 0..table.len() {
        if y.get(l) {
            table[l].xor(&z);
        }
    }
    if y.popcount() & 1 == 1 {
        table.push(z);
    }
    table = proper(table);
    }
    table
}
