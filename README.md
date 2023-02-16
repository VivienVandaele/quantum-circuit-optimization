Hadamard and internal Hadamard gates optimization based on the paper [Optimal Hadamard gate count for Clifford+T synthesis of Pauli rotations sequences](https://arxiv.org/abs/2302.07040).

### Installation
[Install rust](https://www.rust-lang.org/tools/install) and run ```cargo +nightly run -r circuits/inputs/tof_3.qc```
where ```circuits/inputs/tof_3.qc``` can be replaced by any .qc file. The optimized circuit will be written in the .qc format in the folder ```circuits/outputs/```.
