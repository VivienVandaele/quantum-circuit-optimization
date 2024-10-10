Rust implementation of quantum circuits optimization algorithms presented in the following papers:
- BBMerge and FastTMerge algorithms from [Optimal number of parametrized rotations and Hadamard gates in parametrized Clifford circuits with non-repeated parameters](https://arxiv.org/abs/2407.07846)
- InternalHOpt algorithm from [Optimal Hadamard gate count for Clifford+T synthesis of Pauli rotations sequences](https://arxiv.org/abs/2302.07040)
- TOHPE and FastTODD algorithms from [Lower T-count with faster algorithms](https://arxiv.org/abs/2407.08695)

### Usage
[Install rust](https://www.rust-lang.org/tools/install) and run ```cargo run -r [OPTIONS] file.qc```
where ```file.qc``` is any .qc file. If your machine supports AVX2, it can be enabled with ```RUSTFLAGS="-C target-cpu=native" cargo run -r [OPTIONS] file.qc```.

```[OPTIONS]``` are (case-insensitive and in no order):
- ```BBMerge``` runs the BBMerge algorithm
- ```FastTMerge``` runs the FastTMerge algorithm
- ```InternalHOpt``` runs the InternalHOpt algorithm
- ```TOHPE``` runs the TOHPE algorithm
- ```FastTODD``` runs the FastTODD algorithm

If no options are provided, then the FastTMerge, InternalHOpt and FastTODD algorithms will be applied.
The gadgetization of internal Hadamard gates will be done whenever the TOHPE or FastTODD algorithms are applied.
The optimized circuit will be written in the .qc format in the folder ```circuits/outputs/```.
