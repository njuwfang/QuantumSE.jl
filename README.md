# QuantumSE.jl

A prototype tool for symbolic execution of quantum programs (QSE) with symbolic stabilizer states.

By symbolizing the phases of stabilizer states, **symbolic stabilizer states**  enable us to use symbolic expressions
to characterize the possible adversarial errors in quantum error correction (QEC) programs.
In this way, QSE with symbolic stabilizer states facilitate the efficient analysis of QEC programs.

See [the arXiv paper](https://arxiv.org/abs/2311.11313) for more details.

## Installation

To install it, run [Julia 1.10+](https://julialang.org/downloads/) REPL and use:

###### Due to julia's upgrade, there may be problems with the Z3 dependency
```
] add https://github.com/njuwfang/QuantumSE.jl.git
```

#### SMT Solver

QuantumSE.jl uses [Bitwuzla](https://github.com/bitwuzla/bitwuzla) as the default solver.

To install our experimental version, please follow the instructions below:
 
##### Manual installation (Ubuntu 22.04, 20.04)

1. Install required dependencies:
    
    ```bash
    sudo apt-get install python3 python3-pip pkg-config m4 libgmp-dev
    ```

    ```bash
    pip install meson ninja
    ```
2. Clone our forked repo of Bitwuzla and cd to it.

    ```bash
    git clone https://github.com/njuwfang/bitwuzla-for-QuantumSE.git && cd bitwuzla-for-QuantumSE
    ```
3. Configure it with `--kissat`
    
    ```bash
    ./configure.py --kissat
    ```
4. cd to `./build` and build
    
    ```bash
    cd ./build && ninja
    ```
5. Install `bitwuzla`.

    ```bash
    sudo ninja install
    ```

## Evaluations

We have evaluated QuantumSE.jl in [the arXiv paper]().

Clone this repo and cd to `./QuantumSE.jl/example`.

```bash
git clone https://github.com/njuwfang/QuantumSE.jl.git && cd ./QuantumSE.jl/example
```

#### Finding Bugs in QEC Programs

1. Repetition codes:
    ```bash
    julia RepetitionCode.jl
    ```

2. Toric codes:
    ```bash
    julia ToricCode.jl
    ```
3. Quantum Tanner codes:
    ```
    julia TannerCode.jl
    ```
The performance results are stored in `.dat` files.

#### Comparing [symQV](https://github.com/fabianbauermarquart/symQV)

1. Clone symQV's repo and install its dependencies.
    ```bash
    git clone https://github.com/fabianbauermarquart/symQV.git && ./symQV/install.sh
    ```
2. Run scripts.
    ```bash
    python3 comp_symqv.py
    ```

    ```bash
    julia comp_quantumse.jl
    ```

The performance results are stored in `.dat` files.

#### Comparing [Stim](https://github.com/quantumlib/Stim)

1. Install Stim (as of 7/11/2023, the latest stable version is 1.12.0) with `pip`.
    ```bash
    pip install stim==1.12.0
    ```
2. Generate benchmark files (in .stim format).
    ```bash
    julia generate_randomcircuits.jl
    ```
3. Benchmark Stim and QuantumSE.jl.
    ```bash
    python3 benchmark_stim.py
    ```

    ```bash
    julia benchmark_quantumse.jl
    ```

    The benchmark results are stored in `.dat` files.
