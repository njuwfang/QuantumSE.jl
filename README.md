# QuantumSE.jl

A prototype tool for symbolic execution of quantum programs (QSE) with symbolic stabilizer states.

By symbolizing the phases of stabilizer states, **symbolic stabilizer states**  enable us to use symbolic expressions
to characterize the possible adversarial errors in quantum error correction (QEC) programs.
In this way, QSE with symbolic stabilizer states facilitate the efficient analysis of QEC programs.

See [the arXiv paper](https://arxiv.org/abs/2311.11313) for more details.

## Evaluations

We have evaluated QuantumSE.jl in [the arXiv paper](https://arxiv.org/abs/2311.11313).

Clone this repo and cd to `./QuantumSE.jl/example`.

```bash
git clone https://github.com/njuwfang/QuantumSE.jl.git && cd ./QuantumSE.jl/example
```

## Installation

To fully utililize the capabilities of `QuantumSE` package, install [Julia 1.10+](https://julialang.org/downloads/) and `QuantumSE` and `Bitwuzla` packages respectively by following the instructions below. Its recommended to familiarize yourself with [Julia's REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), especially "julia", "pkg" and "shell" modes; and [Julia's environments](https://docs.julialang.org/en/v1/manual/code-loading/#Environments)

### QuantumSE package

#### Option 1: Install QuantumSE from source

In julia's REPL pkg mode(activate by pressing `]`), Add the relevant packages

```
(@v1.10) pkg> add z3_jll#z3-v4.12.4+0 Z3#b0ff00a https://github.com/njuwfang/QuantumSE.jl.git
```

> NOTE: Due to julia's upgrade, there may be problems with the Z3 dependency

#### Option 2: Activate QuantumSE by cloning source code

Clone this repo and cd to `./QuantumSE.jl/`.
```bash
git clone https://github.com/njuwfang/QuantumSE.jl.git && cd ./QuantumSE.jl
```


Then in julia's REPL pkg mode, execute:

```bash
(@v1.10) pkg> activate .
(QuantumSE) pkg> instantiate
```

### SMT Solver

QuantumSE.jl uses [Bitwuzla](https://github.com/bitwuzla/bitwuzla) as the default solver.

To install our experimental version, please follow the instructions below:
 
### Manual installation (Ubuntu 22.04, 20.04)

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

> NOTE: If using ubuntu in WSL, open Julia REPL under root privilege and execute these commands under shell mode of Julia REPL.


### Finding Bugs in QEC Programs

> NOTE: If running ubuntu in WSL, after installing `QuantumSE` and `Bitwuzla` in julia REPL, run `julia> include("<path to example>\TannerCode.jl")` and similarly for others.

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

### Comparing [symQV](https://github.com/fabianbauermarquart/symQV)

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

### Comparing [Stim](https://github.com/quantumlib/Stim)

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
