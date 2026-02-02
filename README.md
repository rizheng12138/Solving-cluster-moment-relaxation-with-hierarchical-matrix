# Solving Cluster Moment Relaxation with Hierarchical Matrix (MATLAB)

This repository provides a MATLAB reference implementation and reproducibility scripts for the paper:

**Solving cluster moment relaxation with hierarchical matrix** (open resource).

The codebase is organized so that most users can directly run the main script to reproduce benchmark results for pre-generated instances (N = 64, 128, 256, 512, 1024).
If you want to test larger sizes or regenerate the datasets and initializations, you can run the provided preparation scripts.

---

## Contents

- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Quick Start](#quick-start)
- [Running Different Problem Sizes](#running-different-problem-sizes)
- [Outputs](#outputs)
- [Regenerating Data and Initialization (Optional)](#regenerating-data-and-initialization-optional)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Repository Structure

- data/ pre-generated datasets and initializations
- results/ saved benchmark outputs
- scripts/ entry scripts (prepare data / init / main)
- src/ core functions (solver, objective, gradients, operators, IO helpers)

**Key entry points**

- `scripts/main.m` — run the solver and save results to `results/`
- `scripts/prepare_data.m` — (optional) generate `data/<N>_data.mat`
- `scripts/init_S0.m` — (optional) generate `data/<N>_S0.mat`

**Core implementation**

- `src/run_solver.m` — main solver routine (Algorithm 3)
- `src/H.m` — hierarchical PSD representation (Eq. (25) in the paper)
- `src/Loss.m`, `src/dLoss.m`, `src/prepareCG.m` — Manopt cost/gradient + caching
- `src/build_VW.m` — constructs sparse operators used in gradients
- `src/data_filename.m`, `src/s0_filename.m`, `src/result_filename.m` — file naming helpers
- `src/set_default_opts.m` — centralized configuration
- `src/get_J.m` — get  `J`  in the primal objective  `tr(JM)`
- `grad_S0.m` — gradient used in  `init_S0.m`

---

## Dependencies

- MATLAB (tested on R2024a)
- **Manopt** (required): [Manopt](https://www.manopt.org/)

> Manopt must be on your MATLAB path for `rlbfgs`, `productmanifold`, etc.

---

## Quick Start

1. Clone this repository.
2. Open MATLAB and set the **current folder** to the repository root.
3. Make sure **Manopt** is on your path.
4. Run:

```matlab
run('scripts/main.m');
```

By default,  `main.m`  runs with  `N = 64`  and writes:

```
results/64_benchmark.mat
```

### Running Different Problem Sizes

Edit the options line in  `scripts/main.m`:

```
opts = set_default_opts('N', 128);
```

Available pre-generated sizes in data/ (as provided in this repo):

```
N = 64, 128, 256, 512, 1024
```

If you choose an  `N`  that is not pre-generated, see:
[Regenerating Data and Initialization](#regenerating-data-and-initialization)

## Outputs

Each run saves a MATLAB `.mat` file:

```
results/<N>_benchmark.mat
```

Typical fields include:

- `eta_P`, `eta_D`, `eta_g`: feasibility and gap measures per outer iteration
- `energy`, `energy_per_site`: objective history / per-site changes
- `iter`, `time`, `relative_err`: summary statistics
- algorithm parameters used in the run (e.g., `mu`, `tau`, `maxiter_opt`, `r_l`)

## Regenerating Data and Initialization (Optional)

Most users do **not** need this section because `data/` already contains:

- `<N>_data.mat` for `N = 64, 128, 256, 512, 1024`
- `<N>_S0.mat` for the same sizes

If you want to test a new `N` (e.g., `2048`) or rebuild files:

### Step 1: Generate data

```matlab
run('scripts/prepare_data.m');
```

This produces:

```
data/<N>_data.mat
```

### Step 2: Generate initialization

```
run('scripts/init_S0.m');
```

This produces:

```
data/<N>_S0.mat
```

### Step 3: Run the solver

```
run('scripts/main.m');
```

These scripts use  `set_default_opts(...)`  to control  `N, h`, and other parameters.

## Troubleshooting

### 1) `Undefined function or variable 'rlbfgs'` / Manopt not found

Manopt is not on your MATLAB path. Install Manopt and add it to your path, e.g.:

```matlab
addpath(genpath('/path/to/manopt'));
```

### 2) Missing  `data/<N>_data.mat`  or  `data/<N>_S0.mat`

You selected an  `N`  that is not pre-generated. Run:

```
run('scripts/prepare_data.m');
run('scripts/init_S0.m');
```

### 3) Output is too verbose

Manopt verbosity can be reduced by setting  `options.verbosity = 0`  inside the relevant script/function.

## Citation

If you use this code in academic work, please cite the paper:

- Yi Wang, Rizheng Huang, Yuehaw Khoo, *Solving cluster moment relaxation with hierarchical matrix*, **Journal of Computational Physics**, 541 (2025) 114331.
  DOI: https://doi.org/10.1016/j.jcp.2025.114331
  ScienceDirect: https://www.sciencedirect.com/science/article/pii/S0021999125006138

### BibTeX

```bibtex
@article{wang2025cluster_moment_hmatrix,
  title   = {Solving cluster moment relaxation with hierarchical matrix},
  author  = {Wang, Yi and Huang, Rizheng and Khoo, Yuehaw},
  journal = {Journal of Computational Physics},
  volume  = {541},
  pages   = {114331},
  year    = {2025},
  doi     = {10.1016/j.jcp.2025.114331},
  url     = {https://www.sciencedirect.com/science/article/pii/S0021999125006138}
}
```

## License

This project is released under the MIT License. See the `LICENSE` file for details.
