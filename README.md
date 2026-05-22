# Echo State Network for Hamiltonian Dynamics Prediction

Codes for the paper: **Inferring measure synchronization in coupled bosonic Josephson junctions with reservoir computing**

## Overview

This repository contains MATLAB codes for inferring measure synchronization (MS) in coupled bosonic Josephson junctions using parameter-aware reservoir computing (RC). The RC model is trained on time series data and can infer critical coupling strengths for MS transitions as well as the variation of system order parameters.

## Repository Structure

```
codes/
├── fig4codes/       # Codes for Figure 4
│   ├── he.m         # ESN training & prediction script
│   ├── para2.mat    # Optimized hyperparameters
│   ├── sl.mat       # Training data (phase space)
│   └── model.txt    # Analytical model reference data
├── fig6codes/       # Codes for Figure 6
│   ├── he.m
│   ├── para.mat
│   ├── sl.mat
│   └── model.txt
└── fig7codes/       # Codes for Figure 7
    ├── he.m
    ├── para.mat
    ├── sl.mat
    └── model.txt
```

## Method

The Echo State Network is a type of recurrent neural network with:
- A randomly initialized, fixed recurrent reservoir
- Only the output weights are trained via ridge regression
- Squared activation for even-indexed reservoir nodes

### ESN Update Rule

```
x(t+1) = (1-a) * x(t) + a * tanh(W_in * u(t) + W * x(t))
```

where:
- `x(t)` — reservoir state
- `u(t)` — input
- `W_in` — input weight matrix
- `W` — reservoir weight matrix (sparse, with controlled spectral radius)
- `a` — leaking rate

### Training

Output weights `W_out` are computed in closed form via ridge regression:

```
W_out = Y * X^T / (X * X^T + reg * I)
```

Transient segments are discarded to remove initial transient effects.

## Requirements

- MATLAB (tested with R2020b or later)
- No additional toolboxes required

## Usage

1. Open MATLAB and navigate to the desired figure directory (e.g., `codes/fig4codes/`)
2. Run `he.m`
3. The script will:
   - Load pre-optimized hyperparameters from `para.mat`
   - Load training data from `sl.mat`
   - Train the ESN
   - Run prediction with varying control parameters
   - Plot results comparing ESN predictions with analytical model

## Hyperparameters

Each figure uses pre-optimized hyperparameters stored in `para.mat` / `para2.mat`:

| Parameter | Description |
|-----------|-------------|
| `eig_rho` | Spectral radius of the reservoir |
| `W_in_a` / `W_in_b` | Input weight scaling |
| `a` | Leaking rate |
| `reg` | Ridge regression regularization |
| `density` | Reservoir connection density |

## Citation

If you use this code, please cite the associated paper.

## License

This code is provided for research purposes.
