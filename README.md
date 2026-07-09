# SLA — Swimming-Limited Aggregation

MATLAB code for agent-based simulations and figure generation accompanying the preprint:

> Sintès, G., Goral, M., López-León, T., Lindner, A., & Tătulea-Codrean, M. (2026). Swimming-limited aggregation of bacteria in liquid crystals. *arXiv* preprint arXiv:2607.05239. https://doi.org/10.48550/arXiv.2607.05239

## Requirements

- MATLAB (any recent release)

## Repository structure

├── Experiments/

│   ├── Codes/             # Scripts to process experimental data

│   │   └── Output/        # Processed experimental data (generated)

│   └── Data/              # Raw experimental data

├── Simulations/

│   └── Codes/             # Scripts to reproduce simulation results

│       └── Functions/     # Required functions

│   └── Data/              # Simulation results (generated)

├── Figures/               # Scripts to reproduce preprint subfigures

│   └── Output/            # Subfigures (generated)

## How to use

Run the following steps in order.

**1. Process the experimental data.**
Run the scripts in `Experiments/Codes`. These analyze the raw data in `Experiments/Data` and write the processed results to `Experiments/Codes/Output`, ready for figure generation.

**2. Generate the simulation results.**
Run the scripts in `Simulations/Codes`. These produce the primary agent-based simulation results and save them to `Simulations/Data`. Sample data required for figure generation is already included in this folder.

**3. Reproduce the figures.**
Run the scripts in `Figures` to reproduce the subfigures from the preprint.

## Citation

If you use this code, please cite our preprint:

> Sintès, G., Goral, M., López-León, T., Lindner, A., & Tătulea-Codrean, M. (2026). Swimming-limited aggregation of bacteria in liquid crystals. *arXiv* preprint arXiv:2607.05239. https://doi.org/10.48550/arXiv.2607.05239
