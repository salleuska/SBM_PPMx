Supplementary code for an ESBM/PPMx project on stochastic block models with covariate-informed cohesion.

## Preliminaries

Install the required R packages first. The project uses `here` to make file paths reproducible from the repository root.

## Repository layout

The core sampler code lives in the legacy folder `R_utilties/` (note the existing folder spelling in the repo).

- `R_utilties/esbm.R` is the main sampler used in this repository. It supports running with or without covariates and allows one or more similarity functions.
- `R_utilties/stirling.cpp` provides the compiled helper routines used by the sampler.
- `R_utilties/similarity_functions.R` defines the supported similarity functions and their arguments.
- `test_sim_ESBM/` contains older sanity-check scripts used while extending the sampler.
- `application_brain/` contains the applied analysis workflow for the brain-network example.
- `simulation/` contains synthetic-data generation, simulation runners, result processing, and plotting scripts.

### Continuous-covariate similarity

The repo includes a Gaussian auxiliary-model similarity for continuous covariates:

- observation model: \(x_i \mid \mu_h \sim N(\mu_h, 1)\)
- prior on the cluster mean: \(\mu_h \sim N(m_0, s_0^2)\)

This cohesion can be calibrated through a learnable or fixed power parameter `alpha`.

## Simulation workflow

All simulation code is under `simulation/`.

### Data generation

`simulation/0_simulateBinaryESBM.R` creates binary-network datasets under several covariate scenarios. Each `.rds` file stores:

- adjacency matrix `Y`
- covariates `x`
- true partition `partition`
- scenario metadata such as number of covariates and whether they are informative or misleading

### Running simulation experiments

The main learn-alpha runner is:

- `simulation/learn_alpha/1_runSimulationBinaryESBM.R`

It accepts command-line arguments such as:

- `data_path`
- `seed`
- `N_iter`
- `model` (`ESBM` or `SBM`)
- `similarity_calibration`
- `alpha_init`, `a_alpha`, `b_alpha`

Outputs are written under:

- `simulation/learn_alpha/<calibration>/results/<nCov>Cov` for ESBM runs
- `simulation/learn_alpha/SBM/results/<nCov>Cov` for network-only baseline runs

The shell scripts in the repository are the intended orchestration layer for launching batches of runs in parallel.

### Processing simulation outputs

Learn-alpha ESBM outputs can be combined with:

- `simulation/learn_alpha/2_resultsSimulationBinaryESBM.R`

Network-only SBM baseline outputs can be combined with:

- `simulation/learn_alpha/2_processSBMresults.R`

These processing scripts compute posterior summaries and clustering diagnostics, including:

- posterior similarity matrices
- SALSO point estimates
- adjusted Rand index (ARI)
- silhouette summaries based on `1 - PSM`

## Notes

This is still a research codebase rather than a polished R package. The simulation and plotting pipelines are active, but documentation and automated testing are lighter than the core methodology work.
