Supplementary code to the paper XXX by ABCD authors

## Preliminaries
Install packages in the `install.R` script. This repo uses the R library `here` for handlings project file paths reproducibly

## Description

The directory `R_utilities` contains the file `esbm_original.R` and `stirling.cppp` from 
[https://github.com/danieledurante/ESBM](https://github.com/danieledurante/ESBM)
which include functions for the posterior sampling of stochastic block model for binary data with a categorical covariate

- **`esbm.R`**  
  A modified version of the ESBM sampler.  
  The key addition is the ability to **pass an arbitrary similarity/cohesion function** as an argument.  
  This allows us to replace the fixed categorical similarity of the original ESBM with PPMx-style covariate‐dependent cohesion terms.

  Arguments for the similarity are provided through a list (`sim_args`), so different similarity models can be supplied without altering the sampler.

- **`similarity_functions.R`**  
  Contains the interface and implementations of the allowed similarity functions. Each function takes a covariate vector/matrix for a cluster and returns a scalar similarity value.  

  **2025-11-12 update:**  
  Implemented a similarity based on a **normal auxiliary model**:  
  - observation model: \( x_i \mid \mu_h \sim N(\mu_h, 1) \)  
  - prior on cluster mean: \( \mu_h \sim N(m_0, s_0^2) \)  
  This provides a continuous-covariate cohesion that can be raised to a power parameter \( \alpha \) for calibration.

### Planned extensions

We intend to extend the code to support:

- Poisson likelihood for weighted networks (with and without self-loops);
- cohesion for continuous covariates derived from full normal-normal marginalization;
- circular covariates (e.g., von Mises-based cohesion);
- weighted similarity functions and mixtures of similarities.

---
- **`esbm.R`**  
  Contains modified version of the ESBM sampler, with the possibility to pass one or more arbitrary similarities functions as an argument. Note that I have not checked that the postprocessing functions work.

- **`similarity_functions.R`**  
  Contains the interface and implementations of the allowed similarity functions. Each function takes a covariate vector/matrix for a cluster and returns a scalar similarity value.  

  **2025-11-12 update:**  
  Implemented a similarity based on a **normal auxiliary model**:  
  - observation model: \( x_i \mid \mu_h \sim N(\mu_h, 1) \)  
  - prior on cluster mean: \( \mu_h \sim N(m_0, s_0^2) \)  
  This provides a continuous-covariate cohesion that can be raised to a power parameter \( \alpha \) for calibration.

The `test_sim_ESBM` folder contains just some code to verify the new `esbm` function was working

### Planned extensions

We intend to extend the code to support:

- Poisson likelihood for weighted networks (with and without self-loops);
- cohesion for continuous covariates derived from full normal-normal marginalization;
- circular covariates (e.g., von Mises-based cohesion);
- weighted similarity functions and mixtures of similarities.

---

## Notes on simulations

All simulation code is located in the `simulation/` directory.

### Data generation

**`0_simulateBinaryESBM.R`**  
Creates synthetic binary networks under different covariate scenarios.  
Output stored as `.rds` includes:
- adjacency matrix `Y`;
- covariates `x`;
- true cluster labels `z_true`;
- scenario metadata (e.g., one informative covariate, two covariates, or no covariate signal).

### Running the simulation study

Two scripts run the ESBM under different settings:

- **`1_runSimulationBinaryESBM_oneCovariate.R`**  
  Uses only one covariate and tests different similarity strengths or models.

- **`1_runSimulationBinaryESBM.R`**  
   Version with both covariantes

These scripts:

1. Parse input arguments (`data_path`, `alpha`, `seed`).  
2. Load the dataset and the sampler.  
3. Run ESBM with the chosen similarity and \( \alpha \).  
4. Save results in `simulation/results/[oneCov/twoCov]`.

Note: I run the scripts from bash in parallel via `runSimOneCov.sh` 



- **`2_resultsSimulationBinaryESBM.R`**
Collects all outputs and computes:
	ARI vs true labels,
	effect of the similarity parameter ( \alpha ) by scenario.


