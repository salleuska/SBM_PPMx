Supplementary code to the paper XXX by ABCD authors

## Preliminaries
Install packages in the `install.R` script. This repo uses the R library `here` for handlings project file paths reproducibly

## Description

The directory `R_utilities` contains the file `esbm_original.R` and `stirling.cppp` from 
[https://github.com/danieledurante/ESBM](https://github.com/danieledurante/ESBM)
which include functions for the posterior sampling of stochastic block model for binary data with a categorical covariate

The file `esbm.R` modifies the original gibbs function to allow passing as an argument a general similarity function and eventual parameters of the similarity funciton. See `R_utilities/similiary_functions.R` for details on the format and implemented functions. 

20251112 - implemented a similarity based on a normal auxiliary model with variance fixed to 1 with normal prior on the mean only. 

**Note** we plan to build on this and add:
- Poisson likelihood for weighted network (with and withouth self-loop)
- cohesion function for continuous covariate (derived from normal distribution)
- cohesion function for covariates defined on a circle 
- weighted cohesion function

## Notes on simulation ideas
