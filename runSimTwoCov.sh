#!/bin/bash

alphas="0 0.5 1 2"

# Two-covariate scenarios: NN, IN, NI, II
scenarios="NN IN NI II"

for sc in $scenarios; do
  parallel -j 8 --verbose \
    Rscript simulation/1_runSimulationBinaryESBM_twoCov.R \
      data_path=simulation/data/binarySBM_2cov_${sc}.rds \
      alpha1={1} \
      alpha2={2} \
      seed=1 \
    ::: $alphas \
    ::: $alphas
done