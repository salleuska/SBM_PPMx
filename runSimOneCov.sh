#!/bin/bash

alphas="0 0.5 1 2 4"

# 1-cov scenarios (x1 used; x2 is plotting noise)
scenarios="neutral informative mislead_random mislead_shifted"

for sc in $scenarios; do
  parallel -j 8 --verbose \
    Rscript simulation/1_runSimulationBinaryESBM_oneCovariate.R \
      data_path=simulation/data/binarySBM_1cov_${sc}.rds \
      alpha={} \
      seed=1 \
    ::: $alphas
done