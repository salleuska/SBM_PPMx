#!/bin/bash

# Application run: 3 covariates (x,y,z) with GN prior
# Vary alpha_x, alpha_y, alpha_z over the same grid

alphas="0 0.5 1 2"

DATA_PATH="application_brain/consensus_scale33_tau50.RData"
SCRIPT="application_brain/1_runBrainApplication.R"   # <-- set to your script name

parallel -j 8 --verbose \
  Rscript "$SCRIPT" \
    data_path="$DATA_PATH" \
    alpha_x={1} \
    alpha_y={2} \
    alpha_z={3} \
    seed={#} \
  ::: $alphas \
  ::: $alphas \
  ::: $alphas