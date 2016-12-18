#!/bin/bash

cd "/homes/mb2215/bitbucket/Thesis/bilinear-rpca/hypersearch"
matlab -nodisplay -nosplash -r "genplot_err_func_lambda_mu_fine($1, $2, $3)"
