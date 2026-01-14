#!/bin/bash

# Arguments: input_csv output_csv nsamples
INPUT_CSV=$1
OUTPUT_CSV=$2
NSAMPLES=$3

# Virtual environment directory
# source activate uqpy

CONDA="/home/fadikar/softwares/anaconda3/condabin/conda"

# Run the Python script
# python flow_sampler.py --input "$INPUT_CSV" --output "$OUTPUT_CSV" --nsamples $NSAMPLES
$CONDA run -n uqpy python flow_sampler.py --input "$INPUT_CSV" --output "$OUTPUT_CSV" --nsamples $NSAMPLES

# Deactivate
# conda deactivate
