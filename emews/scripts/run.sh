#! /bin/bash
R_PATH=$1
INPUT_PARAMS=$2

READLINE_LIB="/lcrc/project/EMEWS/improv/sfw/gcc-13.2.0/readline-8.2/lib"
GETTEXT_LIB="/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-13.2.0/gettext-0.22.3-ycr7uf7/lib"
export LD_LIBRARY_PATH=$READLINE_LIB:$GETTEXT_LIB:$LD_LIBRARY_PATH
export PATH=/lcrc/project/EMEWS/improv/sfw/gcc-13.2.0/R-4.3.2/lib64/R/bin:$PATH

export MKL_NUM_THREADS=$CFG_MKL_THREADS
echo "MKL_NUM_THREADS: $MKL_NUM_THREADS"

Rscript $R_PATH/run_experiment.R $INPUT_PARAMS
