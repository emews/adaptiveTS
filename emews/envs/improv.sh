module load gcc/13.2.0 intel-mkl/2020.4.304-gcc-13.2.0

READLINE_LIB="/lcrc/project/EMEWS/improv/sfw/gcc-13.2.0/readline-8.2/lib"
GETTEXT_LIB="/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-13.2.0/gettext-0.22.3-ycr7uf7/lib"

export LD_LIBRARY_PATH=$READLINE_LIB:$GETTEXT_LIB:$LD_LIBRARY_PATH

R=/lcrc/project/EMEWS/improv/sfw/gcc-13.2.0/R-4.3.2/lib64/R
PGSQL=/lcrc/project/EMEWS/improv/sfw/gcc-13.2.0/postgresql-12.17
export PATH=$R/bin:$PGSQL/bin:$PATH

# export MKL_NUM_THREADS=1
# export OMP_NUM_THREADS=$MKL_NUM_THREADS
