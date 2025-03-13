
#!/bin/bash

if [ "$#" -ne 2 ]; then
  script_name=$(basename $0)
  echo "Usage: ${script_name} walltime exp_id"
  exit 1
fi

WALLTIME=$1
JOB_NAME=$2
THIS=$( cd $( dirname $0 ) ; /bin/pwd )
OUTPUT_FILE=$THIS/${JOB_NAME}_output.txt

echo $OUTPUT_FILE

# | qsub -
cat <<EOF | qsub -
#!/bin/bash -l

# UG Section 2.5, page UG-24 Job Submission Options
# Note: Command line switches will override values in this file.

# ----------------- PBS Directives ----------------- #
# These options are mandatory in LCRC; qsub will fail without them.

# select: number of nodes. Adjust these values as per your job's requirement. In the below example, 4 nodes are requested.
# ncpus: number of cores per node. (use 128 unless you have a reason otherwise)
# mpiprocs: number of MPI processes (use the same value as ncpus unless you have a reason otherwise ).
# walltime: a limit on the total time from the start to the completion of a job
#PBS -A EMEWS
#PBS -l select=1:ncpus=1

#PBS -l walltime=$WALLTIME


# Queue for the job submission
#PBS -q compute

# Job name (first 15 characters are displayed in qstat)
#PBS -N $JOB_NAME

# Output options: 'oe' to merge stdout/stderr into stdout, 'eo' for stderr, 'n' to not merge.
#PBS -j oe
#PBS -o $OUTPUT_FILE

# Email notifications: 'b' at begin, 'e' at end, 'a' on abort. Remove 'n' for no emails.
#PBS -m be
#PBS -M ncollier@anl.gov

# --------------- Script Execution Part --------------- #
# The following is an example setup for an MPI job.

echo "Job ID: \$PBS_JOBID"
echo "Running on host: \$(hostname)"
echo "Running on nodes: \$(cat \$PBS_NODEFILE)"

source $THIS/../envs/improv.sh

Rscript $THIS/../R/run_citycovid_experiment.R $JOB_NAME $THIS/../R/adaptiveTS_cfg.yaml
EOF