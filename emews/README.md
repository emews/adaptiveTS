# EMEWS Runs #

Prior to any runs the environment needs to be loaded with

```bash
$ source /lcrc/project/EMEWS/improv/envs/improv_swift_11062024.sh
```

## Sweep ##

The sweep workflow reads a csv file (a *UPF* file), creating a command line to pass to `run_experiment.R` for each
row in the csv file. See `data/upfs/config.csv` for an example.

To run a sweep, use the `swift/improv_run_sweep.sh` script. This script takes two arguments:

1. A sweep experiment id. This will be used to create directory in `emews/experiments`
in which the experiment will be run.
2. A configuration file. This file will be used to configure the workflow.

For example,

```bash
$ cd swift
$ ./improv_run_sweep.sh experiment_1 ../data/cfgs/improv_sweep.cfg
```

This will create an experiment_1 directory in `emews/experiments` and run a sweep workflow using
the configuration in `emews/data/cfgs/improv_sweep.cfg`.

The workflow configuration file has the following format:

```bash
CFG_WALLTIME=00:20:00
NODES=2
CFG_PPN=64
CFG_PROCS=$(( NODES * CFG_PPN ))
# SHOULD BE ~ 128 / CFG_PPN
CFG_MKL_THREADS=2

CFG_QUEUE=compute
CFG_PROJECT=EMEWS

CFG_UPF_FILE=upfs/upf_100.csv
```

* CFG_WALLTIME - the workflow job's walltime
* NODES - the number of nodes to allocate to the job
* CFG_PPN - the number of procs per node. The maximum is 128 on Improv. 
* CFG_PROCS - the total number of processes allocated to the job. This line should not be edited.
* CFG_MKL_THREADS - the number of threads to allocate to each R instance. See below for more details.
* CFG_QUEUE - the machine queue to run
* CFG_PROJECT - the project to charge for the run
* CFG_UPF_FILE - the sweep csv file to use

The results of a run will be written the experiment directory's `results` directory.
The file name will be prefixed with the `exp_id` value from the csv upf file.

### PPN ###

Setting the CFG_PPN value, sets the number of R instances (how many rows in the UPF file) to run on each node. The workflow itself
requires 2, so this must always be at least 3. There's probably some sweet spot here that makes best use of
threads while not using too many nodes. For now, maybe `CFG_PPN` of 64 and `CFG_MKL_THREADS`=2 is a good starting point.

### MKL THREADS ###

The R distribution on Improv uses Intel's MKL library which in turn uses parallel
threads to speed up computation. By default, each R instance run on a node in our workflow sees
all the available threads on a node regardless of how many other R instances are running. So,
on Improv each R instance will see 128 threads and try to use all these resulting in thread contention
and much slower runtimes. The `CFG_MKL_THREADS` value is exported as `MKL_NUM_THREADS` which tells
the MKL library how many threads to use rather than the default 128. On Improv this value should be
128 / CFG_PPN where CFG_PPN will set the number of R instances run on each node.