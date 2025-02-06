#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=32
#PBS -k oe
#PBS -N goad
#PBS -q core32

# An example script to run GOAD on a PBS cluster.

# echo ------------------------------------------------------
# echo -n 'GOAD is running on node '; cat $PBS_NODEFILE
# echo ------------------------------------------------------
# echo PBS: qsub is running on $PBS_O_HOST
# echo PBS: originating queue is $PBS_O_QUEUE
# echo PBS: executing queue is $PBS_QUEUE
# echo PBS: working directory is $PBS_O_WORKDIR
# echo PBS: execution mode is $PBS_ENVIRONMENT
# echo PBS: job identifier is $PBS_JOBID
# echo PBS: job name is $PBS_JOBNAME
# echo PBS: node file is $PBS_NODEFILE
# echo PBS: current home directory is $PBS_O_HOME
# echo PBS: PATH = $PBS_O_PATH
# echo ------------------------------------------------------

# Set the path to the GOAD directory
GOAD_DIR=/home/username/goad 

# Change working directory to the GOAD directory for access to the config files
cd $GOAD_DIR

$GOAD_DIR/target/release/goad \
-w 0.532 \
--geo ./examples/data/hex.obj \
> log
