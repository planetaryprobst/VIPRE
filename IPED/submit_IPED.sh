#!/bin/bash
#--------------------------------------------------------------------------
#   Bash script to submit IPED for parallel computing using a PBS system.
#   Receives input from submit_main.sh. Submits computation_parallel.py per
#   missing trajectory ID.
#
#   Takes as input the body name, the entry altitude and the homepath of
#   the software directory IPED.
#
#--------------------------------------------------------------------------
#   Call:
#   BODY=$BODY,HENTRY=$HENTRY,HOMEPATH=$HOMEPATH submit_IPED.sh
#--------------------------------------------------------------------------
#
#   -----
#   Input
#   -----
#   body        str         body name
#   hEntry      (1)         entry altitude
#   hompepath   str         path to software directory, e.g. home/user/IPED
#
#   ------
#   Output
#   ------
#
#   text files
#
#*************************************************************************#
# Language: MATLAB R2019b (OSX)
# Author: Alena Probst
# History:
# Version |    Date    |     Name      | Change history
# v1.0    | 08.12.2021 |  A. Probst    | First release
#*************************************************************************#

#PBS -q array
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=12:00:00
#PBS -m abe
#PBS -v summary=true
#PBS -M user@testemail.com


### Load modules into your environment
module load matlab/2019a
module load python36
 
### Run code
### NOTE: The PBS_ARRAY_INDEX environment variable is substituted
###       with the job array index, whose range is defined with
###       the "#PBS -J 1-1000" directive above.


### Change directory
cd $HOMEPATH/IPED-python/

### Run executable
### python 3 computation_parallel.py body(str) hEntry(km) homepath(str) PBS_idx
python3 computation_parallel.py "$BODY" "$HENTRY" "$HOMEPATH" "$PBS_ARRAY_INDEX"
