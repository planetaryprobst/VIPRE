#!/bin/bash
#--------------------------------------------------------------------------
#   Bash script to submit IPED for single core computing using a PBS
#   system. Calls setup.py to perform the computations.
#
#   Takes as input the body name, the entry altitude and the homepath of
#   the software directory IPED.
#
#--------------------------------------------------------------------------
#   Call:
#   ./submit_IPED_singleCore.sh 'body' hEntry prints 'homepath'
#--------------------------------------------------------------------------
#
#   -----
#   Input
#   -----
#   body        str         body name
#   hEntry      (1)         entry altitude
#   prints      (1)         1 for screen prints, 0 for no outputs
#   hompepath   str         path to software directory, e.g. home/user/IPED
#
#   ------
#   Output
#   ------
#
#   -
#
#*************************************************************************#
# Language: MATLAB R2019b (OSX)
# Author: Alena Probst
# History:
# Version |    Date    |     Name      | Change history
# v1.0    | 08.12.2021 |  A. Probst    | First release
#*************************************************************************#


#PBS -N IPED_Neptune
#PBS -q array
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=240:00:00
#PBS -m abe
#PBS -v summary=true
#PBS -M user@testemail.com

BODY = $1
HENTRY = $2
PRINTS = $3
HOMEPATH = $4
  
### Change directory
cd $HOMEPATH/IPED-python/

### Run executable
### python 3 script.py body(str) hEntry(km)
python3 setup.py "$BODY" "$HENTRY$ "$PRINTS" "$HOMEPATH"
