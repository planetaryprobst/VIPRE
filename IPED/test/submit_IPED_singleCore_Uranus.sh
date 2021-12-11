#!/bin/bash
#--------------------------------------------------------------------------
#   Test case: Uranus
#   Submits IPED for single core computing using a PBS system.
#   Calls setup.py to perform the computations.
#
#--------------------------------------------------------------------------
#   Call:
#   ./submit_IPED_singleCore_Uranus.sh 'homepath'
#--------------------------------------------------------------------------
#
#   -----
#   Input
#   -----
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


#PBS -N IPED_Uranus
#PBS -q array
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=240:00:00
#PBS -m abe
#PBS -v summary=true
#PBS -M user@testemail.com
    
HOMEPATH = $1
  
### Change directory
cd $HOMEPATH/IPED-python/

### Run executable
### python 3 script.py body(str) hEntry(km)
python3 setup.py 'Uranus' 700 0 "$HOMEPATH"
