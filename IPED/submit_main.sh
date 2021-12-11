#!/bin/bash
#--------------------------------------------------------------------------
#   Bash script to submit IPED for parallel computing using a PBS system.
#   Calls setup_parallel.py to prepare computations, and relays the output
#   to the submit_IPED.sh file for parallel computing.
#
#   Takes as input the body name, the entry altitude and the homepath of
#   the software directory IPED.
#
#--------------------------------------------------------------------------
#   Call:
#   ./submit_main.sh 'body' hEntry 'homepath'
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
#   -
#
#*************************************************************************#
# Language: MATLAB R2019b (OSX)
# Author: Alena Probst
# History:
# Version |    Date    |     Name      | Change history
# v1.0    | 08.12.2021 |  A. Probst    | First release
#*************************************************************************#

### Set up environment variables
BODY=$1
HENTRY=$2
HOMEPATH = $3

### Change directory
cd $HOMEPATH/IPED-python/

### Run setup
### python 3 setup_parallel.py body(str) hEntry(km) homepath(str)
LEN=`python3 setup_parallel.py $BODY $HENTRY $HOMEPATH`

if [ $LEN -gt 0 ]
then 
    LEN=$(expr $LEN - 1)
fi

cd $HOMEPATH

qsub -J 0-$LEN -N IPED_$BODY -v BODY=$BODY,HENTRY=$HENTRY,HOMEPATH=$HOMEPATH submit_IPED.sh

### Remove Missing Data File
#HENTRY=$( printf '%05d' $HENTRY )
#rm ./src/data/$BODY$HENTRY"_missingIDs.txt"
