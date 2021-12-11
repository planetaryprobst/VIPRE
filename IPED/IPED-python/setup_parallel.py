#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#--------------------------------------------------------------------------
#
#   setup_parallel.py
#
#   Setup for parallel / multi-core computations.
#--------------------------------------------------------------------------
#
#   Script called by submit_main.sh to define the number of missing IDs that 
#   still need to be computed.
#
#   The script creates a Lookup file that compares existing output files with 
#   the trajectory set and identifies new trajectories. 
#
#   The missing IDs are handed over to the script computation_parallel.py to 
#   create an output file per unique trajectory.
#
#--------------------------------------------------------------------------
#   Call:
#   python3 setup_parallel.py 'body' hEntry 'homePath' 
#--------------------------------------------------------------------------
#
#   -----
#   Input
#   -----
#   body      str      -          name of body
#   hEntry    (1)     km          entry altitude of probe in atmosphere 
#   homePath  str     path        home directory path to software package:
#                                 home/user/IPED 
#
#   ------
#   Output
#   ------
#   missingIDs  lst     -        trajectory IDs as strings in a list
#  
#*************************************************************************#
# Language: Python 3 (OSX) using Matlab 2019b
# Author: Alena Probst
# History:
# Version |    Date    |     Name      | Change history
# v1.0    | 06.12.2021 |  A. Probst    | First release
#*************************************************************************#

"""


###################
# IMPORT MODULES
###################

# UTILS
###################
import src.utils.dataProcessing as dP


# DATA PROCESSING
###################
import sys
# from os import listdir, walk
from os.path import isfile
import ast

##############################################################################

def readFile(filename):

    try:
        file = open(filename,'r')
        file = file.read()
        out = ast.literal_eval(file)
    except:
        out = []
    
    return out
##############################################################################

###################
# USER INPUT
###################

# body name as str or JPL ID as int
body = sys.argv[1] 
# entry altitude in km
hEntry = int(sys.argv[2]) 
# home directory path
homePath = sys.argv[3]


###################
# PATH DEFINITIONS
###################

# path to software package
swPath = homePath + 'IPED-python/'
# path to data folder
dataPath = swPath + 'src/data/'


###################
# INITIALIZATION
###################
    

# Does Lookup File exist?
filename = swPath + 'Lookup-' + body + '.txt'

if isfile(filename):
    lookup = True
else:
    lookup = False

# reading traj data file
data = dP.processTrajFile(body)

# creating dict, with the pre-processed data (key = ID) incl. vHyp
dct = dP.createTrajID(body,data)

# creating dict, with the vID and vInf allocation based on dct
vInf = dP.createVinfID(dct)

str_hEntry = dP.hEntry2str(hEntry)

# creating or reading lookup file
if ~lookup:
    lookup = open(filename,'w')
    lookup.write( str(vInf) )
    lookup.close()   
    
lookup = readFile(filename)


# Check which vIDs still need to be calculated
missingIDs = []

for vID in lookup:
     
    dirpath = dataPath + body + '/' + str_hEntry + '/'
    filepath = dirpath + body + str_hEntry + '-vID' + str(vID) + '.txt'    
    
    if not isfile(filepath):
        missingIDs.append(vID)


if bool(missingIDs):        
    filename = dataPath + body + str_hEntry + '_missingIDs.txt'
    missingIDsFile = open(filename,'w')
    missingIDsFile.write( str(missingIDs) )
    missingIDsFile.close() 
    
print(len(missingIDs))