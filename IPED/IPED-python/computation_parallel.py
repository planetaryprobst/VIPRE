#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#--------------------------------------------------------------------------
#   
#   computation_parallel.py
#
#   Algorithms for parallel / multi-core computations.
#--------------------------------------------------------------------------
#
#   Script called by submit_IPED.sh to call the matlab engine to start the 
#   IPED computations. Computes the entry conditions for all safe entry
#   opportunities of a planetary entry probe being released from an
#   interplanetary trajectory. 
#     
#   The script uses a Lookup file created by setup_parallel.py as 
#   Treference. he trajectory IDs computed are based on the output of 
#   setup_parallel.py The PBS_idx specifies the unique ID computed with 
#   this script. 
#   
#   The script calls the Matlab engine to do the computations and saves the
#   output data in text files.
#
#--------------------------------------------------------------------------
#   Call:
#   python3 computation_parallel.py 'body' hEntry 'homePath' PBS_idx
#--------------------------------------------------------------------------
#
#   -----
#   Input
#   -----
#   body      str      -          name of body
#   hEntry    (1)     km          entry altitude of probe in atmosphere 
#   homePath  str     path        home directory path to software package:
#                                 home/user/IPED 
#   PBS_idx   (1)     -           idx to identify the trajectory ID to be 
#                                 computed here
#
#   ------
#   Output
#   ------
#   txt files
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

# MATLAB
###################
import matlab.engine

# import matlabComputation as mComp 
# import matplotlib.pyplot as plt
# from matplotlib.colors import LinearSegmentedColormap
# from mpl_toolkits import mplot3d

# UTILS
###################
import src.utils.dataProcessing as dP

# RANDOM
###################
# import time
#from datetime import datetime  
# import pprint as pp

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

class vInfIDs:
    
    # Initializer / Instance Attributes
    def __init__(self,vID,body):
        
        self.vID = vID
        self.body = body
        self.vInf = lookup[vID][0][0]
        self.arrivalDate = lookup[vID][0][1]
        self.tIDs = lookup[vID][1:]
        self.data = {}
        

    # check if data file is available
    def loadData(self,hEntry,mypath):
        
        str_hEntry = dP.hEntry2str(hEntry)
        dirpath = mypath + self.body + '/' + str_hEntry + '/'
        filepath = dirpath + self.body + str_hEntry + '-vID' + str(self.vID) + '.txt'    
        
        if isfile(filepath):
            self.data[hEntry] = dP.dataPreProcessing(filepath)
        # else:
        #     self.data = False
        
        

##############################################################################


###################
# USER INPUT
###################

# body name as str or JPL ID as int
body = sys.argv[1] 
# entry altitude in km
hEntry = int(sys.argv[2])
# home directory path
homePath = sys.argv([3])
# PBS index
PBS_idx = int(sys.argv[4]) 


###################
# PATH DEFINITIONS
###################

# path from home to IPED directory
homePath = '/home/aprobst/IPED/'
# path to software package
swPath = homePath + 'IPED-python/'
# path to data folder
dataPath = swPath + 'src/data/'
# path to SPICE kernels
spicePath = homePath + 'NAIF/naif.jpl.nasa.gov/pub/naif/generic_kernels/'


###################
# INITIALIZATION
###################


str_hEntry = dP.hEntry2str(hEntry)

# Read Lookup File 
filename = swPath + 'Lookup-' + body + '.txt'    
lookup = readFile(filename)

# Read missing vIDs
filename = dataPath + body + str_hEntry + '_missingIDs.txt'
missingIDs = readFile(filename)

# vID allocation
vID = missingIDs[PBS_idx]


# Create trajecory object
vIDs = {}
vIDs[vID] = vInfIDs(vID,body)


# Start Matlab engine
eng = matlab.engine.start_matlab()

success = eng.pathDef(homePath)

# loading SPICE kernels
success = eng.loadKernels(spicePath)

if not success:
    sys.exit('Error in SPICE Kernel loading process.')    
   

# start computation for vID

# preparation of data    
v_inf = matlab.double(vIDs[vID].vInf)
epoch = vIDs[vID].arrivalDate
hEntry = float(hEntry)

datafileName = body + dP.hEntry2str(hEntry) + '-vID' + str(vID)

# Start Matlab computation
out = eng.IPED_pythonInterface(body,datafileName,v_inf,epoch,hEntry,dataPath)
            
       