#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#--------------------------------------------------------------------------
#   IPED - Impact of the location of the Planetary Entry probe on
#   spacecraft and mission design
#
#   Setup for single-core and local computations.
#--------------------------------------------------------------------------
#
#   A software tool to compute entry conditions for all safe entry
#   opportunities of a planetary entry probe being released from an
#   interplanetary trajectory. 
#   
#   The entry conditions are computed for all safe entries available for 
#   each unique interplanetary trajectory in a trajectory set. An entry  
#   trajectory is safe if it doesn't cross a ring area. An entry trajectory 
#   is unique when the arrival date and geometry is unique in the set.  
#   Identical arrival geometries are characterized by the same unique ID.
#    
#   The script creates a Lookup file that compares existing output files  
#   with the trajectory set, identifies new trajectories and computes the  
#   missing output files to compelte the output data set by creating and 
#   adding new IDs. 
#
#   It creates an output file per unique trajectory. The file names  
#   contains the unique file ID.
#   
#   IPED is part of the software package VIPRE, consisting of VAPRE and 
#   IPED, a software to Visualize the Impact of the PRobe Entry location
#   on the spacecraft and mission design.
#   
#--------------------------------------------------------------------------
#   Call:
#   python3 setup.py 'body' hEntry prints 'homePath' 
#--------------------------------------------------------------------------
#
#   -----
#   Input
#   -----
#   body      str      -          name of body
#   hEntry    (1)     km          entry altitude of probe in atmosphere 
#   prints    (1)     -           boolean for status as screen prints 
#                                 0 for False, 1 for True
#   homePath  str     path        optional: home directory path to software
#                                 package: home/user/IPED/ 
#                                 If empty, relative paths are used.
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

def readLookupFile(filename):

    try:
        lookup = open(filename,'r')
        lookup = lookup.read()
        out = ast.literal_eval(lookup)
    except:
        out = []
    
    return out


##############################################################################

class Trajectory:
    
    # Initializer / Instance Attributes
    def __init__(self, body, tID, vID):
        
        self.body = body
        self.tID = tID
        self.traj = dct[tID]
        self.vID = vID


#############################

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
    def loadData(self,hEntry,dataPath):
        
        str_hEntry = dP.hEntry2str(hEntry)
        dirpath = dataPath + self.body + '/' + str_hEntry + '/'
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
# boolean: true for screen prints / status output
prints = bool(int(sys.argv[3]))

# home directory path
if len(sys.argv) < 5:
    homePath = '../'    
else:        
    homePath = sys.argv[4]

###################
# PATH DEFINITIONS
###################

# path to software package
swPath = homePath + 'IPED-python/'
# path to data folder
dataPath = swPath + 'src/data/'
# path to SPICE kernels
spicePath = homePath + 'NAIF/naif.jpl.nasa.gov/pub/naif/generic_kernels/'


if prints:
    print('#*************************************************************************#')
    print('#')
    print('#    Welcome to IPED')
    print('#')
    print('#*************************************************************************#')
    print('#')
    print('#    You are computing the following case:')
    print('#')
    print('#    body:           ' + body)
    print('#    entry altitude: ' + str(hEntry) + ' km')
    print('#')
    print('#*************************************************************************#')
    print('#')
    print('#    Status Output:')
    print('#')

###################
# INITIALIZATION
###################

# Does Lookup File exist?
filename = 'Lookup-' + body + '.txt'

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

# # Filtering for Jupiter Trajectories
# Jup_traj = []

# for line in data:
#     if 599 in line[14][::2]:
#         Jup_traj.append(1)
#     else:
#         Jup_traj.append(0)
        
# with open('Jup_traj_file.txt', 'w') as file:
#     for entry in Jup_traj:
#         file.write("%i\n" % entry)

# # Jup_traj_params = []  
# # Jup_traj_params.append([Jup_traj[0][2] Jup_traj[0][13] Jup_traj[0][12] 'Time of Flight, days'])

# # for line in Jup_traj[:2]:
# #     ToF = line[1] - line[0] 
# #     mass = line[13]
# #     deltaV = line[12]
# #     v_inf = math.sqrt((line[2] - line[8])**2 + (line[3] - line)
# #     Jup_traj_params.append([])
        
    
# # flyBys2 = [ dct[key][14][::2] for key in dct]
# # keys = [x for x in dct]

# # Ven = [i for i in range(len(flyBys2)) if 299 in flyBys2[i]]
# # Ear = [i for i in range(len(flyBys2)) if 399 in flyBys2[i]]
# # Mar = [i for i in range(len(flyBys2)) if 499 in flyBys2[i]]
# # Jup = [i for i in range(len(flyBys2)) if 599 in flyBys2[i]]


# creating or reading lookup file
if ~lookup:
    lookup = open(filename,'w')
    lookup.write( str(vInf) )
    lookup.close()   
    
if prints:
    print('#    Lookup file opened or created')   
    print('#    ...')
    
lookup = readLookupFile(filename)

if prints:
    print('#    Lookup file read')   
    print('#    ...')


# Create trajecory objects
trajectories = {}
vIDs = {}
missingData = []

for vID in lookup:
    
    # key = vID
    vIDs[vID] = vInfIDs(vID,body)
    # # load data into object  
    # vIDs[vID].loadData(hEntry,dataPath)        
    
    # if not bool(vIDs[vID].data):
    #     missingData.append(vID)
    
    dirPath = dataPath + body + '/' + str_hEntry + '/'
    filePath = dirPath + body + str_hEntry + '-vID' + str(vID) + '.txt'    
    
    if not isfile(filePath):
        missingData.append(vID)
            
    nIDs = len(missingData)
    # for tID in lookup[vID][1:]:
        
    #     trajectories[tID] = Trajectory(body,tID,vID)
     
if prints:    
    print('#    ' + str(nIDs) + ' missing trajectories and IDs identified')   
    print('#    ...')     
 

if bool(missingData):
    
    if prints:
        print('#    Initializing Matlab computation of missing IDs')   
        print('#    ...')    
    
    # starting Matlab engine
    eng = matlab.engine.start_matlab()
    
    if prints:
        print('#    Matlab engine started')   
        print('#    ...')      
    
    # defining Matlab paths
    success = eng.pathDef(homePath)
    
    if not success:
        sys.exit('#     ???   Matlab paths did NOT load ')    
    
    if prints:
        print('#    Matlab paths defined')   
        print('#    ...')      
    
    # loading SPICE kernels
    success = eng.loadKernels(spicePath)
    
    if not success:
        sys.exit('Error in SPICE Kernel loading process.')    
    
    if prints:
        print('#    SPICE kernels successfully loaded')   
        print('#    ...')     
        print('#*************************************************************************#')
        print('#')

    for vID in missingData:
        
        if prints:
            print('#    No ' + str(missingData.index(vID) + 1) + ' out of ' + str(nIDs) + ' - trajectory ID ' +  str(vID))
            print('#    ...')
        # preparation of data    
        v_inf = matlab.double(vIDs[vID].vInf)
        epoch = vIDs[vID].arrivalDate
        hEntry = float(hEntry)
        # tID = vIDs[vID].tIDs
        datafileName = body + dP.hEntry2str(hEntry) + '-vID' + str(vID)
        # Matlab
        out = eng.IPED_pythonInterface(body,datafileName,v_inf,epoch,hEntry,dataPath,prints)
            
else:
    if prints:
        print('#*************************************************************************#')
        print('#')
        print('#    No missing trajectories and IDs could be identified')   
        print('#')     

if prints:
    print('#    Entry data set complete.')
    print('#*************************************************************************#')
