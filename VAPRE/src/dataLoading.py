#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#--------------------------------------------------------------------------
#   dataLoading.py
#--------------------------------------------------------------------------
#
#   Collection of functions that load the IPED data input files and make 
#   them available for processing.
#
#--------------------------------------------------------------------------
#
#   Use:
#   import src.dataLoading as dL 
#
#*************************************************************************#
# Language: Python 3 (OSX) using Matlab 2019b
# Author: Alena Probst
# History:
# Version |    Date    |     Name      | Change history
# v1.0    | 06.12.2021 |  A. Probst    | First release
#*************************************************************************#
"""

# UTILS
###################

import src.utils.dataFunctions as dF # former dP or aDP
import src.dataProcessing as dP


def load_trajectoryFile(body):
    return dF.processTrajFile(body)

def load_entryFile(body,vID,hEntry):
    # reading body data file
    return dF.readDataFile(vID,body,hEntry)

def load_trajectoryData(body):
    
    # reading trajectory file data
    trajFile = dF.processTrajFile(body)
    # generate trajectory IDs
    dct,vInf = dP.generate_trajectoryIDs(body,trajFile)
    # generate trajectory data
    launchDate,ToF,vHyp,arrMass,deltaV,vIDs,tIDs = dP.generate_trajectoryData(trajFile,dct)
    # pack variables in list   
    data = [dct,vInf,launchDate,ToF,vHyp,arrMass,deltaV,vIDs,tIDs]
    
    return data