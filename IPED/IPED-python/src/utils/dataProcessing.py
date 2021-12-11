#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#--------------------------------------------------------------------------
#   dataProcessing.py
#--------------------------------------------------------------------------
#
#   Collection of functions for data processing used by setup.py, 
#   setup_parallel.py and computations_parallel.py
#
#--------------------------------------------------------------------------
#
#   Use:
#   import dataProcessing as dP 
#
#*************************************************************************#
# Language: Python 3 (OSX) using Matlab 2019b
# Author: Alena Probst
# History:
# Version |    Date    |     Name      | Change history
# v1.0    | 06.12.2021 |  A. Probst    | First release
#*************************************************************************#

"""

import numpy as np

#############################

def readFile(filepath):
    
    file = open(filepath,'r')
    lines = file.readlines()
    file.close()
    
    return lines

#############################

def hEntry2str(hEntry):
    
    hEntry = str(int(hEntry))
    while len(str(hEntry)) < 5:
        hEntry = '0' + hEntry
        
    return hEntry

#############################

def dataPreProcessing(filepath):
    
    # read the file
    lines = readFile(filepath)
    # create output list
    data = []
    
    # processing the header
    header = lines[0].strip()
    header = header.split()
    data.append(header)
    
    # processing the data from file
    for line in lines[1:]:
        # removing new line sign at the end of each line
        line = line.strip()
        # split the line to a list
        line = line.split()
        
        line = [float(x) for x in line ]
        line[1] = int(line[1])
        line[2] = int(line[2])
        line[9] = int(line[9])
        
        # append line to data list
        data.append(line)
        
        
    return data

#############################

def processTrajFile(body):
    
    filepath = 'src/data/'+body+'.csv'
    lines = readFile(filepath)
    # processing the data from trajectory file
    count = 0
    data = []
    
    for line in lines:
        # removing new line sign at the end of each line
        line = line.strip()
        # split the line to a list
        line = line[0:-1].split(',')
        
        # for the header    
        if count == 0:
            # summarizing the flybys to a list
            line[-2] = line[-2:]
            # slicing the list so that there are not duplicate values
            line = line[:-1]
            # defining the number of elements in the list
            col = len(line) 
            
        # for each line except the header
        else:
            # first col-1 elements are floating point numbers
            line[:col-1] = [float(x) for x in line[:col-1]]
            # every following elements alternate between int and float
            line[col-1::2] = [int(x) for x in line[col-1::2]]
            line[col::2] = [float(x) for x in line[col::2]]
            # summarizing the number of flybys into one list element
            line[col-1] = line[col-1:]
            # slicing the list so that there are not duplicate values
            line = line[:col]
        
        # append line to data list
        data.append(line)
        count += 1
        
    return data

#############################

def createTrajID(body,data):
        
    # saving header            
    dct = {'header': data[0]}
    dct['header'].append('vHyperbolic (IAU body inertial RF)')

    for entry in data[1:]:
        # creating trajectory IDs for every entry
        ID = body + str(int(entry[0] * 100)) + \
                str(int(entry[1] * 100)) + \
                    str(entry[12])[0] + str(entry[12])[2:6]      
                    
        if ID not in dct.keys():
            dct[ID] = entry
            
            dct[ID].append((np.array(entry[2:5])).tolist())
            
    return dct
 
#############################
    
def createVinfID(dct):
    # shall create an ID for dict entries with the same vHyperbolic (dct[15]) 
    # and arrival date (dct[1])
    seen = []
    vInf = {}
    count = 0
    
    # for each key, value pair
    for key, value in dct.items():
        
        lst = [value[15],value[1]]
        
        if key != 'header':
        
            # check if value already has been seen
            if lst not in seen:
                # if not, add it to the list seen
                seen.append(lst)
                # use the current ID count as key for the new dict and safe the  
                # value in a new dict vInf, along with the key of original dct
                vInf[count] = [lst,key]
                # add vInf ID count in original dct (maybe unnecessary)
                dct[key].append(count)
                # increase count by one
                count += 1
            
            # if value has been seen already
            else:
                
                index = seen.index(lst)
                # add vInf ID count in original dct (maybe unnecessary)
                dct[key].append(index)
                      
                # look for ID count of that value and append trajectory ID to 
                # corresponding vInf entry
                vInf[index].append(key)   
            
    return vInf
