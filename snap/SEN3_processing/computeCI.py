# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:37:51 2020

Main script to call CI_computation in for loop

@author: Abhinav
"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import CI_computation
import os

from snappy import ProductIO

#path
direc='E:/EPA/lake_erie_HAB/remote_sensing_data/MERIS_RR/2011'

# read product
# list all the products with extension N1

fname_list=os.listdir(direc)
for fname in fname_list:
    direc1 = direc+'/'+fname
    if os.path.isdir(direc1):
        fname1 = os.listdir(direc1)
        CI_computation.CI_compute(direc1,fname1[0])
        
