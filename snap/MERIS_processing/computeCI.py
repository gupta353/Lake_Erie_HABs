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
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/gupta353_MERIS_full_resolution_L2_2011_001_2020-05-21T00-51-46'

# read product
# list all the products with extension N1

fname_list=os.listdir(direc)
for fname in fname_list:
    if fname.endswith(".N1"):
        CI_computation.CI_compute(direc,fname)
        
