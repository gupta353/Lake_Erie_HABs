# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:37:51 2020

Main script to call secchi_depth_computation in for loop

@author: Abhinav
"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import secchi_depth_computation
import os

from snappy import ProductIO

#path
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel/2020'

# read product
# list all the products with extension '.dim'

fname_list=os.listdir(direc)
for fname in fname_list:
    if fname.endswith(".dim"):
        secchi_depth_computation.sd_compute(direc,fname)
        
