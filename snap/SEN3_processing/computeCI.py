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
import shutil

from snappy import ProductIO

#path
direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel/2016'


# read product
# list all the products with extension
fname_list = os.listdir(direc)
#fname_list = ['S3A_OL_2_WFR_20160822T152431_104.SEN3']
for fname in fname_list:
    if os.path.splitext(fname)[1] == '.SEN3':
        CI_computation.CI_compute(direc,fname)

        
