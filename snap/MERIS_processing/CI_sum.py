# -*- coding: utf-8 -*-
"""
Created on Wed July 13 15:31:00 2020

Computation of total CI using composite_CI images as inputs

@author: Abhinav
"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import CI_computation
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy

from snappy import ProductIO

#path
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/gupta353_MERIS_full_resolution_L2_2011_001_2020-05-21T00-51-46/composite_product'

# read product
# list all the products with extension N1
fname_list=os.listdir(direc)
product_app = []
CI_sum = []                 # sum of CI
num_pixels = []             # Total number of pixels
for fname in fname_list:
    if fname.endswith(".dim"):
        file_path=direc+'/'+fname
        product=ProductIO.readProduct(file_path)
        CI=product.getBand('composite_CI_final')
        Width=CI.getRasterWidth()
        Height=CI.getRasterHeight()
        CI_data = np.zeros(Width*Height, dtype=np.float32)
        CI.readPixels(0,0,Width,Height,CI_data)
        CI_data.shape = (Height, Width)

        # compute number of pixels with positive CI
        CI_data_tmp = copy.deepcopy(CI_data)
        nan_ind = np.isnan(CI_data_tmp)
        CI_data_tmp[nan_ind] = -1
        inds = np.argwhere(CI_data_tmp>0)
        num_pixels.append(len(inds))

        # compute the sum of CI over all pixels
        product_app.append(fname)
        CI_sum.append(np.nansum(CI_data))
        
# write data to a text-file
num = len(CI_sum)
filename=direc+'/'+'total_CI.txt'
fid=open(filename,'w')
fid.write('Product_name'+'\tTotal_CI'+'\tNumber_of_pixels_positive_CI\n')
for ind in range(0,num):
    fid.write(product_app[ind]+'\t'+str(CI_sum[ind])+'\t'+str(num_pixels[ind])+'\n')
fid.close()
        

        
