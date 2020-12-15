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
import datetime

from snappy import ProductIO
from datetime import date

#path
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel/2016/composite_product'

# read product
# list all the products with extension dim
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

# create a list of data to be written
num = len(CI_sum)
data = []
for ind in range(0,num):
    split_txt = product_app[ind].split('_')
    begin_date = split_txt[4]
    dt = datetime.datetime.strptime(begin_date,'%Y%m%d')
    datenum = date.toordinal(dt)
    date_prop = date.fromordinal(datenum)
    date_str = str(date_prop.year)+'-'+str('%02d' %date_prop.month)+'-'+str('%02d' %date_prop.day)
    data.append([product_app[ind],CI_sum[ind],num_pixels[ind],date_str])
    data.sort(key = lambda i: i[3])
    
# write data to a text-file

filename=direc+'/'+'total_CI.txt'
fid=open(filename,'w')
fid.write('Product_name'+'\t\tTotal_CI'+'\tNumber_of_pixels_positive_CI'+'\tbegin_date\n')
for ind in range(0,num):
    fid.write(data[ind][0]+'\t\t'+str(data[ind][1])+'\t'+str(data[ind][2])+'\t'+str(data[ind][3])+'\n')
fid.close()
        

        
