# -*- coding: utf-8 -*-
"""
Created on Mon July 31 19:11:00 2020

Computation of average secchi depth over all pixels using composite secchi depth images as inputs

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
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/gupta353_MERIS_full_resolution_L2_2006_001_2020-05-19T22-11-38/composite_sd_product'

# read product
# list all the products with extension dim
fname_list=os.listdir(direc)
product_app = []
SD_sum = []                 # sum of SD
num_pixels = []             # Total number of pixels
for fname in fname_list:
    if fname.endswith(".dim"):
        file_path=direc+'/'+fname
        product=ProductIO.readProduct(file_path)
        SD=product.getBand('composite_secchi_depth_final')
        Width=SD.getRasterWidth()
        Height=SD.getRasterHeight()
        SD_data = np.zeros(Width*Height, dtype=np.float32)
        SD.readPixels(0,0,Width,Height,SD_data)
        SD_data.shape = (Height, Width)

        # compute number of pixels with positive SD
        SD_data_tmp = copy.deepcopy(SD_data)
        nan_ind = np.isnan(SD_data_tmp)
        SD_data_tmp[nan_ind] = -1
        inds = np.argwhere(SD_data_tmp>=0)
        num_pixels.append(len(inds))

        # compute the sum of SD over all pixels
        product_app.append(fname)
        SD_sum.append(np.nansum(SD_data))

# create a list of data to be written
num = len(SD_sum)
data = []
for ind in range(0,num):
    split_txt = product_app[ind].split('_')
    begin_date = split_txt[3]
    dt = datetime.datetime.strptime(begin_date,'%Y%m%d')
    datenum = date.toordinal(dt)
    date_prop = date.fromordinal(datenum)
    date_str = str(date_prop.year)+'-'+str('%02d' %date_prop.month)+'-'+str('%02d' %date_prop.day)
    data.append([product_app[ind],SD_sum[ind],num_pixels[ind],date_str,SD_sum[ind]/num_pixels[ind]])
    data.sort(key = lambda i: i[3])
    
# write data to a text-file
filename=direc+'/'+'total_secchi_depth.txt'
fid=open(filename,'w')
fid.write('Product_name'+'\tTotal_secchi_depth(m)'+'\tNumber_of_pixels_positive_SD'+'\tbegin_date\n'+'\taverage_secchi_depth(m)')
for ind in range(0,num):
    fid.write(data[ind][0]+'\t'+str(data[ind][1])+'\t'+str(data[ind][2])+'\t'+str(data[ind][3])+'\t'+str(data[ind][4])+'\n')
fid.close()
        

        
