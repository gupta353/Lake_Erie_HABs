# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:37:51 2020

Test file to undestand SEN3 data format 

@author: Abhinav
"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import CI_computation
import os
import gpfOP
import numpy as np

from snappy import ProductIO
from snappy import jpy
from snappy import GPF
from snappy import HashMap
from snappy import Mask


#path
direc = 'E:/EPA/lake_erie_HAB/remote_sensing_data/MERIS_RR/2011'
n1 = 'EN1_MDSI_MER_RR__2P_20110505T161611_20110505T170004_048000_0141_20180818T230607_0100'
n2 = 'ENV_ME_2_RRG____20110505T161611_20110505T170004_________________2632_102_141______DSI_R_NT____.SEN3'
filename = direc+'/'+n1+'/'+n2+'/xfdumanifest.xml'

product = ProductIO.readProduct(filename)

# read Lake Erie mask 
mask_file='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Lake_Erie_mask/subset_0_of_Lake_Erie_mask.dim'
mask_product=ProductIO.readProduct(mask_file)

# create a subset of the product according to Lake Erie mask
wkt = "POLYGON((-83.7117 41.3894, -81.9583 41.1292, -81.7111 42.0236, -83.4917 42.2736, -83.7117 41.3894))"
product=gpfOP.subsetAG(product,wkt)
gpfOP.plot_image(product,'M08_rho_w')

# record widht and height of the subsetted product
Height=product.getSceneRasterHeight()
Width=product.getSceneRasterWidth()

# subset masked product to match the subsetted product
lat=product.getRasterDataNode('latitude')
lat_data=np.zeros(Width*Height, dtype=np.float32)
lat.readPixels(0,0,Width,Height,lat_data)
lat_data.shape=(Height,Width)

lat_top_left=lat_data[0,0]
lat_bottom_left=lat_data[Height-1,0]
lat_bottom_right=lat_data[Height-1,Width-1]
lat_top_right=lat_data[0,Width-1]

long=product.getRasterDataNode('longitude')
long_data=np.zeros(Width*Height, dtype=np.float32)
long.readPixels(0,0,Width,Height,long_data)
long_data.shape=(Height,Width)

long_top_left=long_data[0,0]
long_bottom_left=long_data[Height-1,0]
long_bottom_right=long_data[Height-1,Width-1]
long_top_right=long_data[0,Width-1]

wkt_new="POLYGON(("+str(long_bottom_left)+' '+str(lat_bottom_left)+', '+str(long_bottom_right)+' '+str(lat_bottom_right)+', '+str(long_top_right)+' '+str(lat_top_right)+', '+str(long_top_left)+' '+str(lat_top_left)+', '+str(long_bottom_left)+' '+str(lat_bottom_left)+"))"
mask_product=gpfOP.subsetAG(mask_product,wkt_new)
#gpfOP.plot_image(mask_product,'LE_Mask')

# Resample masked product
param=HashMap()
param.put('targetWidth',Width)
param.put('targetHeight',Height)
param.put('upsampling','Nearest')
mask_product = GPF.createProduct('Resample', param, mask_product)

# Merge product and mask product
#product=gpfOP.MergeAG(product,mask_product,'NaN')
sourceProducts= HashMap()
sourceProducts.put('masterProduct', product)
sourceProducts.put('slaveProduct', mask_product)
parameters = HashMap()
parameters.put('geographicError','NaN')
product = GPF.createProduct('Merge', parameters, sourceProducts)

# create a band out of water mask of product
newBandName = 'water_mask'
datatype = 'float32'
expression = 'CO_DO_WATER? 1 : 0'
noDataVal='nan'
water_mask_prod=gpfOP.BandMathsAG(product,newBandName,datatype,expression,noDataVal)

# merge the product 'water_mask_prod' with product
product=gpfOP.MergeAG(product,water_mask_prod,'NaN')

# Use Bandmaths to create intersection of water mask and lake Erie mask
newBandName = 'water_LE_mask'
datatype = 'float32'
expression = 'water_mask*LE_Mask'
noDataVal='nan'
water_LE_prod=gpfOP.BandMathsAG(product,newBandName,datatype,expression,noDataVal)

# merge water_LE_prod_with product
product=gpfOP.MergeAG(product,water_LE_prod,'NaN')

# Use BandMaths operator to compute CI
newBandName = 'CI'
datatype = 'float32'
expression = 'water_LE_mask==1? (M07_rho_w-M08_rho_w) + (M09_rho_w - M07_rho_w)*(680-664)/(708-664): NaN'
noDataVal='nan'
CI_prod=gpfOP.BandMathsAG(product,newBandName,datatype,expression,noDataVal)

# assign 0 to all the values of CI that are less than zero
newBandName = 'CI'
datatype = 'float32'
expression = 'CI<0 ? 0: CI'
noDataVal='nan'
CI_prod=gpfOP.BandMathsAG(CI_prod,newBandName,datatype,expression,noDataVal)
gpfOP.plot_image(CI_prod,'CI')
# merge CI_prod_with product
product=gpfOP.MergeAG(product,CI_prod,'NaN')

CI=product.getBand('CI')
Width=CI.getRasterWidth()
Height=CI.getRasterHeight()
CI_data = np.zeros(Width*Height, dtype=np.float32)
CI.readPixels(0,0,Width,Height,CI_data)
CI_data.shape = (Height, Width)

CI_sum2 = np.nansum(CI_data)
