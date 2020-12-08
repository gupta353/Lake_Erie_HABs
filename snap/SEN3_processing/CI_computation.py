# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:37:51 2020

computes CI for pixels that correspond to the Lake Erie and are cloud free
Ref: Mishra and Stumpf et al. (2019).  Measurement of Cyanobacterial bloom magnitude using satellite remote sensing.
     Stumpf et al. (2012). Interannual variability of Cyanobacterial bloom in lake Erie.

@author: Abhinav Gupta

"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import geopandas as gpd
import os
import gpfOP

from snappy import ProductIO
from snappy import jpy
from snappy import GPF
from snappy import HashMap
from snappy import Mask

SubsetOp = jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')

def CI_compute(direc, fname):
    # path
    #direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/MERIS'

    # read product
    file_path=direc+'/'+fname+'/xfdumanifest.xml'
    product=ProductIO.readProduct(file_path)

    # read Lake Erie mask 
    mask_file='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Lake_Erie_mask/subset_0_of_Lake_Erie_mask.dim'
    mask_product=ProductIO.readProduct(mask_file)

    # create a subset of the product according to Lake Erie mask
    wkt = "POLYGON((-83.7117 41.3894, -81.9583 41.1292, -81.7111 42.0236, -83.4917 42.2736, -83.7117 41.3894))"
    product=gpfOP.subsetAG(product,wkt)

    # if product is empty after subset: break
    if not list(product.getBandNames()):
        return
    
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
    
    # Resample masked product
    param=HashMap()
    param.put('targetWidth',Width)
    param.put('targetHeight',Height)
    param.put('upsampling','Nearest')
    mask_product = GPF.createProduct('Resample', param, mask_product)

    # Merge product and mask product
    product=gpfOP.MergeAG(product,mask_product,'NaN')
    
    # create a band out of inland water mask of product
    newBandName = 'water_mask'
    datatype = 'float32'
    expression = 'WQSF_lsb.INLAND_WATER and WQSF_lsb.LAND ? 1 : 0'
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
    expression = 'water_LE_mask==1? (Oa08_reflectance-Oa10_reflectance) + (Oa11_reflectance - Oa08_reflectance)*(680-664)/(708-664): NaN'
    noDataVal='nan'
    CI_prod=gpfOP.BandMathsAG(product,newBandName,datatype,expression,noDataVal)

    # assign 0 to all the values of CI that are less than zero
    newBandName = 'CI'
    datatype = 'float32'
    expression = 'CI<0 ? 0: CI'
    noDataVal='nan'
    CI_prod=gpfOP.BandMathsAG(CI_prod,newBandName,datatype,expression,noDataVal)

    # merge CI_prod_with product
    product=gpfOP.MergeAG(product,CI_prod,'NaN')

    # write product to a new file
    wfname=os.path.splitext(fname)[0]
    write_filename=direc+'/'+wfname+'.dim'
    ProductIO.writeProduct(product,write_filename,'BEAM-DIMAP')
