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
    #fname='MER_FRS_2PPBCM20110801_160652_000000313105_00112_49264_0001'
    file_path=direc+'/'+fname
    product=ProductIO.readProduct(file_path)

    # read Lake Erie mask 
    mask_file='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Lake_Erie_mask/subset_0_of_Lake_Erie_mask.dim'
    mask_product=ProductIO.readProduct(mask_file)

    # create a subset of the product according to Lake Erie mask
    wkt = "POLYGON((-83.56861 41.5027, -82.08055 41.25472, -81.81277 41.9925, -83.32 42.27277, -83.56861 41.5027))"
    geometry = WKTReader().read(wkt)

    op = SubsetOp()
    op.setSourceProduct(product)
    op.setGeoRegion(geometry)
    product = op.getTargetProduct()

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
    geometry_new = WKTReader().read(wkt_new)

    op1 = SubsetOp()
    op1.setSourceProduct(mask_product)
    op1.setGeoRegion(geometry_new)
    mask_product = op1.getTargetProduct()
    
    # Resample masked product
    param=HashMap()
    param.put('targetWidth',Width)
    param.put('targetHeight',Height)
    param.put('upsampling','Nearest')
    mask_product = GPF.createProduct('Resample', param, mask_product)

    # Merge product and mask product
    sourceProducts= HashMap()
    sourceProducts.put('masterProduct', product)
    sourceProducts.put('slaveProduct', mask_product)
    parameters = HashMap()
    parameters.put('geographicError','NaN')
    product = GPF.createProduct('Merge', parameters, sourceProducts)

    # create a band out of water mask of product
    BandDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

    targetBand = BandDescriptor()
    targetBand.name = 'water_mask'
    targetBand.type = 'float32'
    targetBand.expression = 'water? 1 : 0'
    targetBand.noDataValue=float("nan")

    targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
    targetBands[0] = targetBand

    params = HashMap()
    params.put('targetBands', targetBands)
    water_mask_prod = GPF.createProduct('BandMaths', params, product)

    # merge the product 'water_mask_prod' with product
    sourceProducts= HashMap()
    sourceProducts.put('masterProduct', product)
    sourceProducts.put('slaveProduct', water_mask_prod)
    parameters = HashMap()
    parameters.put('geographicError','NaN')
    product = GPF.createProduct('Merge', parameters, sourceProducts)

    # Use Bandmaths to create intersection of water mask and lake Erie mask
    targetBand = BandDescriptor()
    targetBand.name = 'water_LE_mask'
    targetBand.type = 'float32'
    targetBand.expression = 'water_mask*LE_Mask'
    targetBand.noDataValue=float("nan")

    targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
    targetBands[0] = targetBand

    params = HashMap()
    params.put('targetBands', targetBands)
    water_LE_prod = GPF.createProduct('BandMaths', params, product)

    # merge water_LE_prod_with product
    sourceProducts= HashMap()
    sourceProducts.put('masterProduct', product)
    sourceProducts.put('slaveProduct', water_LE_prod)
    parameters = HashMap()
    parameters.put('geographicError','NaN')
    product = GPF.createProduct('Merge', parameters, sourceProducts)

    # Use BandMaths operator to compute CI
    BandDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

    targetBand = BandDescriptor()
    targetBand.name = 'CI'
    targetBand.type = 'float32'
    targetBand.expression = 'LE_Mask == 1 and water? (reflec_7-reflec_8) + (reflec_9 - reflec_7)*(680-664)/(708-664): NaN'
    targetBand.noDataValue=float("nan")

    targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
    targetBands[0] = targetBand

    params = HashMap()
    params.put('targetBands', targetBands)
    CI_prod = GPF.createProduct('BandMaths', params, product)

    # assign 0 to all the values of CI that are less than zero
    targetBand = BandDescriptor()
    targetBand.name = 'CI'
    targetBand.type = 'float32'
    targetBand.expression = 'CI<0 ? 0: CI'
    targetBand.noDataValue=float("nan")

    targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
    targetBands[0] = targetBand

    params = HashMap()
    params.put('targetBands', targetBands)
    CI_prod = GPF.createProduct('BandMaths', params, CI_prod)

    # merge CI_prod_with product
    sourceProducts= HashMap()
    sourceProducts.put('masterProduct', product)
    sourceProducts.put('slaveProduct', CI_prod)
    parameters = HashMap()
    parameters.put('geographicError','NaN')
    product = GPF.createProduct('Merge', parameters, sourceProducts)


    # write product to a new file
    wfname=os.path.splitext(fname)[0]
    write_filename=direc+'/'+'CI_product/'+wfname+'.dim'
    ProductIO.writeProduct(product,write_filename,'BEAM-DIMAP')
