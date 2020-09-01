# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:37:51 2020

computes SD for pixels that correspond to the Lake Erie and are cloud free
Ref: Zolfaghari and Duguay, (2016). Estimation of water-quality parameters in Lake Erie from MERIS using Mixed Effect Models  
     
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

def sd_compute(direc, fname):
    # path
    #direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/MERIS'

    # read product
    file_path=direc+'/'+fname
    product=ProductIO.readProduct(file_path)

    # read Lake Erie mask 
    mask_file='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Lake_Erie_mask/subset_0_of_Lake_Erie_mask.dim'
    mask_product=ProductIO.readProduct(mask_file)

    # create a subset of the product according to Lake Erie mask
    wkt = "POLYGON((-83.56861 41.5027, -82.08055 41.25472, -81.81277 41.9925, -83.32 42.27277, -83.56861 41.5027))"
    product=gpfOP.subsetAG(product,wkt)

    # if product is empty after subset: break
    if not list(product.getBandNames()):
        return

    # Use BandMaths operator to compute secchi depth
    newBandName = 'secchi_depth'
    datatype = 'float32'
    expression = '(LE_Mask == 1 and water) and (reflec_6>0 and reflec_4>0)? exp10(-1.04*(reflec_6/reflec_4) + 0.99): NaN'
    noDataVal='nan'
    sd_prod=gpfOP.BandMathsAG(product,newBandName,datatype,expression,noDataVal)

    # remove insane values of secchi depth
    
    # merge sd_prod_with product
    product=gpfOP.MergeAG(product,sd_prod,'NaN')

    # write product to a new file
    wfname=os.path.splitext(fname)[0]
    write_filename=direc+'/'+'sd_product/'+wfname+'.dim'
    ProductIO.writeProduct(product,write_filename,'BEAM-DIMAP')
