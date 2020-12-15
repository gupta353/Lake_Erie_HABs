# -*- coding: utf-8 -*-
"""
Created on Tue Jun 09 14:51:00 2020

computation of fraction of cloud cover

@author: Abhinav Gupta

"""

import snappy
import numpy as np
import gpfOP

from snappy import ProductIO
from snappy import HashMap
from snappy import GPF

# computation of cloud cover for the raw L2 product
def cloud_cover(product):

    if not list(product.getBandNames()):
        return 'NaN'
 
    ## compute cloud cover over the Lake 
    
    water_LE_mask = product.getBand('water_LE_mask')
    Height=product.getSceneRasterHeight()
    Width=product.getSceneRasterWidth()
    water_LE_mask_data=np.zeros(Width*Height, dtype=np.float32)
    water_LE_mask.readPixels(0,0,Width,Height,water_LE_mask_data)
    water_LE_mask_data.shape = (Height,Width)

    water_LE_sum=np.nansum(water_LE_mask_data)

    # compute total number of pixels of lake Erie
    LE_mask = product.getBand('LE_mask')
    LE_mask_data=np.zeros(Width*Height, dtype=np.float32)
    LE_mask.readPixels(0,0,Width,Height,LE_mask_data)

    total_LE_pixels=np.nansum(LE_mask_data)

    cloud_cover=1-(water_LE_sum/total_LE_pixels)
    return cloud_cover

# computation of cloud cover for the composite product (cloud cover over teh band 'composite_CI_final')
def cloud_cover_cmp(product):

    if not list(product.getBandNames()):
        return 'NaN'

   # apply bandmaths
    expr = 'nan(composite_CI_final)? 0:1'             # NaNs in the Lake Erie are replaced by 0
    product_cl = gpfOP.BandMathsAG(product,'cloud_cover','float32',expr,'NaN')

    cloud_cover = product_cl.getBand('cloud_cover')
    Width = cloud_cover.getRasterWidth()
    Height = cloud_cover.getRasterHeight()
    cloud_cover_data=np.zeros(Width*Height, dtype=np.float32)
    cloud_cover.readPixels(0,0,Width,Height,cloud_cover_data)
    cloud_cover_data.shape = (Height,Width)

    # number of no-cloud pixels in Lake Erie
    sum_noclcov = np.sum(cloud_cover_data)
           
    # total number of pixels in Lake Erie
    LE_mask = product.getBand('LE_Mask')
    LE_mask_data=np.zeros(Width*Height, dtype=np.float32)
    LE_mask.readPixels(0,0,Width,Height,LE_mask_data)
    total_LE_pixels=np.nansum(LE_mask_data)

    return 1-(sum_noclcov/total_LE_pixels)
     
