# -*- coding: utf-8 -*-
"""
Created on Tue Jun 09 14:51:00 2020

computation of fraction of cloud cover 

@author: Abhinav Gupta

"""

import snappy
import numpy as np

from snappy import ProductIO

def cloud_cover(product):

    if not list(product.getBandNames()):
        return 'NaN'
    
    # compute cloud cover over the Lake Erie
    water_LE_mask = product.getBand('water_LE_mask')
    Height=product.getSceneRasterHeight()
    Width=product.getSceneRasterWidth()
    water_LE_mask_data=np.zeros(Width*Height, dtype=np.float32)
    water_LE_mask.readPixels(0,0,Width,Height,water_LE_mask_data)
    water_LE_mask_data.shape = (Height,Width)

    water_LE_sum=np.sum(water_LE_mask_data)

    # compute total number of pixels of lake Erie
    LE_mask = product.getBand('LE_mask')
    LE_mask_data=np.zeros(Width*Height, dtype=np.float32)
    LE_mask.readPixels(0,0,Width,Height,LE_mask_data)

    total_LE_pixels=np.sum(LE_mask_data)

    cloud_cover=1-(water_LE_sum/total_LE_pixels)
    return cloud_cover
