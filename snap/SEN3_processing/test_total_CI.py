"""
This script computes teh sum of CI for two Full and reduced resolution MERIS images

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

# Reduced resolution
#direc = 'E:/EPA/lake_erie_HAB/remote_sensing_data/MERIS_RR/2011/EN1_MDSI_MER_RR__2P_20110509T152910_20110509T161303_048057_0198_20180819T010117_0100'
#fname = 'ENV_ME_2_RRG____20110509T152910_20110509T161303_________________2632_102_198______DSI_R_NT____.dim'
#filename = direc+'/'+fname
#product = ProductIO.readProduct(filename)

# store data in an array
#CI=product.getBand('CI')
#Width=CI.getRasterWidth()
#Height=CI.getRasterHeight()
#CI_data = np.zeros(Width*Height, dtype=np.float32)
#CI.readPixels(0,0,Width,Height,CI_data)
#CI_data.shape = (Height, Width)

#CI_sum1 = np.nansum(CI_data)

# full_resolution
direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/gupta353_MERIS_full_resolution_L2_2011_001_2020-05-21T00-51-46'
fname = 'MER_FRS_2PPBCM20110806_162322_000000343105_00184_49336_0068.dim'
filename = direc+'/'+fname
product = ProductIO.readProduct(filename)

# store data in an array
CI=product.getBand('CI')
Width=CI.getRasterWidth()
Height=CI.getRasterHeight()
CI_data = np.zeros(Width*Height, dtype=np.float32)
CI.readPixels(0,0,Width,Height,CI_data)
CI_data.shape = (Height, Width)

CI_sum2 = np.nansum(CI_data)
