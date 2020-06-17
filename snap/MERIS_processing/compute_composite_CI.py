"""
Created on Mon Jun 16 17:23:00 2020

Creates a composite image
Ref: Stumpf et al. (2012)

@author: Abhinav Gupta
"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import CI_computation
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import datetime
import math
import gpfOP
from snappy import ProductIO
from snappy import jpy
from snappy import HashMap
from snappy import GPF
from datetime import date

BandDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

#path
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/gupta353_MERIS_full_resolution_L2_2002_001_2020-05-21T14-50-31/composite_product'
fname='MER_FRS_2PPBCM_2002720_2002729.dim'
filename=direc+'/'+fname
product=ProductIO.readProduct(filename)

# read all the CI bands
bandnames=list(product.getBandNames())
CI_bandnames=[var for var in bandnames if var.split('_')[0]=='CI' and len(var.split('_'))>1 ]

# Apply bandmaths in for loop
newBandName = 'composite_CI_0'
datatype = 'float32'
expression = CI_bandnames[0]
CI_prod=gpfOP.BandMathsAG(product,newBandName,datatype,expression)

# merge
geographicError='NaN'
compProduct=gpfOP.MergeAG(product,CI_prod,geographicError)

for comp_ind in range(0,len(CI_bandnames)):
    # Apply BandMath
    newBandName = 'composite_CI_'+str(comp_ind+1)
    datatype = 'float32'
    expression = 'if (not nan(composite_CI_'+str(comp_ind)+') and not nan('+CI_bandnames[comp_ind]+')) then max(composite_CI_'+str(comp_ind)+','+CI_bandnames[comp_ind]+')'+' else (if not nan(composite_CI_'+str(comp_ind)+') then composite_CI_'+str(comp_ind)+' else '+CI_bandnames[comp_ind]+')'
    CI_prod=gpfOP.BandMathsAG(compProduct,newBandName,datatype,expression)

    # Merge
    compProduct=gpfOP.MergeAG(compProduct,CI_prod,geographicError)

# plot image    
# gpfOP.plot_image(compProduct,newBandName)

