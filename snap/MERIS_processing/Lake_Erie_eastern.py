# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 21:03:00 2020

Eastern part of lake

@author: Abhinav
"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import cloud_cover_perc
import os

from snappy import ProductIO
from snappy import jpy

SubsetOp = jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')

#path
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Lake_Erie_mask'
fname='subset_0_of_Lake_Erie_mask.dim'
filename=direc+'/'+fname

product=ProductIO.readProduct(filename)

# create subset
wkt = "POLYGON((-83.59 41.51, -82.71 41.33, -82.4075 42.12, -83.308 42.305, -83.59 41.51))"
geometry = WKTReader().read(wkt)

op = SubsetOp()
op.setSourceProduct(product)
op.setGeoRegion(geometry)
subproduct = op.getTargetProduct()

# write the subsetted product
wfname=os.path.splitext(fname)[0]+'_eastern'
write_filename=direc+'/'+wfname+'.dim'
ProductIO.writeProduct(subproduct,write_filename,'BEAM-DIMAP')
