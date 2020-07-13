# -*- coding: utf-8 -*-
"""
Created on Tue Jun 09 20:40:00 2020

Main script to call cloud_cover_perc in for loop
Cloud cover is computed both over entire lake Erie mask and Eastern part of the lake

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
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/MERIS_2009_full_extent_product/CI_product'

# compute cloud cover in for loop
fname_list=os.listdir(direc)
cloud_cover=[]
for fname in fname_list:
    if fname.endswith("dim"):
        file_path=direc+'/'+fname
        product=ProductIO.readProduct(file_path)
        cloud_cover.append(fname)
        cloud_cover.append(cloud_cover_perc.cloud_cover(product))

# write cloud cover data in textfile
num=len(cloud_cover)
filename=direc+'/'+'cloud_cover.txt'
fid=open(filename,'w')
fid.write('product_name'+'\tcloud_cover_fraction\n')
for ind in range(0,num,2):
    fid.write(cloud_cover[ind]+'\t'+str(cloud_cover[ind+1])+'\n')
fid.close()

# compute cloud cover over a subset of product (Eastern part of Lake Erie)
fname_list=os.listdir(direc)
sub_cloud_cover=[]

wkt = "POLYGON((-83.59 41.51, -82.71 41.33, -82.4075 42.12, -83.308 42.305, -83.59 41.51))"
geometry = WKTReader().read(wkt)

for fname in fname_list:
        if fname.endswith("dim"):
            file_path=direc+'/'+fname
            product=ProductIO.readProduct(file_path)

            op = SubsetOp()
            op.setSourceProduct(product)
            op.setGeoRegion(geometry)
            subproduct = op.getTargetProduct()
            
            # compute cloud cover
            sub_cloud_cover.append(fname)
            sub_cloud_cover.append(cloud_cover_perc.cloud_cover(subproduct))
            
# write cloud cover data in textfile
num=len(sub_cloud_cover)
filename=direc+'/'+'subset_cloud_cover.txt'
fid=open(filename,'w')
fid.write('product_name'+'\tcloud_cover_fraction\n')
for ind in range(0,num,2):
    fid.write(sub_cloud_cover[ind]+'\t'+str(sub_cloud_cover[ind+1])+'\n')
fid.close()
