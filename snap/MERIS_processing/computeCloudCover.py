# -*- coding: utf-8 -*-
"""
Created on Tue Jun 09 20:40:00 2020

Main script to call cloud_cover_perc in for loop
Cloud cover is computed both over entire lake Erie mask and Western part of the lake
This script is written to compute cloud cover in composite products

@author: Abhinav
"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import cloud_cover_perc
import os
import datetime

from datetime import datetime
from datetime import date
from snappy import ProductIO
from snappy import jpy

SubsetOp = jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')

#path
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/gupta353_MERIS_full_resolution_L2_2011_001_2020-05-21T00-51-46/composite_product'

# compute cloud cover in for loop
fname_list=os.listdir(direc)
cloud_cover=[]
for fname in fname_list:
    if fname.endswith("dim"):
        file_path=direc+'/'+fname
        product=ProductIO.readProduct(file_path)
        split_text = fname.split('_')
        begin_date = split_text[3]
        dt = datetime.strptime(begin_date,'%Y%m%d')
        datenum = date.toordinal(dt)
        date_prop = date.fromordinal(datenum)
        date_str = str(date_prop.year)+'-'+str('%02d' %date_prop.month)+'-'+str('%02d' %date_prop.day)
        cloud_cover.append([fname,cloud_cover_perc.cloud_cover_cmp(product),date_str,datenum])
        cloud_cover.sort(key = lambda i: i[3])
        
# write cloud cover data in textfile
num=len(cloud_cover)
filename=direc+'/'+'cloud_cover.txt'
fid=open(filename,'w')
fid.write('product_name'+'\tcloud_cover_fraction'+'\tbegin_date\n')
for ind in range(0,num):
    fid.write(cloud_cover[ind][0]+'\t'+str(cloud_cover[ind][1])+'\t'+str(cloud_cover[ind][2])+'\n')
fid.close()

"""
# compute cloud cover over a subset of product (Western part of Lake Erie)
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
            split_text = fname.split('_')
            begin_date = split_text[3]
            dt = datetime.strptime(begin_date,'%Y%m%d')
            datenum = date.toordinal(dt)
            date_prop = date.fromordinal(datenum)
            date_str = str(date_prop.year)+'-'+str('%02d' %date_prop.month)+'-'+str('%02d' %date_prop.day)
            
            cloud_cover.append([fname,cloud_cover_perc.cloud_cover(product),date_str,datenum])
            cloud_cover.sort(key = lambda i: i[3])
            sub_cloud_cover.append([fname,cloud_cover_perc.cloud_cover_cmp(subproduct),date_str,datenum])
            sub_cloud_cover.sort(key = lambda i: i[3])
            
# write cloud cover data in textfile
num=len(sub_cloud_cover)
filename=direc+'/'+'subset_cloud_cover.txt'
fid=open(filename,'w')
fid.write('product_name'+'\tcloud_cover_fraction'+'\tbegin_date\n')
for ind in range(0,num):
    fid.write(sub_cloud_cover[ind][0]+'\t'+str(sub_cloud_cover[ind][1])+'\t'+str(sub_cloud_cover[ind][2])+'\n')
fid.close()
"""
