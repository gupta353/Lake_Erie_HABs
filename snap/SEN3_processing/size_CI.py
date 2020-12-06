"""
Created on Mon Jun 15 15:47:00 2020

filter out the imges based on height and width

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
from snappy import ProductIO
from snappy import jpy

SubsetOp = jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')


#path
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/gupta353_MERIS_full_resolution_L2_2011_001_2020-05-21T00-51-46'

# set thresholds (minimum values) for heights and widths
height_thresh=290
width_thresh=490

# exclude the products according to heights and width
fname_filter_hw=[]
fname_list=os.listdir(direc)
for fname in fname_list:
    if fname.endswith(".dim"):
        file_path=direc+'/'+fname
        product=ProductIO.readProduct(file_path)

        CI=product.getBand('CI')
        Width=CI.getRasterWidth()
        Height=CI.getRasterHeight()
        print(fname+'\t'+'Height='+str(Height),'Width='+str(Width))

        if Height<height_thresh or Width<width_thresh:
            fname_filter_hw.append(fname)

# write fname_filter to a textfile
wfname='list_of_prodcts_to_be_filtered_out.txt'
wfilename=direc+'/'+wfname
wfid=open(wfilename,'w')
for count in range(0,len(fname_filter_hw)):
    wfid.write(fname_filter_hw[count]+'\n')
wfid.close()
