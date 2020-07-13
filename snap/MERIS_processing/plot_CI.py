# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:37:51 2020

plot of a product band

@author: Abhinav
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

#path
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/gupta353_MERIS_full_resolution_L2_2011_001_2020-05-21T00-51-46'

# read product
# list all the products with extension N1
fname_list=os.listdir(direc)
for fname in fname_list:
    if fname.endswith(".dim"):
        file_path=direc+'/'+fname
        product=ProductIO.readProduct(file_path)
        CI=product.getBand('CI')
        Width=CI.getRasterWidth()
        Height=CI.getRasterHeight()
        CI_data = np.zeros(Width*Height, dtype=np.float32)
        CI.readPixels(0,0,Width,Height,CI_data)
        CI_data.shape = (Height, Width)

        plt.figure(figsize=(8, 8)) # adjusting the figure window size
        fig = plt.imshow(CI_data, cmap = cm.jet) #matplotlib settings for the current image
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        plt.colorbar(fraction=0.05,shrink=0.6)

        name=fname.split('.')
        save_name=name[0]+'.png'
        savefile=direc+'/'+save_name
        plt.savefig(savefile,dpi=300,quality=100)
        plt.close()
