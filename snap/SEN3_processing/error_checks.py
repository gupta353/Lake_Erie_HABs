"""
error checks

"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import os
import snappy

from snappy import ProductIO

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel/2016/composite_product'
fname_list = os.listdir(direc)

## check size of the CI band used to create composite product
#fname_list = ['S3A_OL_2_WFR_201689_2016818.dim']
for fname in fname_list:
    if fname.endswith('.dim'):
        filename = direc+'/'+fname
        product = ProductIO.readProduct(filename)

        bands = list(product.getBandNames())
        CI_bands = [x for x in bands if 'CI_73' in x]

        for bnd in CI_bands:
            CI = product.getBand(bnd)
            width = CI.getRasterWidth()
            height = CI.getRasterHeight()
            print("Widht: "+str(width)+'\t'+"Height: "+str(height))

