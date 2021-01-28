"""

This script computes height and width of each of the products and records them into a textfile

"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import os
import snappy

from snappy import ProductIO

direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel/2016'
fname_list = os.listdir(direc)

# read product and compute height and width of each product
prod_size = []
for fname in fname_list:
    if fname.endswith('.dim'):
        filename = direc+'/'+fname
        product = ProductIO.readProduct(filename)

        CI = product.getBand('CI')
        width = CI.getRasterWidth()
        height = CI.getRasterHeight()
        prod_size.append([fname,width,height])

# write data to a text file
sname = 'product_size.txt'
filename = direc+'/'+sname
fid = open(filename,'w+')
fid.write('product\t\tWidth\tHeight\n')
for ind in range(len(prod_size)):
    fid.write(prod_size[ind][0]+'\t\t'+str(prod_size[ind][1])+'\t'+str(prod_size[ind][2])+'\n')
fid.close()
