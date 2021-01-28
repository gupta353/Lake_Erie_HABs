"""
Created on Thus Aug 27 15:09:00 2020

Merge all the products that fall into a time-window and save them to create a composite image of secchi disk depth 

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
import compute_composite_CI
import re
import gpfOP

from snappy import ProductIO
from snappy import jpy
from snappy import HashMap
from snappy import GPF
from datetime import date
from compute_composite_SDD import computeCompositeSDD

SubsetOp = jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')
NodeDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.MergeOp$NodeDescriptor')

#path
direc='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel/2020'
# target widht and target heights of resamplig
targetWidth=500
targetHeight=300

# read the list of products to be filtered out
filename = direc+'/list_of_products_to_be_filtered_out.txt'
fid = open(filename,'r')
filter_prods = fid.read().splitlines()
fid.close()

# time-window period to create a composite image
time_window_per=10

# being and end dates
begin_date='2020-05-01'
end_date='2020-10-31'

# list the min-max dates of each time-window including begin_date
begin_date_obj=datetime.datetime.strptime(begin_date,'%Y-%m-%d')
begin_datenum=date.toordinal(date(begin_date_obj.year,begin_date_obj.month,begin_date_obj.day))
end_date_obj=datetime.datetime.strptime(end_date,'%Y-%m-%d')
end_datenum=date.toordinal(date(end_date_obj.year,end_date_obj.month,end_date_obj.day))

num_time_windows=math.floor((end_datenum-begin_datenum+1)/time_window_per)  # total number of time-windows

time_window=[]
for ind in range(0,num_time_windows):
    time_window.append([begin_datenum+ind*time_window_per,begin_datenum+(ind+1)*time_window_per-1])
    

# list the products and remove the products listed in filter_prods
fname_list=os.listdir(direc)
for prod in filter_prods:
    fname_list.remove(prod)    

# remove all the products with extension other than '.dim'
fname_list = [fname for fname in fname_list if fname.endswith('.dim')]

# list the dates of each product
product_date=[]
product_datenum=[]
product_engdate=[]
for fname in fname_list:
    if fname.endswith("dim"):
        x=fname.split('_')[4]
        y=datetime.datetime.strptime(x.split('T')[0],'%Y%m%d')
        z=date.toordinal(y)
        product_datenum.append(z)
        product_date.append([fname,z])
        #product_engdate.append(x.split('2PPBCM')[1])

# Create a product containing sd band (filled with) with maximum size out of the products to be listed in the same time-window
for tw_ind in range(0,num_time_windows):
    ind=[i for i, datenum_value in enumerate(product_datenum) if(datenum_value>=time_window[tw_ind][0] and datenum_value<=time_window[tw_ind][1])] # find the indices of the product that are contained in this time-window

    if len(ind) != 0:
        widths = []
        heights = []
        fname_tw = []
        for i in ind:
            fname = fname_list[i]
            fname_tw.append(fname)
            filename = direc+'/'+fname
            product = ProductIO.readProduct(filename)
            CI = product.getBand('secchi_depth_over_pos_CI')
            widths.append(CI.getRasterWidth())
            heights.append(CI.getRasterHeight())

            max_width = max(widths)
            max_height = max(heights)

            i = [i for i in range(len(widths)) if widths[i]==max_width]

            # create a NaN band in the product with largest width
            fname = fname_tw[i[0]]
            filename = direc+'/'+fname
            product = ProductIO.readProduct(filename)
            newBandName = 'secchi_depth_over_pos_CI'
            datatype = 'float32'
            expression = 'NaN'
            noDataVal='nan'
            composite_prod = gpfOP.BandMathsAG(product,newBandName,datatype,expression,noDataVal);

            # create a Lake_erie_water band
            newBandName = 'secchi_depth_over_pos_CI'
            datatype = 'float32'
            expression = 'WQSF_lsb_INLAND_WATER?1:0'
            noDataVal='nan'
            prod_tmp = gpfOP.BandMathsAG(product,newBandName,datatype,expression,noDataVal);

            #composite_prod = gpfOP.MergeAG(composite_prod,prod_tmp,'NaN')
            
        ## mosaic composite_prod with other products in the time-window
        #  mosaic parameters
        westBound = -83.60220985614447
        eastBound = -81.78879308692983
        northBound = 42.31086488661365
        southBound = 41.1893175594222
        pixelSizeY = 0.00270
        pixelSizeX = 0.00357

        # mosaic water band
        varName = 'LE_mask'
        varExp = 'secchi_depth_over_pos_CI'
        product1 = ProductIO.readProduct(filename)
        LE_mask_prd = gpfOP.mosaicAG(composite_prod,prod_tmp,westBound,eastBound,northBound,southBound,pixelSizeX,pixelSizeY,varName,varExp)

        count = 0
        for fname in fname_tw:
            count = count + 1
            filename = direc+'/'+fname
            varName = 'secchi_depth_over_pos_CI_'+str(count)
            varExp = 'secchi_depth_over_pos_CI'
            product1 = ProductIO.readProduct(filename)
            prd_tmp = gpfOP.mosaicAG(composite_prod,product1,westBound,eastBound,northBound,southBound,pixelSizeX,pixelSizeY,varName,varExp)

            newBandName = 'secchi_depth_over_pos_CI_'+str(count)
            datatype = 'float32'
            expression = 'secchi_depth_over_pos_CI_'+str(count)+'_count ==1?secchi_depth_over_pos_CI_'+str(count)+':NaN'
            noDataVal = 'NaN'
            prd_tmp = gpfOP.BandMathsAG(prd_tmp,newBandName,datatype,expression,noDataVal)

            composite_prod = gpfOP.MergeAG(composite_prod,prd_tmp,'NaN')

        # compute composite image
        composite_prod = gpfOP.MergeAG(composite_prod,LE_mask_prd,'NaN')
        compositeProduct=computeCompositeSDD(composite_prod)
                
        # write the composite product
        bdate_obj=datetime.date.fromordinal(time_window[tw_ind][0])
        bdate=str(bdate_obj.year)+str(bdate_obj.month)+str(bdate_obj.day)
        edate_obj=datetime.date.fromordinal(time_window[tw_ind][1])
        edate=str(edate_obj.year)+str(edate_obj.month)+str(edate_obj.day)

        wfname='S3A_OL_2_WFR_'+bdate+'_'+edate+'.dim'
        wfilename=direc+'/composite_sd_product/'+wfname
        ProductIO.writeProduct(compositeProduct,wfilename,'BEAM-DIMAP')

"""
# resample, merge products and compute composite image in each time-window
for tw_ind in range(0,num_time_windows):
    ind=[i for i, datenum_value in enumerate(product_datenum) if(datenum_value>=time_window[tw_ind][0] and datenum_value<=time_window[tw_ind][1])] # find the indices of the product that are contained in this time-window
    
    if len(ind) != 0:
        fname=product_date[ind[0]][0]
        filename=direc+'/'+fname
        composite_prod=ProductIO.readProduct(filename)
        
        lst_bnd=ind;
        count=-1
        for prod_ind in ind:
            count=count+1
            fname=product_date[prod_ind][0]
            filename=direc+'/'+fname
            prod_tmp=ProductIO.readProduct(filename)

            # resample the product
            param=HashMap()
            param.put('targetWidth',targetWidth)
            param.put('targetHeight',targetHeight)
            param.put('upsampling','Nearest')
            param.put('downsampling','Mean')
            prod_tmp = GPF.createProduct('Resample', param, prod_tmp)

            # Merge products
            include_2 = NodeDescriptor()
            include_2.setProductId('slaveProduct')
            include_2.setName('secchi_depth_over_pos_CI')
            lst_bnd[count]='secchi_depth_over_pos_CI_'+str(product_engdate[prod_ind])
            include_2.setNewName(lst_bnd[count])
                
            include_bands = jpy.array('org.esa.snap.core.gpf.common.MergeOp$NodeDescriptor',1)
            include_bands[0]=include_2

            sourceProducts= HashMap()
            sourceProducts.put('masterProduct', composite_prod)
            sourceProducts.put('slaveProduct', prod_tmp)
            parameters = HashMap()
            parameters.put('geographicError','NaN')
            parameters.put('includes',include_bands)
            composite_prod = GPF.createProduct('Merge', parameters, sourceProducts)

        # compute composite image
        compositeProduct=computeCompositeSDD(composite_prod)
            
        # write the composite product
        bdate_obj=datetime.date.fromordinal(time_window[tw_ind][0])
        bdate=str(bdate_obj.year)+str(bdate_obj.month)+str(bdate_obj.day)
        edate_obj=datetime.date.fromordinal(time_window[tw_ind][1])
        edate=str(edate_obj.year)+str(edate_obj.month)+str(edate_obj.day)

        wfname='MER_FRS_2PPBCM_'+bdate+'_'+edate+'.dim'
        wfilename=direc+'/composite_sd_product_1/'+wfname
        ProductIO.writeProduct(compositeProduct,wfilename,'BEAM-DIMAP')
"""
