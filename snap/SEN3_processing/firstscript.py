 # -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:37:51 2020

@author: Administrator
"""
import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import geopandas as gpd
import os
import gpfOP

from snappy import ProductIO
from snappy import jpy
from snappy import GPF
from snappy import HashMap
from snappy import Mask

SubsetOp = jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')

"""
direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel/2016'
fname = 'S3A_OL_2_WFR_20161010T155429_160.SEN3'
filename = direc+'/'+fname+'/'+'xfdumanifest.xml'

product=ProductIO.readProduct(filename)


mask_file='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Lake_Erie_mask/subset_0_of_Lake_Erie_mask.dim'
mask_product=ProductIO.readProduct(mask_file)

# read product and display size of the product
#reflec_8=product.getBand('reflec_8');
#Width8=reflec_8.getRasterWidth();
#Height8=reflec_8.getRasterHeight();
#print("reflec_8 size :" +str(Width8)+","+str(Height8))

#reflec_9=product.getBand('reflec_9');
#Width9=reflec_9.getRasterWidth();
#Height9=reflec_9.getRasterHeight();
#print("reflec_9 size :" +str(Width9)+","+str(Height9))

#reflec_8_data = np.zeros(Width8*Height8, dtype=np.float32)
#reflec_8.readPixels(0,0,Width8,Height8,reflec_8_data)
#reflec_8_data.shape = (Height8, Width8)

#reflec_9_data = np.zeros(Width9*Height9, dtype=np.float32)
#reflec_9.readPixels(0,0,Width9,Height9,reflec_9_data)
#reflec_9_data.shape = (Height9, Width9)

#diff_data=reflec_9_data-reflec_8_data

#plt.figure(figsize=(8, 8)) # adjusting the figure window size
#fig = plt.imshow(diff_data, cmap = cm.gray) #matplotlib settings for the current image
#fig.axes.get_xaxis().set_visible(False)
#fig.axes.get_yaxis().set_visible(False)
#plt.show()

# get a spatial subset of the product and mask product according to geometry defined below
wkt = "POLYGON((-83.7117 41.3894, -81.9583 41.1292, -81.7111 42.0236, -83.4917 42.2736, -83.7117 41.3894))"
geometry = WKTReader().read(wkt)

op = SubsetOp()
op.setSourceProduct(product)
op.setGeoRegion(geometry)
product = op.getTargetProduct()

op1 = SubsetOp()
op1.setSourceProduct(mask_product)
op1.setGeoRegion(geometry)
mask_product = op1.getTargetProduct()

Oa10=product.getBand('Oa10_reflectance')
Widthp=Oa10.getRasterWidth();
Heightp=Oa10.getRasterHeight();
print("Oa10 size in subset product: "+str(Widthp)+','+str(Heightp))

LE_Mask=mask_product.getBand('LE_Mask')
Width=LE_Mask.getRasterWidth();
Height=LE_Mask.getRasterHeight();
print("Mask size in subset product: "+str(Width)+','+str(Height))


# obtain wkt of the of subsetted product
lat=product.getRasterDataNode('latitude')
lat_data=np.zeros(Width*Height, dtype=np.float32)
lat.readPixels(0,0,Width,Height,lat_data)
lat_data.shape=(Height,Width)

lat_top_left=lat_data[0,0]
lat_bottom_left=lat_data[Height-1,0]
lat_bottom_right=lat_data[Height-1,Width-1]
lat_top_right=lat_data[0,Width-1]

long=product.getRasterDataNode('longitude')
long_data=np.zeros(Width*Height, dtype=np.float32)
long.readPixels(0,0,Width,Height,long_data)
long_data.shape=(Height,Width)

long_top_left=long_data[0,0]
long_bottom_left=long_data[Height-1,0]
long_bottom_right=long_data[Height-1,Width-1]
long_top_right=long_data[0,Width-1]

wkt_new="POLYGON(("+str(long_bottom_left)+' '+str(lat_bottom_left)+', '+str(long_bottom_right)+' '+str(lat_bottom_right)+', '+str(long_top_right)+' '+str(lat_top_right)+', '+str(long_top_left)+' '+str(lat_top_left)+', '+str(long_bottom_left)+' '+str(lat_bottom_left)+"))"
geometry_new = WKTReader().read(wkt_new)

op1 = SubsetOp()
op1.setSourceProduct(mask_product)
op1.setGeoRegion(geometry_new)
mask_product = op1.getTargetProduct()

LE_mask=mask_product.getBand('LE_Mask')
Width=LE_mask.getRasterWidth();
Height=LE_mask.getRasterHeight();
print("Mask size in subset product: "+str(Width)+','+str(Height))

LE_mask_data = np.zeros(Width*Height, dtype=np.float32)
LE_mask.readPixels(0,0,Width,Height,LE_mask_data)
LE_mask_data.shape = (Height, Width)

plt.figure(figsize=(8, 8)) # adjusting the figure window size
fig = plt.imshow(LE_mask_data, cmap = cm.gray) #matplotlib settings for the current image
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
plt.show()

# Read Lake Erie shapefile
#shp_file = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/Lake_Erie_shape_file/hydro_p_LakeErie/hydro_p_LakeErie.shp'
#shp = gpd.read_file(shp_file)
#wkt = str(shp['geometry'][1])
#geometry=WKTReader().read(wkt)


# get a spatial subset of a product using shapefile
#par=HashMap()
#par.put('copymetadata',True)
#par.put('geoRegion',geometry)
#product = GPF.createProduct('Subset',par,product)
#mask_product = GPF.createProduct('Subset',par,mask_product)

#reflec_8=product.getBand('reflec_8')
#Width=reflec_8.getRasterWidth();
#Height=reflec_8.getRasterHeight();
#print("reflec_8 size in sub_product: "+str(Width)+','+str(Height))

#reflec_8_data = np.zeros(Width*Height, dtype=np.float32)
#reflec_8.readPixels(0,0,Width,Height,reflec_8_data)
#reflec_8_data.shape = (Height, Width)

#plt.figure(figsize=(8, 8)) # adjusting the figure window size
#fig = plt.imshow(reflec_8_data, cmap = cm.gray) #matplotlib settings for the current image
#fig.axes.get_xaxis().set_visible(False)
#fig.axes.get_yaxis().set_visible(False)
#plt.show()

# Resample operator
param=HashMap()
param.put('targetWidth',Widthp)
param.put('targetHeight',Heightp)
param.put('upsampling','Nearest')
mask_product = GPF.createProduct('Resample', param, mask_product)

LE_Mask=mask_product.getBand('LE_Mask')
Width=LE_Mask.getRasterWidth();
Height=LE_Mask.getRasterHeight();
print("Mask size in subset product: "+str(Width)+','+str(Height))

LE_Mask_data = np.zeros(Width*Height, dtype=np.float32)
LE_Mask.readPixels(0,0,Width,Height,LE_Mask_data)
LE_Mask_data.shape = (Height, Width)

# Merge operator
sourceProducts= HashMap()
sourceProducts.put('masterProduct', product)
sourceProducts.put('slaveProduct', mask_product)
parameters = HashMap()
parameters.put('geographicError','NaN')
mg_prd = GPF.createProduct('Merge', parameters, sourceProducts)

# BandMaths operator
BandDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

targetBand = BandDescriptor()
targetBand.name = 'CI'
targetBand.type = 'float32'
targetBand.expression = 'LE_Mask == 1? (Oa08_reflectance-Oa10_reflectance) + (Oa11_reflectance - Oa08_reflectance)*(680-664)/(708-664): NaN'
targetBand.noDataValue=float("nan")
targetBand.validExpression='CI>0'

targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
targetBands[0] = targetBand

params = HashMap()
params.put('targetBands', targetBands)
result = GPF.createProduct('BandMaths', params, mg_prd)

#targetBand = BandDescriptor()
#targetBand.name = 'CI'
#targetBand.type = 'float32'
#targetBand.expression = 'CI>0 ? CI: NaN'
#targetBand.noDataValue=float("nan")

#targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
#targetBands[0] = targetBand

#params = HashMap()
#params.put('targetBands', targetBands)
#result = GPF.createProduct('BandMaths', params, result)

CI=result.getBand('CI')
WidthCI=CI.getRasterWidth()
HeightCI=CI.getRasterHeight()
CI_data = np.zeros(WidthCI*HeightCI, dtype=np.float32)
CI.readPixels(0,0,WidthCI,HeightCI,CI_data)
CI_data.shape = (HeightCI, WidthCI)

plt.figure(figsize=(8, 8)) # adjusting the figure window size
fig = plt.imshow(CI_data, cmap = cm.gray) #matplotlib settings for the current image
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
plt.show()
"""
# check if frcation of clouds is greater than a threshold (not working yet)
# sum_le_mask=np.nansum(LE_Mask_data);

#BandDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

#targetBand = BandDescriptor()
#targetBand.name = 'water_mask'
#targetBand.type = 'float32'
#targetBand.expression = 'water? 1 : 0'
#targetBand.noDataValue=float("nan")

#targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
#targetBands[0] = targetBand

#params = HashMap()
#params.put('targetBands', targetBands)
#result1 = GPF.createProduct('BandMaths', params, product)

#sourceProducts= HashMap()
#sourceProducts.put('masterProduct', result1)
#sourceProducts.put('slaveProduct', mask_product)
#parameters = HashMap()
#parameters.put('geographicError','NaN')
#mg_prd2 = GPF.createProduct('Merge', parameters, sourceProducts)

#targetBand = BandDescriptor()
#targetBand.name = 'water_LE_mask'
#targetBand.type = 'float32'
#targetBand.expression = 'water_mask*LE_Mask'
#targetBand.noDataValue=float("nan")

#targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
#targetBands[0] = targetBand

#params = HashMap()
#params.put('targetBands', targetBands)
#result1 = GPF.createProduct('BandMaths', params, mg_prd2)

#water_LE_mask=result1.getBand('water_LE_mask')
#Width=water_LE_mask.getRasterWidth()
#Height=water_LE_mask.getRasterHeight()
#water_LE_mask_data = np.zeros(Width*Height, dtype=np.float32)
#water_LE_mask.readPixels(0,0,Width,Height,water_LE_mask_data)
#water_LE_mask_data.shape = (Height, Width)

#cloud_free_region=np.nansum(water_LE_mask_data)

#plt.figure(figsize=(8, 8)) # adjusting the figure window size
#fig = plt.imshow(water_LE_mask_data, cmap = cm.gray) #matplotlib settings for the current image
#fig.axes.get_xaxis().set_visible(False)
#fig.axes.get_yaxis().set_visible(False)
#plt.show()

# Size of algae 
#CI_size=np.nansum(CI_data)


# code for mosaicing
"""
direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel/2016'
fname = 'S3A_OL_2_WFR_20160801T160925_81.dim'
filename = direc+'/'+fname
product1 = ProductIO.readProduct(filename)

fname = 'S3A_OL_2_WFR_20160803T151728_83.dim'
filename = direc+'/'+fname
product2 = ProductIO.readProduct(filename)

# mosaicing by using the basic operator

products = jpy.array('org.esa.snap.core.datamodel.Product', 2)
Variable = jpy.get_type('org.esa.snap.core.gpf.common.MosaicOp$Variable')
vars = jpy.array('org.esa.snap.core.gpf.common.MosaicOp$Variable', 1)

products[0] = product1
products[1] = product2

parameters = HashMap()
parameters.put('westBound',-83.60220985614447)
parameters.put('eastBound',-81.78879308692983)
parameters.put('northBound',42.31086488661365)
parameters.put('southBound',41.1893175594222)
parameters.put('pixelSizeY',0.002)
parameters.put('pixelSizeX',0.002)
#parameters.put('resampling','Nearest')
parameters.put('variables', vars)

vars[0] = Variable('CI_demo','CI')       
Mosaic = GPF.createProduct('Mosaic', parameters, products)

newBandName = 'CI_demo'
datatype = 'float32'
expression = 'CI_demo_count ==1?CI_demo:NaN'
noDataVal = 'NaN'
Mosaic = gpfOP.BandMathsAG(Mosaic,newBandName,datatype,expression,noDataVal)

# mosaicing by using gpfOP

westBound = -83.60220985614447
eastBound = -81.78879308692983
northBound = 42.31086488661365
southBound = 41.1893175594222
pixelSizeY = 0.002
pixelSizeX = 0.002
varName = 'CI_demo'
varExpression = 'CI'
Mosaic = gpfOP.mosaicAG(product1,product2,westBound,eastBound,northBound,southBound,pixelSizeX,pixelSizeY,varName,varExpression)


reflec_8=Mosaic.getBand('CI_demo')
Width=reflec_8.getRasterWidth();
Height=reflec_8.getRasterHeight();
print("reflec_8 size in sub_product: "+str(Width)+','+str(Height))

reflec_8_data = np.zeros(Width*Height, dtype=np.float32)
reflec_8.readPixels(0,0,Width,Height,reflec_8_data)
reflec_8_data.shape = (Height, Width)

plt.figure(figsize=(8, 8)) # adjusting the figure window size
fig = plt.imshow(reflec_8_data, cmap = cm.gray) #matplotlib settings for the current image
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
plt.show()
"""
#save_name='im4.png'
#savefile=direc+'/'+save_name
#plt.savefig(savefile,dpi=300,quality=100)
#plt.close()
"""
wfname='mosaic_demo'
write_filename=direc+'/'+wfname+'.dim'
ProductIO.writeProduct(Mosaic,write_filename,'BEAM-DIMAP')
"""

# Merge LE_mask and composite product to compute cloud cover
direc = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel/2016/composite_product'
fname = 'S3A_OL_2_WFR_201651_2016510.dim'
filename = direc+'/'+fname
product1 = ProductIO.readProduct(filename)

fname = 'subset_0_of_Lake_Erie_mask.dim'
filename = 'D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/Lake_Erie_mask/'+fname
product2 = ProductIO.readProduct(filename)

prd = gpfOP.MergeAG(product1,product2,'NaN')
wfname='merge_demo'
write_filename=direc+'/'+wfname+'.dim'
ProductIO.writeProduct(prd,write_filename,'BEAM-DIMAP')
