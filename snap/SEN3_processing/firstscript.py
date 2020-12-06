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

from snappy import ProductIO
from snappy import jpy
from snappy import GPF
from snappy import HashMap
from snappy import Mask

SubsetOp = jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')


file_path='D:/Research/EPA_Project/Lake_Erie_HAB/Data/remote_sensing_data/MERIS/MER_FRS_2PPBCM20090518_153522_000000312079_00097_37726_0015.N1'
product=ProductIO.readProduct(file_path)

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

# get a spatial subset of the product
wkt = "POLYGON((-83.56861 41.5027, -82.08055 41.25472, -81.81277 41.9925, -83.32 42.27277, -83.56861 41.5027))"
geometry = WKTReader().read(wkt)

op = SubsetOp()
op.setSourceProduct(product)
op.setGeoRegion(geometry)
product = op.getTargetProduct()

op1 = SubsetOp()
op1.setSourceProduct(mask_product)
op1.setGeoRegion(geometry)
mask_product = op1.getTargetProduct()

reflec_8=product.getBand('reflec_8')
Widthp=reflec_8.getRasterWidth();
Heightp=reflec_8.getRasterHeight();
print("reflec_8 size in subset product: "+str(Widthp)+','+str(Heightp))

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
#param=HashMap()
#param.put('targetWidth',Widthp)
#param.put('targetHeight',Heightp)
#param.put('upsampling','Nearest')
#mask_product = GPF.createProduct('Resample', param, mask_product)

#LE_Mask=mask_product.getBand('LE_Mask')
#Width=LE_Mask.getRasterWidth();
#Height=LE_Mask.getRasterHeight();
#print("Mask size in subset product: "+str(Width)+','+str(Height))

#LE_Mask_data = np.zeros(Width*Height, dtype=np.float32)
#LE_Mask.readPixels(0,0,Width,Height,LE_Mask_data)
#LE_Mask_data.shape = (Height, Width)

# Merge operator
#sourceProducts= HashMap()
#sourceProducts.put('masterProduct', product)
#sourceProducts.put('slaveProduct', mask_product)
#parameters = HashMap()
#parameters.put('geographicError','NaN')
#mg_prd = GPF.createProduct('Merge', parameters, sourceProducts)

# BandMaths operator
#BandDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

#targetBand = BandDescriptor()
#targetBand.name = 'CI'
#targetBand.type = 'float32'
#targetBand.expression = 'LE_Mask == 1 and water? (reflec_7-reflec_8) + (reflec_9 - reflec_7)*(680-664)/(708-664): NaN'
#targetBand.noDataValue=float("nan")
#targetBand.validExpression='CI>0'

#targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
#targetBands[0] = targetBand

#params = HashMap()
#params.put('targetBands', targetBands)
#result = GPF.createProduct('BandMaths', params, mg_prd)

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

#CI=result.getBand('CI')
#WidthCI=CI.getRasterWidth()
#HeightCI=CI.getRasterHeight()
#CI_data = np.zeros(WidthCI*HeightCI, dtype=np.float32)
#CI.readPixels(0,0,WidthCI,HeightCI,CI_data)
#CI_data.shape = (HeightCI, WidthCI)

#plt.figure(figsize=(8, 8)) # adjusting the figure window size
#fig = plt.imshow(CI_data, cmap = cm.gray) #matplotlib settings for the current image
#fig.axes.get_xaxis().set_visible(False)
#fig.axes.get_yaxis().set_visible(False)
#plt.show()

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
