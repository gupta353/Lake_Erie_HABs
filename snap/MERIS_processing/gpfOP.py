"""
Created on Wed Jun 17 13:00:00 2020

GPF operatorsd are implemented in this module

@author: Abhinav Gupta
"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from snappy import ProductIO
from snappy import jpy
from snappy import HashMap
from snappy import GPF

# BandMath operator
def BandMathsAG(product,newBandName,datatype,expression):

    BandDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')    
    targetBand = BandDescriptor()
    targetBand.name = newBandName
    targetBand.type = datatype
    targetBand.expression = expression

    targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
    targetBands[0] = targetBand

    params = HashMap()
    params.put('targetBands', targetBands)
    outputProduct = GPF.createProduct('BandMaths', params, product)
    return outputProduct

# Merge operator
def MergeAG(masterProduct,slaveProduct,geographicError):

    sourceProducts= HashMap()
    sourceProducts.put('masterProduct', masterProduct)
    sourceProducts.put('slaveProduct', slaveProduct)
    parameters = HashMap()
    parameters.put('geographicError',geographicError)
    product = GPF.createProduct('Merge', parameters, sourceProducts)
    return product

# plot image
def plot_image(product,bandname):

    band=product.getBand(bandname)
    Width=band.getRasterWidth()
    Height=band.getRasterHeight()
    band_data=np.zeros(Width*Height, dtype=np.float32)
    band.readPixels(0,0,Width,Height,band_data)
    band_data.shape = (Height, Width)

    plt.figure(figsize=(8, 8)) # adjusting the figure window size
    fig = plt.imshow(band_data, cmap = cm.gray) #matplotlib settings for the current image
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
    plt.show()

# subset operator
def subsetAG(product,wkt):

    SubsetOp = jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
    WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')

    geometry = WKTReader().read(wkt)

    op = SubsetOp()
    op.setSourceProduct(product)
    op.setGeoRegion(geometry)
    subproduct = op.getTargetProduct()
    return subproduct
