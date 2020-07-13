"""
Created on Mon Jun 16 17:23:00 2020

Creates a composite image
Ref: Stumpf et al. (2012)

@author: Abhinav Gupta
"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import gpfOP

from snappy import ProductIO

def computeCompositeCI(compositeProduct):
    
    # read all the CI bands
    bandnames=list(compositeProduct.getBandNames())
    CI_bandnames=[var for var in bandnames if var.split('_')[0]=='CI' and len(var.split('_'))>1 ]

    # Apply bandmaths in for loop
    newBandName = 'composite_CI_0'
    datatype = 'float32'
    expression = CI_bandnames[0]
    noDataVal='nan'
    CI_prod=gpfOP.BandMathsAG(compositeProduct,newBandName,datatype,expression,noDataVal)

    # merge
    geographicError='NaN'
    compProduct=gpfOP.MergeAG(compositeProduct,CI_prod,geographicError)

    for comp_ind in range(0,len(CI_bandnames)):
        # Apply BandMath
        newBandName = 'composite_CI_'+str(comp_ind+1)
        datatype = 'float32'
        expression = 'if (not nan(composite_CI_'+str(comp_ind)+') and not nan('+CI_bandnames[comp_ind]+')) then max(composite_CI_'+str(comp_ind)+','+CI_bandnames[comp_ind]+')'+' else (if not nan(composite_CI_'+str(comp_ind)+') then composite_CI_'+str(comp_ind)+' else '+CI_bandnames[comp_ind]+')'
        CI_prod=gpfOP.BandMathsAG(compProduct,newBandName,datatype,expression,noDataVal)

        # Merge
        compProduct=gpfOP.MergeAG(compProduct,CI_prod,geographicError)

    return compProduct

