"""
Created on Mon Jun 16 17:23:00 2020

Creates a composite image of secchi depth
cf. Stumpf et al. (2012) composite image construction

@author: Abhinav Gupta
"""

import sys
sys.path.append('C:\\Users\\Administrator\\.snap\\snap-python') # or sys.path.insert(1, '<snappy-dir>')
import snappy
import gpfOP

from snappy import ProductIO

def computeCompositeSDD(compositeProduct):
    
    # read all the CI bands
    bandnames=list(compositeProduct.getBandNames())
    SDD_bandnames=[var for var in bandnames if var.split('_')[0]=='secchi' and len(var.split('_'))>5 ]

    # Apply bandmaths in for loop
    newBandName = 'composite_secchi_depth_over_pos_CI_0'
    datatype = 'float32'
    expression = SDD_bandnames[0]
    noDataVal='nan'
    SDD_prod=gpfOP.BandMathsAG(compositeProduct,newBandName,datatype,expression,noDataVal)

    # merge
    geographicError='NaN'
    compProduct=gpfOP.MergeAG(compositeProduct,SDD_prod,geographicError)

    for comp_ind in range(0,len(SDD_bandnames)):
        # Apply BandMath
        if comp_ind == len(SDD_bandnames)-1:
            newBandName = 'composite_secchi_depth_over_pos_CI_final'
        else:
            newBandName = 'composite_secchi_depth_over_pos_CI_'+str(comp_ind+1)
        datatype = 'float32'
        expression = 'if (not nan(composite_secchi_depth_over_pos_CI_'+str(comp_ind)+') and not nan('+SDD_bandnames[comp_ind]+')) then max(composite_secchi_depth_over_pos_CI_'+str(comp_ind)+','+SDD_bandnames[comp_ind]+')'+' else (if not nan(composite_secchi_depth_over_pos_CI_'+str(comp_ind)+') then composite_secchi_depth_over_pos_CI_'+str(comp_ind)+' else '+SDD_bandnames[comp_ind]+')'
        SDD_prod=gpfOP.BandMathsAG(compProduct,newBandName,datatype,expression,noDataVal)

        # Merge
        compProduct=gpfOP.MergeAG(compProduct,SDD_prod,geographicError)

    return compProduct

