""" -----------------------------------------------
Friction functions
Fassett implementation
======= """
import numpy as np
from model_params import gravity, ks, constantn 


def frictionfactor_larsenandlamb_n():
    # NOTE that this returns n, not sqrt eight div cf 
    n=(gravity**(-0.5))*(1.0/8.1)*(ks**(-1.0/6.0))    #Larsen and Lamb, eq. 3 of methods.  Verified by CF.
    
    return n

def frictionfactor_wilsonetal_sq8cf(depth, grainD):
    #the friction factor here cf is (f_d/8), where fd is a darcy friction factor.  See e.g., Wilson et al. 2004 for empiricisms.
    if grainD > 0.064:
        sqrteightdivfc=5.62*np.log10(depth/grainD)+4   #Chosen from Wilson et al. eq 15.  boulders
    elif grainD > 0.002:
        sqrteightdivfc=5.75*np.log10(depth/grainD)+3.514   #Chosen from Wilson et al. eq 14.  gravel
    elif grainD > 0.0001:
        sqrteightdivfc=8.46*(depth/grainD)**0.1005   #Chosen from Wilson et al. eq 13.  sand
    else:
        raise ValueError('The grain diameter is set too small.')
    
    return sqrteightdivfc

def frictionfactor_larsenandlamb_sq8cf(depth):
    sqrteightdivfc=(gravity**-0.5)*(depth**(1.0/6.0))/frictionfactor_larsenandlamb_n()
        
    return sqrteightdivfc

def frictionfactor_constant_sq8cf(depth):
    sqrteightdivfc=(gravity**-0.5)*(depth**(1.0/6.0))/constantn
        
    return sqrteightdivfc                   