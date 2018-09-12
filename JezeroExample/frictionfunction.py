""" -----------------------------------------------
Friction function
Fassett implementation
======= """
import numpy as np

def frictionfactor(depth, grainD):
    #the friction factor here cf is (f_d/8), where fd is a darcy friction factor.  See e.g., Wilson et al. 2004 for empiricisms.
    if grainD > 0.064:
        sqrteightdivfc=5.62*np.log10(depth/grainD)+4   #Chosen from Wilson et al. eq 15.  boulders
    elif grainD > 0.002:
        sqrteightdivfc=5.75*np.log10(depth/grainD)+3.514   #Chosen from Wilson et al. eq 14.  gravel
    elif grainD > 0.0001:
        sqrteightdivfc=8.46*np.log10((depth/grainD)**0.1005)   #Chosen from Wilson et al. eq 13.  sand
    else:
        raise ValueError('The grain diameter is set too small.')
    
    return sqrteightdivfc
                