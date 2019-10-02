""" -----------------------------------------------
Changing Friction
Fassett implementation
======= """


from anuga.operators.base_operator import Operator
from anuga import Region


import math
import numpy as np
from frictionfunctions import *
from model_params import dthresh, gravity, grainD, constantn, frictionscheme 

class friction_operator(Operator, Region)  :
       
    def __init__(self,
                 domain,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
    
   

        Operator.__init__(self, domain, description, label, logging, verbose)

        Region.__init__(self, domain,
                        indices=indices,
                        polygon=polygon,
                        center=center,
                        radius=radius,
                        verbose=verbose)

           
        # set some alaises            
        self.depth = self.domain.quantities['height'].centroid_values


    def __call__(self):
        """
        Applies shear transport to those triangles defined in indices by the specified polygon
        indices == [], then don't apply anywhere
        otherwise apply for the specific indices within the erosion area polygon
        """
        
        updated = True             # flag called OK
        
        if self.indices is not []:     # empty list no polygon - return
   
            #-------------------------------------------
            # model parameters
            #-------------------------------------------   
            #dthresh  
            #grainD 
            #gravity
            nd       = constantn
            
            # Now lets work where water actually is.
            ind = (self.depth >= dthresh)                # indices of triangles in polygon (where depth is above depth threshold: depth>dthresh)
            depth = self.depth[ind]
            
            if frictionscheme=="wilsonetal":
                sqrteightdivfc=frictionfactor_wilsonetal_sq8cf(depth, grainD)
            elif frictionscheme=="larsenandlamb":
                sqrteightdivfc=frictionfactor_larsenandlamb_sq8cf(depth)
            else:
                sqrteightdivfc=frictionfactor_constant_sq8cf(depth)
            
            oneovern = sqrteightdivfc*(gravity**0.5)*(depth**-(1.0/6.0))                   #This is 1/(Manning friction factor) -- (1/n) -- for a given sqrt(8/fc).
            n = 1.0/oneovern                                                         #The inverse of the previous line is the actual n.
            newfriction = np.zeros_like(self.depth)+nd
            newfriction[ind] = n
            self.domain.set_quantity('friction', newfriction, location='centroids', smooth=True)
           
        return (updated)
                
    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        
        return True
