""" -----------------------------------------------
Changing Friction
Fassett implementation
======= """


from anuga.operators.base_operator import Operator
from anuga import Region


import math
import numpy as np
from frictionfunction import *
np.set_printoptions(threshold='nan')


class friction_operator(Operator, Region)  :
       
    def __init__(self,
                 domain,
                 threshold=0.0,                 
                 grainD=0.001,
                 gravity=3.711,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
    
   

        Operator.__init__(self, domain, description, label, logging, verbose)
        self.grainD=grainD
        self.G=gravity
 

        Region.__init__(self, domain,
                        indices=indices,
                        polygon=polygon,
                        center=center,
                        radius=radius,
                        verbose=verbose)

        #-------------------------------------------
        # set properties
        #-------------------------------------------
        self.dthresh  = 0.001      # minimum depth to consider (m)
        #self.grainD   = 0.05      # Used to define locally, now defining globally
        #self.G        = self.gravity     # "  "  " " " " " "" " " " "  " " "       
        self.ndefault = 0.0545     # this is only used for zero depth cells, so probably irrelevant.
            
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
            # get some useful model parameters
            #-------------------------------------------   
            dthresh  = self.dthresh 
            grainD   = self.grainD 
            G        = self.G
            nd       = self.ndefault
            
            # Now lets work where water actually is.
            ind = (self.depth >= dthresh)                # indices of triangles in polygon (where depth is above depth threshold: depth>dthresh)
            depth = self.depth[ind]
            
            sqrteightdivfc=frictionfactor(depth, grainD)
            oneovern = sqrteightdivfc*(G**0.5)*(depth**-(1.0/6.0))
            n = 1.0/oneovern
            newfriction = np.zeros_like(self.depth)+nd
            newfriction[ind] = n
            self.domain.set_quantity('friction', newfriction, location='centroids')
           
        return (updated)
                
    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        
        return True
        
                        

