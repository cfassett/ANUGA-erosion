""" -----------------------------------------------
Bedload transport 
Fassett implementation
======= """


from anuga.operators.base_operator import Operator
from anuga import Region
from frictionfunction import *
import math
import numpy as np
np.set_printoptions(threshold='nan')


class shear_erosion_operator(Operator, Region)  :
       
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
        self.Wd       = 1000       # water density kg/m3
        self.Sd       = 2700       # sed grain density kg/m3.  
        #self.grainD   = 0.05      # Used to define locally, now defining globally
        #self.G        = self.gravity     # "  "  " " " " " "" " " " "  " " "       
        self.dthresh  = 0.01        # minimum depth to consider (m)
        self.porosity = 0.200      # I would expect realistic values 0 to 0.5. Could split out an expected porosity for the sediment that is eroded and deposited, if more complication is desired.
        self.tau_crit = 0.05       # should probably be something like 0.04-0.08
        self.m        = 1.5        # exponent on sediment transport increase with excess dimensionless shear stress
        self.Ke       = 4          # leading coefficient on qb*.  Roughly 4 is the Wong and Parker 2006 corrected MPM (their eq 24)
        self.maxrate  = 0.02       # This is artificial and is designed to keep the erosion from running away aphysically.  This corresponds to ~1 m / minute, which is still ridiculously fast.     
        
            
        # set some alaises            
        self.stage = self.domain.quantities['stage'].centroid_values        
        self.depth = self.domain.quantities['height'].centroid_values
        self.elev = self.domain.quantities['elevation'].centroid_values
        self.xmom = self.domain.quantities['xmomentum'].centroid_values
        self.ymom = self.domain.quantities['ymomentum'].centroid_values

    def __call__(self):
        """
        Applies shear transport to those triangles defined in indices by the specified polygon
        indices == [], then don't apply anywhere
        otherwise apply for the specific indices within the erosion area polygon
        """
        
        updated = True             # flag called OK
        assert self.domain.get_using_discontinuous_elevation()  # make sure the model is allowing for discontinuous elevations at inner timesteps
        if self.indices is not []:     # empty list no polygon - return


    
            #-------------------------------------------
            # get some useful model parameters
            #-------------------------------------------   
            Wd       = self.Wd
            Sd       = self.Sd
            R        = (Sd/Wd)-1
            G        = self.G
            grainD   = self.grainD
            dthresh  = self.dthresh        
            Ke       = self.Ke      
            tau_crit = self.tau_crit
            m        = self.m 
            porosity = self.porosity
            maxrate  = self.maxrate
                       
            dt = self.get_timestep()
    
            #-----------------------------------------------------------------------------------------
            # Compute erosion depths during the timestep and update centroid elevations accordingly 
			# Note this operator is called for each seperate erosion polygon.
            #-----------------------------------------------------------------------------------------
            
            
			
            ind = (self.depth >= dthresh)                # indices of triangles in polygon (where depth is above depth threshold: depth>dthresh)
            if len(ind)>0:
                height = self.stage_c - self.elev_c  # store for later use
                cellarea = self.areas[ind]
                depth = self.depth[ind]
                xmom = self.xmom[ind]
                ymom  = self.ymom[ind]
            
                             
                vel=(np.sqrt((xmom**2+ymom**2))/(depth+1.0e-8))
                sqrteightdivfc=frictionfactor(depth, grainD)   # sqrt(8/fc) where fc is the Darcy-Weisbach friction factor.  I use Wilson et al. 2004's equations.
                cf=(1.0/sqrteightdivfc)**2.0                   
                tau_b=Wd*cf*(vel**2)
                tau_star=tau_b/(Wd*R*G*grainD)

                excessdimshear=(tau_star-tau_crit)
                excessdimshear[excessdimshear<0]=0  

                qb_star=Ke*(excessdimshear**m)

                qb_dim=qb_star*((R*G)**0.5)*(grainD**1.5)  #qb is the transport volume per unit width per unit time [m^3/s]/m.
               
                qb_dim_x=np.zeros_like(self.depth)
                qb_dim_y=np.zeros_like(self.depth)
                dz=np.zeros_like(self.depth)
                velx=xmom/(depth+1.0e-8)
                vely=ymom/(depth+1.0e-8)
                tbscalingx= Wd*cf*vel*velx
                tbscalingy= Wd*cf*vel*vely
                qb_dim_x[ind]= qb_dim*(tbscalingx/(tau_b+1e-8))   # This method of partitioning the sediment  
                qb_dim_y[ind]= qb_dim*(tbscalingy/(tau_b+1e-8))   # follows Gary Parker, eq. 2.12 (CISMnot2.pdf) 
                self.domain.quantities['sedvolx'].set_values(qb_dim_x,location='centroids')
                self.domain.quantities['sedvoly'].set_values(qb_dim_y,location='centroids')
                self.sedxforgrad=self.domain.quantities['sedvolx']
                self.sedxforgrad.compute_gradients()
                self.sedyforgrad=self.domain.quantities['sedvoly']
                self.sedyforgrad.compute_gradients()
                self.gradX=self.sedxforgrad.x_gradient
                self.gradY=self.sedyforgrad.y_gradient			
                dzin=-(1/(1-porosity))*(self.gradX[ind]+self.gradY[ind])*dt      #exner.  force [ind] to only look in wetted cells.
                dzin[dzin>(maxrate*dt)]=maxrate*dt                       #rate limiters for stability
                dzin[dzin<(-maxrate*dt)]=-maxrate*dt                     #   "    " 
                dz[ind]=dzin                                             #assign dzs to broader domain.
                conservationcheck=np.sum(dz*self.areas)
                if conservationcheck>0:
                    sumerodedvolume=-np.sum(dz[dz<0]*self.areas[dz<0])
                    sumshouldhaveeroded=sumerodedvolume+conservationcheck
                    if sumerodedvolume>0:
                        dz[dz<0]=dz[dz<0]*(sumshouldhaveeroded/sumerodedvolume)
                    else:
                        dz=0
                else:
                    sumdepositedvolume=np.sum(dz[dz>0]*self.areas[dz>0])
                    sumshouldhavedeposited=sumdepositedvolume-conservationcheck            
                    if sumdepositedvolume>0:
                        dz[dz>0]=dz[dz>0]*(sumshouldhavedeposited/sumdepositedvolume)
                    else:
                        dz=0
                conservationcheck=np.sum(dz*self.areas)
                
                newelevs=np.zeros_like(self.elev)
                newelevs=self.elev+dz
                self.domain.set_quantity('elevation', newelevs, location='centroids') 
         
                newstage=np.zeros_like(self.stage)
                newstage=self.stage
                newstage=newelevs + height
                self.domain.set_quantity('stage', newstage, location='centroids') 
                self.domain.set_quantity('height', newstage-newelevs, location='centroids')
                self.domain.distribute_to_vertices_and_edges()         

        self.domain.update_ghosts(['elevation', 'stage', 'sedvolx','sedvoly'])
              

        return (updated)      
        
    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        
        return True
        
                        

