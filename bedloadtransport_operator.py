""" -----------------------------------------------
Bedload transport 
Fassett implementation

9/30/2019
Not convinced that this version handles different cell areas correctly for conservation.  
Also need to check the maxrate implementation around the same thing.
Need to crosscheck this with the Coholich version

======= """


from anuga.operators.base_operator import Operator
from anuga import Region
from frictionfunctions import *
import math
import numpy as np
from model_params import Wd, Sd, porosity, grainD, frictionscheme, gravity, dthresh, maxrate, bed_tau_crit, m, Ke


class bedloadtransport_operator(Operator, Region)  :
       
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
            # set some useful model parameters and aliases
            #-------------------------------------------   

            R        = (Sd/Wd)-1.0
            G        = gravity
            tau_crit = bed_tau_crit
                       
            dt = self.get_timestep()
            maxdz    = maxrate*dt
            #-----------------------------------------------------------------------------------------
            # Compute erosion depths during the timestep and update centroid elevations accordingly 
			# Note this operator is called for each seperate erosion polygon.
            #-----------------------------------------------------------------------------------------
            
            
			
            ind = (self.depth >= dthresh)                # indices of triangles in polygon (where depth is above depth threshold: depth>dthresh)
            if len(ind)>0:
                height = self.stage_c[ind] - self.elev_c[ind]  # store for later use
                cellarea = self.areas[ind]
                depth = self.depth[ind]
                xmom = self.xmom[ind]
                ymom  = self.ymom[ind]
            
                             
                vel=(np.sqrt((xmom**2.0+ymom**2.0))/(depth+1.0e-8)) #velocity in Anuga is magnitude of the momentum vector/depth.  small eta 1.0e-8 is to avoid dividing by zero.
                
                if frictionscheme=="wilsonetal":
                    sqrteightdivfc=frictionfactor_wilsonetal_sq8cf(depth, grainD)      #use frictionfactor function to get friction (cf).
                elif frictionscheme=="larsenandlamb":
                    sqrteightdivfc=frictionfactor_larsenandlamb_sq8cf(depth)
                else:
                    sqrteightdivfc=frictionfactor_constant_sq8cf(depth)
                cf=(1.0/sqrteightdivfc)**2.0                   
                tau_b=Wd*cf*(vel**2.0)                            #Shear Stress
                tau_star=tau_b/(Wd*R*G*grainD)                  #Dimensionless shear stress

                excessdimshear=(tau_star-tau_crit)
                excessdimshear[excessdimshear<0]=0  

                qb_star=Ke*(excessdimshear**m)

                qb_dim=qb_star*((R*G)**0.5)*(grainD**1.5)       #qb is the transport volume per unit width per unit time [m^3/s]/m.
               
                qb_dim_x=np.zeros_like(self.depth)
                qb_dim_y=np.zeros_like(self.depth)

                velx=xmom/(depth+1.0e-8)
                vely=ymom/(depth+1.0e-8)
                tbscalingx= Wd*cf*vel*velx                      #see Gary Parker's CISMnot2 notes; partitioning of shear stress into X and Y components
                tbscalingy= Wd*cf*vel*vely
                qb_dim_x[ind]= qb_dim*(tbscalingx/(tau_b+1e-8))   # This method of partitioning the sediment flux
                qb_dim_y[ind]= qb_dim*(tbscalingy/(tau_b+1e-8))   # follows Gary Parker, eq. 2.12 (CISMnot2.pdf) 
                self.domain.quantities['sedvolx'].set_values(qb_dim_x,location='centroids')  #These are quantities mostly for convenience (because it is useful to use the 
                self.domain.quantities['sedvoly'].set_values(qb_dim_y,location='centroids')  #ANUGA built in compute_gradients given the uneven mesh).
                self.sedxforgrad=self.domain.quantities['sedvolx']
                self.sedxforgrad.compute_gradients()
                self.sedyforgrad=self.domain.quantities['sedvoly']
                self.sedyforgrad.compute_gradients()
                self.gradX=self.sedxforgrad.x_gradient
                self.gradY=self.sedyforgrad.y_gradient			
   
                dz=-(1/(1-porosity))*(self.gradX+self.gradY)*dt      #exner.
                dz[self.depth<dthresh]=0.0
                numoferoding=np.count_nonzero(dz<0)
                numofdeposit=np.count_nonzero(dz>0)
                if (numoferoding>0) and (numofdeposit>0):            #stop if there isn't both erosion and deposition.
                    averageachange = self.areas[abs(dz)>0].mean()
                    dz[dz>0]=dz[dz>0]*averageachange/self.areas[dz>0]
                    dz[dz<0]=dz[dz<0]*averageachange/self.areas[dz<0]
                    conserror=np.sum(dz[dz>0]*self.areas[dz>0])+np.sum(dz[dz<0]*self.areas[dz<0])  #volume extra dep (+) or eroded (-)




    # was trying to split conservation fix between erosion and deposition 50/50 gave up.
    #               dz[dz<0]=dz[dz<0]-0.5*conserror/(self.areas[dz<0]*numoferoding)  OLD
                    dz[dz>0]=dz[dz>0]-conserror/(self.areas[dz>0]*numofdeposit)


                    
                    if np.max(np.abs(dz))>maxdz:                       #rate limiter.  Maxrate is the most that a cell can change in a timestep.  This fights a potential stability problems.
                        dz *= maxdz/np.max(np.abs(dz))                         
                    self.elev_c[ind] = self.elev_c[ind] + dz[ind]    
                    self.domain.set_quantity('elevation', self.elev_c, location='centroids')                            
                    # Make sure stage is corrected once erosion happens to conserve fluid volume
                    self.stage_c[ind] = self.elev_c[ind] + height                                       
                    self.domain.set_quantity('stage', self.stage_c, location='centroids')
                    self.domain.distribute_to_vertices_and_edges()
 
              

        return (updated)      
              
        
    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        
        return True
        
                        

