""" -----------------------------------------------
Sediment transport in suspension
Fassett implementation
======= """


from anuga.operators.base_operator import Operator
from anuga import Region
from frictionfunctions import *

import math
import numpy as np
from model_params import Wd, Sd, porosity, gravity, grainD, susp_tau_crit, gamma0, maxconc, mu, maxrate, dthresh, frictionscheme

class suspendedtransport_operator(Operator, Region):
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
        self.conc = self.domain.quantities['concentration'].centroid_values         
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
             
            # Move suspended sediment with fluid 
            self.sediment_flux()
    
            #-------------------------------------------
            # set some useful model parameters and names
            #-------------------------------------------   
            R        = (Sd/Wd)-1.0       #relative density
            G        = gravity         #alias           
            tau_crit = susp_tau_crit   #alias            
            settlingvelocity = self.settlingvelocityfunc(Wd,Sd,G,mu,grainD)
                      
            dt = self.get_timestep()
    
            #-----------------------------------------------------------------------------------------
            # Compute erosion depths during the timestep and update centroid elevations accordingly 
			# Note this operator is called for each seperate erosion polygon.
            #-----------------------------------------------------------------------------------------
            


			
            ind = (self.depth >= dthresh)                  # indices of triangles in polygon (where depth is above depth threshold: depth>dthresh)
            height = self.stage_c[ind] - self.elev_c[ind]  # store for later use
            cellarea = self.areas[ind]
            depth = self.depth[ind]
            xmom = self.xmom[ind]
            ymom  = self.ymom[ind]
            conc = self.conc[ind]                       
            vel=(np.sqrt((xmom**2.0+ymom**2.0))/(depth+1.0e-8))   #velocity in Anuga is magnitude of the momentum vector/depth.  small eta 1.0e-8 is to avoid dividing by zero.
            
            if frictionscheme=="wilsonetal":
                sqrteightdivfc=frictionfactor_wilsonetal_sq8cf(depth, grainD)      #use frictionfactor function to get friction (cf).
            elif frictionscheme=="larsenandlamb":
                sqrteightdivfc=frictionfactor_larsenandlamb_sq8cf(depth)
            else:
                sqrteightdivfc=frictionfactor_constant_sq8cf(depth)
                
            cf=(1.0/sqrteightdivfc)**2.0                      
            tau_b=Wd*cf*(vel**2.0)                              #Shear Stress
            tau_star=tau_b/(Wd*R*G*grainD)                    #Dimensionless shear stress

            entrainmentrate_star=np.zeros_like(tau_star)
            entrainmentrate_star[tau_star>tau_crit]=0.65*gamma0*(((tau_star[tau_star>tau_crit]/tau_crit)-1.0)/(1.0+gamma0*((tau_star[tau_star>tau_crit]/tau_crit)-1.0)))   #Smith and McLean Entraiment (dimensionless)
            entrainmentrate=settlingvelocity*entrainmentrate_star
            
            shearvelocity=np.sqrt(tau_b/Wd)
            rousenumber=settlingvelocity/(0.41*shearvelocity+1.0e-8)  #0.41 is von Karman's constant.  small eta to avoid divide by zero
            
            dstar=np.zeros_like(rousenumber)+1.0                      #dstar is the ratio of sediment in the bottom layer, versus the average, calculated in a Rouse-ian way.  (using Rouse concentrations and polynomial appx.)
            dstar[rousenumber<0.8] =3.057*(rousenumber[rousenumber<0.8]**2)+3.4803*rousenumber[rousenumber<0.8]+1
            dstar[(rousenumber>=0.8) & (rousenumber<2.0)]=-1.4802*rousenumber[(rousenumber>=0.8) & (rousenumber<2.0)]**2+11.181*rousenumber[(rousenumber>=0.8) & (rousenumber<2.0)]-2.3601
            dstar[(rousenumber>=2.0) & (rousenumber<4.0)]=-0.9557*rousenumber[(rousenumber>=2.0) & (rousenumber<4.0)]**2+8.0402*rousenumber[(rousenumber>=2.0) & (rousenumber<4.0)]+1.8783
            dstar[(rousenumber>=4.0) & (rousenumber<6.0)]=-0.1753*rousenumber[(rousenumber>=4.0) & (rousenumber<6.0)]**2+2.1908*rousenumber[(rousenumber>=4.0) & (rousenumber<6.0)]+12.908
            dstar[rousenumber>=6.0]=19.76
            
            bottomconcentration=dstar*conc                                    
            depositionrate=settlingvelocity*bottomconcentration
            
            #artificial rate limiter:
            depositionrate[(entrainmentrate-depositionrate)>maxrate]=0.0
            entrainmentrate[(entrainmentrate-depositionrate)>maxrate]=maxrate
            entrainmentrate[(depositionrate-entrainmentrate)>maxrate]=0.0
            depositionrate[(depositionrate-entrainmentrate)>maxrate]=maxrate
            
           
            # limiting condition where everything will deposit
            boom = (depositionrate-entrainmentrate)>(conc*depth)/dt         
            depositionrate[boom]=(conc[boom]*depth[boom]/dt)  # Do not deposit more sediment than is in the column
            entrainmentrate[boom]=0.0
            
            sedvolincell=conc*depth*cellarea                                  #m^3
            changeinsedvol=(entrainmentrate-depositionrate)*cellarea*dt       #m^3 
            newsedvolincell=sedvolincell+changeinsedvol                       #m^3
            newconcentrations=newsedvolincell/(depth*cellarea)                #%
           
            #handle trying to entrain more than maxconc by not doing so...
            extra=np.zeros_like(newconcentrations)
            extra[newconcentrations>maxconc]=newconcentrations[newconcentrations>maxconc]-maxconc  
            newconcentrations[newconcentrations>maxconc]=maxconc
            extraz=extra*depth
            
            dz=-(1.0/(1.0-porosity))*((entrainmentrate-depositionrate)*dt-extraz)
                       
            concentrations=self.conc
            concentrations[ind]=newconcentrations   
            self.domain.set_quantity('concentration', concentrations, location='centroids') 
            
            newelevs=self.elev
            newelevs[ind]=newelevs[ind]+dz
            self.domain.set_quantity('elevation', newelevs, location='centroids') 
                  
            newstage=self.stage
            newstage[ind]=newelevs[ind] + height
            self.domain.set_quantity('stage', newstage, location='centroids') 
            self.domain.distribute_to_vertices_and_edges()
                                  

        return (updated)   
    
    def settlingvelocityfunc(self,Wd,Sd,G,mu,grainD):
        # From Dietrich, 1982: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/WR018i006p01615
        
        #oldsettlingvelocity = (0.055)*(grainD**2.0)*(Sd-Wd)*G*((mu)**-1)   Settling velocity from Stokes.  Note that the more complicated equations below (from Dietrich 1982)
        #                                                                   converge to this old settling velocity for small particles.
        # CAUTION: Dietrich says this equation is invalid for very large dstar (non-dimensional grain size).  For Mars spherical grains at 2700 kg/m3, this is hit for >0.085m [8.5cm]
        #          The reasons are discussed p1617-1618. 
        
        dstar=(Sd-Wd)*G*(grainD**3)/(Wd*((mu/Wd)**2))            #Dietrich equation 6.  Note the kinematic velocity here (mu/Wd) rather than dynamic viscosity.
        
        dstarlogten=math.log10(dstar)
        if dstar>0.05:           
            logwstar=-3.76715+1.92944*dstarlogten-0.09815*(dstarlogten**2)-0.00575*(dstarlogten**3)+0.00056*(dstarlogten**4)    #Dietrich eq9
            wstar=10**(logwstar)
        else:
            wstar=(dstar**2)/5832    #Dietrich equation 8, Stokes result (old settling velocity)
        
        w=(((Sd-Wd)*G*(mu/Wd)*wstar)/Wd)**(1.0/3.0)     #Dietrich equation 5                    #Note the kinematic viscosity here.
            
        return w   
         
        
            
            
            
    def sediment_flux(self):
        """
        Modified from AnugaSED -- Mariela Perignon.
          --altered to maintain conservation of sediment
          --Does not handle inflow at dirichlet boundary cells (unlike original code).
        
        Calculates the flux of sediment between cells based on the flux of water
        to calculate the new concentration at centroids
        
        Assumes that sediment moves at the same speed as the flow.  Imperfect assumption 
        depending on how sediment is distributed in water column.
        
        Negative edge fluxes are inwards and must use the concentration of the neighbour
        Positive edge fluxes are outwards and use the concentration of the cell
        """
        depth = self.domain.quantities['height'].centroid_values
        ind = (depth >= dthresh )        
        depth_e = self.domain.quantities['height'].edge_values  
        xmom_e = self.domain.quantities['xmomentum'].edge_values
        ymom_e = self.domain.quantities['ymomentum'].edge_values
        cellarea = self.areas
        normals = self.domain.normals
        neighbours = self.domain.neighbours
        edgelengths = self.domain.edgelengths
        conc = self.domain.quantities['concentration'].centroid_values
        num_cells = len(depth)
        dt=self.get_timestep()
        
        totalsedvol=np.sum(conc*depth*cellarea)   #define this for conservation / whether to bother doing the rest
        if totalsedvol>0:
            normal_vels = np.zeros((num_cells,3))
            normal_vels[ind,:] = ((normals[ind,0::2] *
                                        xmom_e[ind,:] +
                                        normals[ind, 1::2] *
                                        ymom_e[ind,:]) /
                                        depth_e[ind,:])
                                        

            
            edge_flux = (depth_e * edgelengths * normal_vels * dt)
           
            neighbour_conc = conc[neighbours]
            
            sed_flux = edge_flux * conc[:,np.newaxis]
            
            sed_flux[edge_flux < 0] = (edge_flux[edge_flux < 0] * neighbour_conc[edge_flux < 0])          

            sed_vol_change = np.sum(-sed_flux, axis=1)
            
            sed_vol_in_cell = conc * depth * cellarea
            new_sed_vol_in_cell = np.maximum(sed_vol_in_cell + sed_vol_change, 0)
            
            new_conc=np.zeros_like(depth)
            new_conc[ind] = (new_sed_vol_in_cell[ind]/(self.depth[ind] * self.areas[ind]))
            aftersedvol=np.sum(new_conc*depth*cellarea)
            conservationfactor=aftersedvol/totalsedvol            
            new_conc = new_conc / conservationfactor
            self.domain.set_quantity('concentration', new_conc, location='centroids')
            
        

                      
    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        
        return False

        
                        

