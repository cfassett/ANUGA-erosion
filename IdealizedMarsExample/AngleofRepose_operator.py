""" -----------------------------------------------
Angle of Repose operator
Fassett implementation

9/30/2019
Note this version is a little loose on how it does sediment conservation, particularly between
cells of different area.  If a little bit is lost in the domain, it can be added back in distributed 
among cells, not locally. 

The Coholich AoR version attempted to correct this, but it is much slower and NOT more stable, some 
for the time being I have reverted to this.

======= """


from anuga.operators.base_operator import Operator
from anuga import Region

import math
import numpy as np
import numpy.ma as ma
from model_params import angleofrepose


class AoR_operator(Operator, Region):
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


    def __call__(self):
        """
        Applies Angle of Repose function to domain
        """
              
        if self.indices is not []:     # empty list no polygon - return
            self.AofR()            

            
        
    def AofR(self):
        """
        Applies AoR collapse, everywhere.  Currently doesn't honor indices
        indices == [], then don't apply anywhere
        otherwise apply for the specific indices within the  polygon
        """
        

        self.relaxationfactor = 0.5    #Max value=0.5; Range=0-0.5.  This controls how quickly to handle sediment collapsing downhill to relax oversteepened slopes.  
                                       
                                       #Could add back in an iterator that makes this operator repeat itself, too, since a single iteration does NOT fully relax all oversteepened slopes. 
                                       #Could set in a way that is cognizant of the timestep, perhaps coupled with the number of iterations.   dt = self.get_timestep(), relaxationfactor=somenumber, number of iterations=some number
        self.arperc=math.tan(math.radians(angleofrepose))  # in %
                    
        # set some alaises            
        self.stage = self.domain.quantities['stage'].centroid_values   
        self.depth = self.domain.quantities['height'].centroid_values
        self.elev = self.domain.quantities['elevation'].centroid_values
        
        assert self.domain.get_using_discontinuous_elevation()  # make sure the model is allowing for discontinuous elevations at inner timesteps
        if self.indices is not []:     # empty list no polygon - return
            #-------------------------------------------
            # get some useful model parameters
            #-------------------------------------------   
            arperc = self.arperc
            relaxationfactor = self.relaxationfactor

            neighbours = self.domain.surrogate_neighbours
            height = self.stage - self.elev                            # save depths to conserve volumes
                           
                                              
            #n0=neighbours[:,0]       #un-needed.  Note that this is how neighbours works, though. (cell,neighbor[1-3]))
            #n1=neighbours[:,1]
            #n2=neighbours[:,2]
            k=len(self.elev)
            self.ident  = np.arange(k)
            #e = np.zeros((k,3))
            #e[:,0] = self.elev[n0]
            #e[:,1] = self.elev[n1]
            #e[:,2] = self.elev[n2]   #The thing below is the numpy-ish way of writing these three lines
            e=self.elev[neighbours]
            
            lxy=np.zeros((k,3))
            #lxy[:,0]=np.sqrt((self.coord_c[:,0]-self.coord_c[n0,0])**2 + (self.coord_c[:,1]-self.coord_c[n0,1])**2)
            #lxy[:,1]=np.sqrt((self.coord_c[:,0]-self.coord_c[n1,0])**2 + (self.coord_c[:,1]-self.coord_c[n1,1])**2 )
            #lxy[:,2]=np.sqrt((self.coord_c[:,0]-self.coord_c[n2,0])**2 + (self.coord_c[:,1]-self.coord_c[n2,1])**2 )  #The thing below is the numpy-ish way of writing these three lines
            lxy=np.sqrt((self.coord_c[:,0][:,np.newaxis]-self.coord_c[neighbours,0])**2+(self.coord_c[:,1][:,np.newaxis]-self.coord_c[neighbours,1])**2)
            with np.errstate(divide='ignore', invalid='ignore'):
                slopeperc=np.where(lxy>0.0,(self.elev[:,np.newaxis]-e)/lxy, 0.0)

            steepestdescentindex = np.argmax(slopeperc, axis=1)        
            
            dslope = slopeperc[self.ident,steepestdescentindex]       # what is the slope to the neighboring cell of steepest descent  (positive slope downhill)
         
            oversteepenedcells = dslope>arperc
            dz=np.zeros_like(self.elev)
            for nonzeroelement in np.argwhere(oversteepenedcells):
                heightchange=-relaxationfactor*lxy[nonzeroelement,steepestdescentindex[nonzeroelement]]*arperc
                volumechange=heightchange*self.areas[nonzeroelement]              # note, negative
                dz[nonzeroelement]=dz[nonzeroelement]+heightchange   #drop cells that are much higher than lowest neighbors down to their AOR*rel
                depositcell=neighbours[int(self.ident[nonzeroelement]),int(steepestdescentindex[nonzeroelement])]
                depositingheightchange=(-volumechange)/(self.areas[depositcell])  #note, positive [volumechange is -]
                dz[depositcell]=dz[depositcell]+depositingheightchange
            
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
            
            newelevs=self.elev+dz
            self.domain.set_quantity('elevation', newelevs, location='centroids') 
     
            newstage=self.stage
            newstage=newelevs + height
            self.domain.set_quantity('stage', newstage, location='centroids') 
            self.domain.set_quantity('height', newstage-newelevs, location='centroids')
            self.domain.distribute_to_vertices_and_edges() 
                                     

    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        
        return False

        
                        

