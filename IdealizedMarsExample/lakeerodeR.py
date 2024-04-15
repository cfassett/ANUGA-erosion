#!/usr/bin/env python
from __future__ import print_function

import os, optparse, subprocess, sys, tempfile

# import osgeo.ogr as ogr
# import osgeo.osr as osr
# import osgeo.gdal as gdal

import copy
import anuga
import numpy as np
import math
from suspendedtransport_operator import suspendedtransport_operator
from bedloadtransport_operator import bedloadtransport_operator
from friction_operator import friction_operator
from AngleofRepose_operator import AoR_operator
from analysisfunctions import finalstats
from model_params import gravity, grainD, constantn


def man(option, opt, value, parser):
    print (parser.usage, file=sys.stderr)
    print ('''\
This program does lake erosion in ANUGA.  
''', file=sys.stderr)
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
    try:
        try:
            usage = "usage: lakeerode.py LakeRadius-in-meters [-m depthmultiplier] [-i initialbreachdepth]\n"
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("-m", dest="dm", help="Depth Multiplier for lake")
            parser.add_option("-i", dest="initbreachdepth", help="Initial breach depth for flood")
            (options, inargs) = parser.parse_args()
            if not inargs: parser.error("need parameters")
            lakeradius=float(inargs[0])
          
            dm=1.0                      #default
            if( options.dm ): dm=float(options.dm)
            h=dm*21.39*math.pow((lakeradius/1000),0.62) #max depth in meters to get V=0.01A**1.31 from Fassett and Head 2008

            initbreachdepth=15.0          #default
            lakebelowconfining=5
            if( options.initbreachdepth): initbreachdepth=float(options.initbreachdepth)
            initbreachdepth=initbreachdepth+lakebelowconfining                            #set this to be a real initial head, not dependent on lakebelowconfining 
            extslope=-0.015
            extslopeint=0.015
 
            #anuga critical parameter
            anuga.g=gravity

 
            #domain parameters
            domainpolygon=[[0,0],[0*lakeradius,2*lakeradius],[2*lakeradius,2*lakeradius],[2*lakeradius,1.2*lakeradius],[10*lakeradius,1.2*lakeradius],[10*lakeradius,0.8*lakeradius],[2*lakeradius,0.8*lakeradius],[2*lakeradius,0],[0,0]]
            outletregion=[[1.9*lakeradius,0.8*lakeradius],[1.9*lakeradius,1.2*lakeradius],[4.5*lakeradius,1.2*lakeradius],[4.5*lakeradius,0.8*lakeradius]]

            boundary_tags={'exterior': [0]}

            high_resolution=(lakeradius/100)**2
            base_resolution=high_resolution*80.0
            initbreachwidth=lakeradius/50.0
            interior_regions = [[outletregion, high_resolution]]


            meshname='lake.msh'
            m = anuga.create_mesh_from_regions(domainpolygon, boundary_tags, maximum_triangle_area=base_resolution, interior_regions=interior_regions, filename=meshname, use_cache=False)
            evolved_quantities =  ['stage', 'xmomentum', 'ymomentum', 'concentration','sedvolx','sedvoly']
            #evolved_quantities =  ['stage', 'xmomentum', 'ymomentum', 'concentration']
            domain = anuga.Domain(meshname, use_cache=False, evolved_quantities = evolved_quantities)
            domain.g = anuga.g  # make sure the domain inherits the package's gravity.

            print ('Number of triangles = ', len(domain))

            def find_nearest(array, value):
                array = np.asarray(array)
                idx = (np.abs(array - value)).argmin()
                return idx
                
            def topography(x, y):
                xc=x-lakeradius 
                yc=y-lakeradius
                
                #sets up lake
                ellipsoid=-(h**2*(1-((xc**2+yc**2)/lakeradius**2)))**0.5
                z=np.where(ellipsoid<0,ellipsoid,0.0)
                
                #sets up exterior slope
                z=np.where(np.logical_and(x>lakeradius*2,x<=lakeradius*4),(x-lakeradius*2)*extslope+abs(y-lakeradius)*extslopeint,z)

                #setup flatlying area
                z=np.where(x>lakeradius*4,(lakeradius*2)*extslope+abs(y-lakeradius)*extslopeint,z)
                

                # add gaussian noise, with zero mean and 0.1 m standard deviation.
                #z=z+np.random.normal(0,0.1,z.shape)
                

                
                # set up breach as a failure at the midpoint in a region
                #mask1 = (x>=(lakeradius*2+(initbreachdepth/extslope)))
                mask1 = (x>=lakeradius)
                #mask2 = (x<=(lakeradius*2-(initbreachdepth/extslope)))
                mask3 = (y>=((lakeradius-(initbreachwidth/2.0))))  
                mask4 = (y<=((lakeradius+(initbreachwidth/2.0))))
                mask5 = (z>-initbreachdepth)
                #mask=mask1 & mask2 & mask3 & mask4 & mask5
                mask=mask1 & mask3 & mask4 & mask5
                z[mask]=-initbreachdepth


                return z

            def initialstage(x,y):
                clow=domain.get_quantity('elevation').get_values(interpolation_points=[[2*lakeradius,lakeradius]])
                istage=clow+initbreachdepth-lakebelowconfining
                
                #outsidelake. Force depth=0 (istage=topo) across rim.
                #istage=np.where(x>lakeradius*2,topography(x,y),istage)
                istage=np.where(x>lakeradius*2,-10000,istage)

                return istage

            name="Marslake"+str(lakeradius)+"_"+str(grainD*1000)+"mm"
            domain.set_name(name)
            # Choose between DE0 (less accurate), DE1 (more accurate, slower), DE2 (even more accurate, slower still)
            domain.set_flow_algorithm('DE0')
            domain.set_CFL(cfl=1)
            domain.set_minimum_allowed_height(0.01)  #raise the default minimum depth from 1 mm to 1 cm to speed the simulation (for computation)
            domain.set_minimum_storable_height(0.01) #raise the default minimum depth from 1 mm to 1 cm to speed the simulation (for storage/visualization)
            domain.set_maximum_allowed_speed(1)      #maximum particle speed that is allowed in water shallower than minimum_allowed_height (back of the envelope suggests 1m/s should be below tau_crit)

 

            np.seterr(invalid='ignore')
            domain.set_quantity('elevation', topography, location='centroids') # elevation is a function
            voltotal=np.sum(domain.quantities['elevation'].centroid_values*domain.areas)
            domain.set_quantity('concentration', 0.00)        # Start with no esd
            domain.set_quantity('friction', 0.0545)      # Constant Manning friction.  Only used to initialize domain (calculated every time step).
            domain.set_quantity('stage', initialstage, location='centroids')
            domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 2})

            np.seterr(invalid='warn')

            # Setup boundary conditions

            Br=anuga.Reflective_boundary(domain)
            domain.set_boundary({'exterior': Br})

            operatezero = suspendedtransport_operator(domain)
            operateone = bedloadtransport_operator(domain)
            operatetwo = friction_operator(domain)
            operatethree = AoR_operator(domain)

            initial=domain.get_quantity('elevation').get_values(interpolation_points=[[2*lakeradius,lakeradius]], location='centroids')
            initialdomain=copy.deepcopy(domain)
            ystep=300
            ftime=86400*10
            count=0
            volthresh=1000
            for t in domain.evolve(yieldstep=ystep, finaltime=ftime):
                print (domain.timestepping_statistics())
                print ('xmom:'+str(domain.get_quantity('xmomentum').get_values(interpolation_points=[[2*lakeradius,lakeradius]], location='centroids')))   
                volcurr=np.sum(domain.quantities['elevation'].centroid_values*domain.areas)
                breacherosion=domain.get_quantity('elevation').get_values(interpolation_points=[[2*lakeradius,lakeradius]], location='centroids')-initial
                print ('erosion: '+str(breacherosion))
                volsed=np.sum(domain.quantities['concentration'].centroid_values*domain.quantities['height'].centroid_values*domain.areas)            
                conservation=(volcurr+volsed-voltotal)/voltotal
                print ('conservation: '+'{:.8%}'.format(conservation))
                if (volsed<volthresh)&(count>0):
                    print ("No sediment moving...ending...")
                    break
                count=count+1
            polyline=[[2*lakeradius,lakeradius-2*initbreachwidth],[2*lakeradius,lakeradius+2*initbreachwidth]]
            time, Q = anuga.get_flow_through_cross_section(name+'.sww',polyline)
            print (Q)

            initname="initial_"+str(lakeradius)+"_"+str(grainD*1000)+"mm"+".asc"
            finname="final_"+str(lakeradius)+"_"+str(grainD*1000)+"mm"+".asc"
            fluxname="flux_"+str(lakeradius)+"_"+str(grainD*1000)+"mm"+".txt"
            np.savetxt(fluxname,Q)
            finalstats(domain,initialdomain,lakeradius)
            np.save('XconcC.npy', domain.quantities['concentration'].centroid_values)
            np.save('XelevC.npy', domain.quantities['elevation'].centroid_values)
            np.save('XxmC.npy', domain.quantities['xmomentum'].centroid_values)
            np.save('XymC.npy', domain.quantities['ymomentum'].centroid_values)
            np.save('XstageC.npy', domain.quantities['stage'].centroid_values)
            np.save('XconcV.npy', domain.quantities['concentration'].vertex_values)
            np.save('XelevV.npy', domain.quantities['elevation'].vertex_values)
            np.save('XxmV.npy', domain.quantities['xmomentum'].vertex_values)
            np.save('XymV.npy', domain.quantities['ymomentum'].vertex_values)
            np.save('XstageV.npy', domain.quantities['stage'].vertex_values)

        except optparse.OptionError as err:
            raise Usage(err.msg)

    except Usage as err:
        print (err.msg, file=sys.stderr)
        # print ("for help use --help", file=sys.stderr)
        return 2

if __name__ == "__main__":
    sys.exit(main())
