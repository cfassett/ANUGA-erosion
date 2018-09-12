#!/usr/bin/env python
import os, optparse, subprocess, sys, tempfile, ogr, osr, gdal

import copy
import pdb
import anuga
import numpy as np
import math
import DiamondSquare as DS
from suspendedtransport import suspension_erosion_operator
from bedloadtransport import shear_erosion_operator
from friction import friction_operator
from AngleofRepose import AoR_operator
from jezeroanalysisfunctions import finalstats

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
This program does lake erosion in ANUGA.  
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
    try:
        try:
            usage = "usage: jezero.py jezerostage*-1 grainsize-in-m \n  Note that 2260 flows but barely does anything, certainly not with a high grain size \n probably 2200 is too much.  note sign convention."
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("-m", dest="dm", help="Depth Multiplier for lake")
            parser.add_option("-i", dest="initbreachdepth", help="Initial breach depth for flood")
            (options, inargs) = parser.parse_args()
            if not inargs: parser.error("need parameters")
            jezerostage=-1.0*float(inargs[0])
            grainD=float(inargs[1])
            # Before doing anything else, set gravity to Mars
            anuga.g=3.711


            #domain parameters
            poly_highres =anuga.read_polygon('PolyHi.csv')
            poly_bound = anuga.read_polygon('PolyBoundBiggerPts.csv')
            base_scale = 25600  #160 m  #HIRES REGION
            lowres = 10000*base_scale   # 10 km 
            
            boundary_tags={'exterior': [0]}

            # Define list of interior regions with associated resolutions
            interior_regions = [[poly_highres,base_scale]]



            meshname='jezero.msh'
            mesh = anuga.create_mesh_from_regions(poly_bound, boundary_tags=boundary_tags, maximum_triangle_area=lowres, interior_regions=interior_regions, verbose=True, filename=meshname, fail_if_polygons_outside=False)
            evolved_quantities =  ['stage', 'xmomentum', 'ymomentum', 'concentration','sedvolx','sedvoly']
            domain = anuga.Domain(meshname, use_cache=False, evolved_quantities = evolved_quantities)
            domain.g = anuga.g  # make sure the domain inherits the package's gravity.
            
            domain.set_quantity("elevation",filename="jbigger.pts",use_cache=False, verbose=True)
            domain.smooth=True

            name="jezero"+str(jezerostage)+"_"+str(grainD*1000)+"mm"            
            domain.set_name(name)            
            initialdomain=copy.deepcopy(domain)
            
            # Choose between DE0 (less accurate), DE1 (more accurate, slower), DE2 (even more accurate, slower still)
            domain.set_flow_algorithm('DE0')
            domain.set_CFL(cfl=1)
            domain.set_minimum_allowed_height(0.01)  #raise the default minimum depth from 1 mm to 1 cm to speed the simulation (for computation)
            domain.set_minimum_storable_height(0.01) #raise the default minimum depth from 1 mm to 1 cm to speed the simulation (for storage/visualization)
            domain.set_maximum_allowed_speed(1)      #maximum particle speed that is allowed in water shallower than minimum_allowed_height (back of the envelope suggests 1m/s should be below tau_crit)

            np.seterr(invalid='ignore')
            elevs=np.load('XelevC.npy')
            stage=np.load('XstageC.npy')
            conc=np.load('XconcC.npy')
            xm=np.load('XxmC.npy')
            ym=np.load('XymC.npy')
            domain.set_quantity('elevation', elevs, location='centroids') 
            domain.set_quantity('concentration', conc, location='centroids')        # Start with no esd
            domain.set_quantity('stage', stage, location='centroids')
            domain.set_quantity('xmomentum', xm, location='centroids')
            domain.set_quantity('ymomentum', ym, location='centroids')

            elevs=np.load('XelevV.npy')
            stage=np.load('XstageV.npy')
            conc=np.load('XconcV.npy')
            xm=np.load('XxmV.npy')
            ym=np.load('XymV.npy')
            domain.set_quantity('elevation', elevs, location='vertices') 
            domain.set_quantity('concentration', conc, location='vertices')        # Start with no esd
            domain.set_quantity('stage', stage, location='vertices')
            domain.set_quantity('xmomentum', xm, location='vertices')
            domain.set_quantity('ymomentum', ym, location='vertices')

            voltotal=np.sum(domain.quantities['elevation'].centroid_values*domain.areas)

            domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 2})

            np.seterr(invalid='warn')
            Br=anuga.Reflective_boundary(domain)
            domain.set_boundary({'exterior': Br})

            operateone = friction_operator(domain,poly_bound,grainD=grainD, gravity=anuga.g)
            operatetwo = suspension_erosion_operator(domain, poly_bound,grainD=grainD, gravity=anuga.g)
            operatethree = shear_erosion_operator(domain,poly_bound,grainD=grainD, gravity=anuga.g)
            operatefour = AoR_operator(domain, poly_bound)  #no grainsize or gravity dependence for the angle of repose change

            initial=domain.get_quantity('elevation').get_values(interpolation_points=[[19133.08548, 1097174.738]], location='centroids')
            ystep=60
            ftime=120
            count=0
            volthresh=1000
            for t in domain.evolve(yieldstep=ystep, finaltime=ftime):
                print domain.timestepping_statistics()
                print 'xmom:'+str(domain.get_quantity('xmomentum').get_values(interpolation_points=[[19133.08548, 1097174.738]], location='centroids'))
                breacherosion=domain.get_quantity('elevation').get_values(interpolation_points=[[19133.08548, 1097174.738]], location='centroids')-initial
                print 'erosion: '+str(breacherosion)


                volcurr=np.sum(domain.quantities['elevation'].centroid_values*domain.areas)
                volsed=np.sum(domain.quantities['concentration'].centroid_values*domain.quantities['height'].centroid_values*domain.areas)            
                conservation=(volcurr+volsed-voltotal)/voltotal
                print 'conservation: '+'{:.8%}'.format(conservation)
                if (volsed<volthresh)&(count>0):
                    print "No sediment moving...ending..."
                    break
                count=count+1
                
                
            polyline=[[40000,1090000],[15000,1120000]]
            time, Q = anuga.get_flow_through_cross_section(name+'.sww',polyline)
            print Q

            initname="initial_"+name+".asc"
            finname="final_"+name+".asc"
            fluxname="flux_"+name+".txt"
            np.savetxt(fluxname,Q)
            anuga.sww2dem(name+'.sww',initname,reduction=0)
            anuga.sww2dem(name+'.sww',finname,reduction=(count-1))
            finalstats(domain,initialdomain)

        except optparse.OptionError, msg:
            raise Usage(msg)

    except Usage, err:
        print >>sys.stderr, err.msg
        # print >>sys.stderr, "for help use --help"
        return 2

if __name__ == "__main__":
    sys.exit(main())
