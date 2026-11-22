#!/usr/bin/env python
from __future__ import print_function

import optparse, sys
import copy
import anuga
import numpy as np
import math
from suspendedtransport_operator import suspendedtransport_operator
from bedloadtransport_operator import bedloadtransport_operator
from friction_operator import friction_operator
from AngleofRepose_operator import AoR_operator
from jezeroanalysisfunctions import finalstats
from model_params import grainD, gravity, constantn


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
            usage = "usage: jezero.py -1*initialstage-in-m \n  Note that -2260 flows but barely does anything.  So range for jezero is ~2260 to 2200. \n"
            parser = optparse.OptionParser(usage=usage)
    
            (options, inargs) = parser.parse_args()
            if not inargs: parser.error("need parameters")
            jezerostage=-1.0*float(inargs[0])            
            # Before doing anything else, set gravity to Mars
            anuga.g=gravity


            #domain parameters
            poly_highres =anuga.read_polygon('JezeroData/PolyHi.csv')
            poly_bound = anuga.read_polygon('JezeroData/PolyBoundBiggerPts.csv')
            base_scale = 25600  #160 m  #HIRES REGION
            lowres = 10000*base_scale   # 10 km 
            
            boundary_tags={'exterior': [0]}

            # Define list of interior regions with associated resolutions
            interior_regions = [[poly_highres,base_scale]]



            meshname='JezeroData/jezero.msh'
            mesh = anuga.create_mesh_from_regions(poly_bound, boundary_tags=boundary_tags, maximum_triangle_area=lowres, interior_regions=interior_regions, verbose=True, filename=meshname, fail_if_polygons_outside=False)
            evolved_quantities =  ['stage', 'xmomentum', 'ymomentum', 'concentration','sedvolx','sedvoly']
            domain = anuga.Domain(meshname, use_cache=False, evolved_quantities = evolved_quantities)
            domain.g = anuga.g  # make sure the domain inherits the package's gravity.
            
            domain.set_quantity("elevation",filename="JezeroData/jbigger.pts",use_cache=False, verbose=True)
            domain.smooth=True

            name="jezero"+str(jezerostage)+"_"+str(grainD*1000.0)+"mm"
            domain.set_name(name)            
            
            # Choose between DE0 (less accurate), DE1 (more accurate, slower), DE2 (even more accurate, slower still)
            domain.set_flow_algorithm('DE0')
            domain.set_CFL(cfl=1.0)
            domain.set_minimum_allowed_height(0.01)  #raise the default minimum depth from 1 mm to 1 cm to speed the simulation (for computation)
            domain.set_minimum_storable_height(0.01) #raise the default minimum depth from 1 mm to 1 cm to speed the simulation (for storage/visualization)
            domain.set_maximum_allowed_speed(1.0)      #maximum particle speed that is allowed in water shallower than minimum_allowed_height (back of the envelope suggests 1m/s should be below tau_crit)

 
            np.seterr(invalid='ignore')
            voltotal=np.sum(domain.quantities['elevation'].centroid_values*domain.areas)
            domain.set_quantity('concentration', 0.00)        # Start with no esd
            domain.set_quantity('friction', constantn)      # Constant Manning friction.  Only used to initialize domain (calculated every time step).
            domain.set_quantity('stage', -10000.0)
            stagepolygon=anuga.read_polygon('JezeroData/stage.csv')
            domain.set_quantity('stage', jezerostage, polygon=stagepolygon)
            domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 2})

            np.seterr(invalid='warn')
            Br=anuga.Reflective_boundary(domain)
            domain.set_boundary({'exterior': Br})

#           fixed_inflow = anuga.Inlet_operator(domain, stagepolygon, Q=1000.0)
            operatezero = suspendedtransport_operator(domain)
            operateone = bedloadtransport_operator(domain)
            operatetwo = friction_operator(domain)
            operatethree = AoR_operator(domain)

            breachpoint=[[19080.340,1097556.781]]
            initialdomain=copy.deepcopy(domain)
            initial=domain.get_quantity('elevation').get_values(interpolation_points=breachpoint, location='centroids')
            print (str(initial))
            ystep=300.0
            ftime=86400.0*5.0
            polyline=[[40000,1090000],[15000,1120000]]
            count=0
            for t in domain.evolve(yieldstep=ystep, finaltime=ftime):
                print (domain.timestepping_statistics())
                print ('xmom:'+str(domain.get_quantity('xmomentum').get_values(interpolation_points=breachpoint, location='centroids')))
                breacherosion=domain.get_quantity('elevation').get_values(interpolation_points=breachpoint, location='centroids')-initial
                print ('erosion: '+str(breacherosion))
                volcurr=np.sum(domain.quantities['elevation'].centroid_values*domain.areas)
                volsed=np.sum(domain.quantities['concentration'].centroid_values*domain.quantities['height'].centroid_values*domain.areas)            
                conservation=(volcurr+volsed-voltotal)/voltotal
                print ('conservation: '+'{:.8%}'.format(conservation))
                count=count+1
                
            time, Q = anuga.get_flow_through_cross_section(name+'.sww',polyline)
            print (Q)

            initname="initial_"+name+".asc"
            finname="final_"+name+".asc"
            fluxname="flux_"+name+".txt"
            np.savetxt(fluxname,Q)
            anuga.sww2dem(name+'.sww',initname,reduction=0)
            anuga.sww2dem(name+'.sww',finname,reduction=(count-1))
            finalstats(domain,initialdomain)
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
