from __future__ import print_function

import anuga, copy
import numpy as np
import math
from scipy import spatial
from suspendedtransport_operator import suspendedtransport_operator
from bedloadtransport_operator import bedloadtransport_operator
from friction_operator import friction_operator
from AngleofRepose_operator import AoR_operator
from model_params import gravity, grainD, constantn


# Before doing anything else, set gravity 
anuga.g=gravity

# Set grainsize parameter for the model 
# Presuming that the walnut shell and sand is both 500 microns (coarse).  Using a weighted average sediment density of
# 0.8*1400.0+ 0.2*2700.0= 1660.0

#grainD=0.0005  see model_params.py

name="CoholichLCLR"

#Geometric parameters for setting up the lake and how full it is (ref in topography, stage)
domainpolygon=[[0,0],[5,0],[5,0.79],[0,0.79],[0,0]]
outletregion=[[0,0],[4,0],[4,0.79],[0,0.79],[0,0]]
boundary_tags={'top': [0], 'right': [1],'bottom': [2], 'left': [3]}
base_resolution=0.0004 #triangle areas.  smaller is higher res.
low_res=1000*base_resolution
interior_regions = [[outletregion, base_resolution]]

#y
craterwidth=0.79
#x
craterlength=0.35
rimheight=0.12
rimwidth=0.14
extslope=0.04
buffslope=0.2
lakeslope=0.4
initbreachwidth=0.08
initbreachdepth=0.08

lakeregion=[[0,0],[0.28,0],[0.28,craterwidth],[0,craterwidth],[0,0]]
nextregion=[[0.28,0],[0.28,craterwidth],[4,craterwidth],[4,0],[0.28,0]]
endregion=[[4,craterwidth],[5,craterwidth],[5,0],[4,0],[4,craterwidth]]


meshname='lake.msh'
m = anuga.create_mesh_from_regions(domainpolygon, boundary_tags, maximum_triangle_area=low_res,interior_regions=interior_regions,filename=meshname, use_cache=False)
evolved_quantities =  ['stage', 'xmomentum', 'ymomentum', 'concentration', 'sedvolx', 'sedvoly']
#evolved_quantities =  ['stage', 'xmomentum', 'ymomentum', 'concentration']
domain = anuga.Domain(meshname, use_cache=False, evolved_quantities = evolved_quantities)
domain.g = anuga.g  # make sure the domain inherits the package's gravity.

print ('Number of triangles = ', len(domain))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
    
def topography(x, y):
    maxx = np.max(x)
    miny, maxy = np.min(y), np.max(y)
    
    #sets slope of flume surface
    z=(4.0-x)*(extslope)
    #sets up lake (lakeslope is to prevent runup to rim being way too steep.)
    z[x<craterlength]=z[x<craterlength]-(x[x<craterlength])*lakeslope
    
    #adds rim
    z=np.where(np.logical_and(x>=craterlength,x<=craterlength+rimwidth),z+rimheight,z)
      
    #notches rim
    mask1 = (x>=craterlength)
    mask2 = (x<=craterlength+rimwidth)
    mask3 = (y>=((craterwidth/2.0-(initbreachwidth/2.0))))  
    mask4 = (y<=((craterwidth/2.0+(initbreachwidth/2.0))))
    #mask5 = (z>-initbreachdepth)
    #mask=mask1 & mask2 & mask3 & mask4 & mask5
    mask=mask1 & mask2 & mask3 & mask4
    z[mask]=z[mask]-initbreachdepth
        
    # add gaussian noise, with zero mean and 0.0005 m standard deviation (a grain width, arbitrary).
    z=z+np.random.normal(0,0.0005,z.shape)
    
    # set up edges as buffers
    maskx = (x>=craterlength)
    masky1 = (y>=0.7) & maskx
    z[masky1]=z[masky1]+buffslope*(y[masky1]-0.70)  
    masky2 = (y<=0.09) & maskx
    z[masky2]=z[masky2]+buffslope*(0.09-y[masky2])
    
    z[z<0.0]=0.0

    return z

def initialstage(x,y):
    istage=0.15
    istage=np.where(x>=craterlength,-1.0,istage)
    #istage=0.0
    return istage


# Choose between DE0 (less accurate), DE1 (more accurate, slower), DE2 (even more accurate, slower still)
domain.set_flow_algorithm('DE2')
domain.set_CFL(cfl=1.0)

np.seterr(invalid='ignore')
domain.set_quantity('elevation', topography, location='centroids') # elevation is a function
topototal=np.sum(domain.quantities['elevation'].centroid_values)
domain.set_quantity('concentration', 0.00)        # Start with no esd
domain.set_quantity('friction', constantn)      # Constant Manning friction.  
#domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 2, 'concentration': 2,'sedvolx':2,'sedvoly':2})
domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 2, 'concentration': 2})

#relax topography before the model run.
Relaxer=AoR_operator(domain)
for i in range(100):
     Relaxer.AofR()
domain.set_quantity('stage', initialstage, location='centroids')
domain.set_quantity('concentration', expression='elevation')   #use this solely for holding initial elevation in the unused concentration quantity for post-processing.
np.seterr(invalid='warn')

# Setup boundary conditions

Br=anuga.Reflective_boundary(domain)
#Bt=anuga.Transmissive_boundary(domain)
#Bd=anuga.Dirichlet_boundary([-10000,0,0])

domain.set_boundary({'left': Br, 'right': Br,'top':Br,'bottom':Br})
#domain.set_store_vertices_uniquely(True)

center=(craterlength/2.0,craterwidth/2.0)
radius=(0.1)
region0 = anuga.Region(domain, center=center, radius=radius)
outletregion0=anuga.Region(domain, polygon=outletregion)
lakeregion0=anuga.Region(domain, polygon=lakeregion)
nextregion0=anuga.Region(domain, polygon=nextregion)
endregion0=anuga.Region(domain, polygon=endregion)


#operatezero = suspension_erosion_operator(domain, outletregion,grainD=grainD, gravity=anuga.g)    #NO SUSPENSION IN THESE EXPERIMENTS
operateone = bedloadtransport_operator(domain)
#operatetwo = friction_operator(domain)   #FRICTION IS CONSTANT
operatethree = AoR_operator(domain)  #no grainsize or gravity dependence for the angle of repose change
#operatefour = stagelimiter_operator(domain,polygon=endregion)   #probably a bad idea
#pumprate=1L/min = 0.001m^3/min = 1.667e-5 m3/s
fixed_inflow = anuga.Inlet_operator(domain, region0 , Q=50.0*0.00001667)
#Add water to domain in lake



count=0
ystep=1.0
ftime=80.0
for t in domain.evolve(yieldstep=ystep, finaltime=ftime):
    print (domain.timestepping_statistics())
    volume = (domain.quantities['stage'].get_integral()-domain.quantities['elevation'].get_integral())
    volume2 = domain.quantities['stage'].get_integral(region=lakeregion0)-domain.quantities['elevation'].get_integral(region=lakeregion0)
    volume3 = domain.quantities['stage'].get_integral(region=nextregion0)-domain.quantities['elevation'].get_integral(region=nextregion0)
    volume4 = domain.quantities['stage'].get_integral(region=endregion0)-domain.quantities['elevation'].get_integral(region=endregion0)
    print (str(t)+","+str(volume)+","+str(volume2)+","+str(volume3)+","+str(volume4))
    if t==21: fixed_inflow.set_Q(1.0*0.00001667)
    if t==28: fixed_inflow.set_Q(0.0)
    print (fixed_inflow.get_Q())

    count=count+1
    
