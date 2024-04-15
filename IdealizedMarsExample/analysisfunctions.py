""" -----------------------------------------------
Analysis functions
Fassett implementation
======= """
import numpy as np
import anuga



def finalstats(domain, initialdomain, lakeradius, outfile='outstats.csv'):
    crossspacing=100
    acrossthebreach=np.arange(0.8*lakeradius,1.2*lakeradius,crossspacing)
    breacherosion=0
    areasum=0
    for y in acrossthebreach:
        ybreacherosion=domain.get_quantity('elevation').get_values(interpolation_points=[[2*lakeradius,y]], location='centroids')-initialdomain.get_quantity('elevation').get_values(interpolation_points=[[2*lakeradius,y]], location='centroids')
        if ybreacherosion<0: areasum=areasum+float(ybreacherosion)*crossspacing    #rectangular integration  :)
        if ybreacherosion<breacherosion: breacherosion=float(ybreacherosion)
    
    areasum=-areasum  #because negative areas seem weird
    
    elevdiff=domain.quantities['elevation'].centroid_values-initialdomain.quantities['elevation'].centroid_values  #negative values are eroded, positive values are deposited
    finalremovedvolume=-np.sum(elevdiff[elevdiff<0]*domain.areas[elevdiff<0])
    
    transportedmask=domain.get_centroid_coordinates()[:,0]>2*lakeradius
    heights=(domain.quantities['stage'].centroid_values-domain.quantities['elevation'].centroid_values)
    transportedvolume=(np.where(heights>0,heights,0)*domain.areas)[transportedmask].sum()

    
    print ('Breach Erosion Depth (m): '+str(-breacherosion))
    print ('Cross-sectional Area (m^2): '+str(areasum))
    print ('Outlet volume (m^3): '+str(finalremovedvolume))
    print ('Drained volume (m^3):' +str(transportedvolume))
    with open(outfile, 'a') as the_file:
        the_file.write(domain.get_name()+','+str(-breacherosion)+','+str(areasum)+','+str(finalremovedvolume)+','+str(transportedvolume)+'\n')
