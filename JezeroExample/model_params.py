# This file controls all of the parameters for a simulation
# The goal of centralizing these to avoid inconsistencies between operators


# DOMAIN PARAMETERS
gravity=3.711      # Mars Gravity       m/s2                  
dthresh=0.01       # Depth Threshold    m  {Anuga numerical  parameter}

# FRICTION SCHEME
frictionscheme="wilsonetal"   # "wilsonetal", "larsenandlamb", "constant"   -->defaults to constant
constantn=0.0545        # Used for constant, default in unwet cells for other schemes

                        #  see Larsen and Lamb, methods, eq. 4.  These parameters only matter if the friction scheme is Larsen and Lamb.
r_d=2.0                 # Hydraulic roughness scaling parameter
r_br=2.0                # ""
sigma_br=5.0            # bed roughness sigma (m) 
ks=r_d*r_br*sigma_br    # k_s in Larsen and Lamb

# ANGLE OF REPOSE PARAMETER
angleofrepose=35.0  # Degrees

# FLUID/SEDIMENT PARAMETERS
grainD=0.01        # Grain Size         m   MAIN PARAMETER FOR CONTROLLING RESISTANCE
                   #                        TO FLOW

Wd=1000.0          # Water Density      kg/m3
Sd=3000.0          # Sediment Density   kg/m3 
porosity=0.2       # Sediment porosity  %

# SEDIMENT TRANSPORT PARAMETERS
# GENERAL/NUMERICAL
maxrate=5.0        # max erosion m/s.  This helps maintain numerical stability by 
                   #    preventing single triangles from deepening aphysically.

# BEDLOAD
bed_tau_crit=0.04  # bedload tau_critical 
m=1.5              # exponent on dimensionless shear stress
Ke=4.0             # leading coefficient on qb*  From Wong and Parker 2006 eq 24

# SUSPENDED 
gamma0   = 0.0024       # suspended sediment gamma0 (Smith and McLean)
susp_tau_crit = 0.04   # Tau_crit for entrainment
mu       = 8.9e-4      # Dynamic viscosity of water
maxconc  = 0.30        # artificial limit on sediment conc. by volume, expressed as %  Using 30%.  
                            # (Note this is a huge number by mass, especially since I'm not back-calculating the effect on the flow!)

# note... there are a bunch of additional hardwired suspended sediment parameters, particularly the Rouse profile.  
# Also, settling velocities in the suspended sediment scheme are from Dietrich 1982.  They are non-dimensionalized but have hardwired coefficients.
