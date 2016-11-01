import numpy as np
from microhh_tools import *

# Read the name list (.ini file)
nl    = Read_namelist()
kmax  = nl.grid.ktot
zsize = nl.grid.zsize

# Case settings
dthetadz = 0.003 # Potential temperature lapse rate (K/m)
z_canopy = 200   # Height of canopy (m)
dz_0     = 10    # Grid spacing lowest grid level (m)
dz_s     = 1.04  # Increase grid spacing per grid level (-)

# Create empty arrays for vertical profiles
z    = np.zeros(kmax) # Full grid level (m)
th   = np.zeros(kmax) # potential temperature (K)
u    = np.zeros(kmax) # u-component wind (m/s)
ug   = np.zeros(kmax) # u-component geostrophic wind (m/s)
acp  = np.zeros(kmax) # Reduction function canopy drag (0..1)

# Define vertical grid
dz    = np.zeros(kmax)
zh    = np.zeros(kmax+1)
dz[0] = dz_0
for k in range(1, kmax):
    dz[k] = dz[k-1] * dz_s    
for k in range(1, kmax+1):
    zh[k] = zh[k-1] + dz[k-1] 
for k in range(kmax):
    z[k] = 0.5 * (zh[k] + zh[k+1])    
# Write the vertical grid size back to the namelist:
print('zsize = {}'.format(zh[-1]))
replace_namelist_var('zsize', zh[-1])

#dz = zsize / kmax
#z = np.linspace(0.5*dz, zsize-0.5*dz, kmax)

# Define the profiles:
for k in range(kmax):
    th  [k] = 300. + dthetadz*z[k]
    u   [k] = 10.
    ug  [k] = 10.

    # Canopy drag reduction from 1 at surface to 0 at z_canopy:
    if (z[k] < z_canopy):
        acp [k] = 1-(z[k]/z_canopy)

# Write the data to *.prof file as input for MicroHH. Profiles which aren't specified 
# (like e.g. `v` or `vg` in this case) are initialised at zero by the model.
proffile = open('patch.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s}\n'.format('z','th','u','ug','acp'))   # header
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E}\n'.format(z[k],th[k],u[k],ug[k],acp[k]))
proffile.close()
