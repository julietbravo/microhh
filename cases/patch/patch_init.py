import numpy as np
import matplotlib.pylab as pl

from microhh_tools import *     # Available in MICROHH_DIR/python/

# Read the name list (.ini file)
nl    = Read_namelist()
kmax  = nl.grid.ktot
zsize = nl.grid.zsize

# Initial profiles
# ----------------------------------------------
# Case settings
dthetadz = 0.003 # Potential temperature lapse rate (K/m)
lambda_q = 2000  # Moisture scale height (Stevens, 2007, (m))

z_canopy = 100   # Height of canopy (m) (only for swcanopy=1 in patch.ini) 
dz_0     = 10    # Grid spacing lowest grid level (m) (only for stretched grid)
dz_s     = 1.04  # Increase grid spacing per grid level (-) (only for stretched grid)

# Create empty arrays for vertical profiles
z    = np.zeros(kmax) # Full grid level (m)
th   = np.zeros(kmax) # liquid water potential temperature (K)
qt   = np.zeros(kmax) # Specific humitidy (g/kg)
u    = np.zeros(kmax) # u-component wind (m/s)
ug   = np.zeros(kmax) # u-component geostrophic wind (m/s)
acp  = np.zeros(kmax) # Reduction function canopy drag (0..1)

# Define vertical grid -> stretched
if  (False):
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

# Define vertical grid -> uniform
if (True):
    dz = zsize / kmax
    z  = np.linspace(0.5*dz, zsize-0.5*dz, kmax)

# Define the profiles:
for k in range(kmax):
    th[k] = 300. + dthetadz * z[k]
    qt[k] = 5e-3 * np.exp(-z[k] / lambda_q)
    u [k] = 0.
    ug[k] = 0.

    # Canopy drag reduction from 1 at surface to 0 at z_canopy:
    if (z[k] < z_canopy):
        acp [k] = 1-(z[k]/z_canopy)

# Write the data to *.prof file as input for MicroHH. Profiles which aren't specified
# (like e.g. `v` or `vg` in this case) are initialised at zero by the model.
proffile = open('patch.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s}\n'.format('z','thl','qt','u','ug','acp'))   # header
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E}\n'.format(z[k],th[k],qt[k],u[k],ug[k],acp[k]))
proffile.close()

# Time dependent surface data
# ----------------------------------------------
time = np.arange(0,36000.01,3600)
wthl = 0.1  * np.sin(np.pi * time / time.max())
wqt  = 1e-4 * np.sin(np.pi * time / time.max())

# Write to *.time file as input for MicroHH.
timefile = open('patch.time','w')
timefile.write('{0:^20s} {1:^20s} {2:^20s}\n'.format('t','sbot[thl]', 'sbot[qt]'))
for t in range(time.size):
    timefile.write('{0:1.14E} {1:1.14E} {2:1.14E}\n'.format(time[t], wthl[t], wqt[t]))
timefile.close()

