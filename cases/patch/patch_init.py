import numpy as np
import matplotlib.pylab as pl

import microhh_tools as mhh     # Available in MICROHH_DIR/python/

# Read the name list (.ini file)
nl    = mhh.Read_namelist()
kmax  = nl.grid.ktot
zsize = nl.grid.zsize

# --------------------------------------
# Create (stretched) vertical grid:
z     = np.zeros(kmax)      # Grid height at cell center
zh    = np.zeros(kmax+1)    # Grid height at cell edges

# For ktot=128
#dz0   = 10.                 # Vertical grid spacing near surface
#dzs   = 1.02                # Increase grid spacing per level with height
#dzmax = 40                  # Maximum vertical grid spacing

# For quick tests without patches, using ktot=64
dz0   = 35.                 # Vertical grid spacing near surface
dzs   = 1.035               # Increase grid spacing per level with height
dzmax = 60                  # Maximum vertical grid spacing

# Calculate grid at cell edges:
dz = dz0
for k in range(1,kmax+1):
    zh[k] = zh[k-1] + dz
    dz   *= dzs
    dz    = min(dz, dzmax)

# Calculate grid at cell center:
for k in range(kmax):
    z[k]  = 0.5 * (zh[k] + zh[k+1])

# Write the vertical domain size back the the .ini file:
print('zsize = {}'.format(zh[-1]))
mhh.replace_namelist_var('zsize', zh[-1])

# --------------------------------------
# Define initial vertical profiles:
th   = np.zeros(kmax)       # liquid water potential temperature (K)
qt   = np.zeros(kmax)       # Specific humitidy (g/kg)
u    = np.zeros(kmax)       # u-component wind (m/s)
ug   = np.zeros(kmax)       # u-component geostrophic wind (m/s)
acp  = np.zeros(kmax)       # Reduction function canopy drag (0..1)

# Thermodynamic quantities:
th[:]  = 290 + 0.0045 * z
qt[:]  = 11e-3 * np.exp(-z/2000)

# Wind:
u[:]   = 0.
ug[:]  = 0.

# Canopy drag over the city:
z_canopy = 100.
for k in range(kmax):
    if (z[k] < z_canopy):
        acp[k] = 1-(z[k]/z_canopy)

# --------------------------------------
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
thls = th[0] + 10 * np.sin(np.pi * time / time.max())
qts  = 0.6 * mhh.qsat(1e5, thls)        # assuming T=theta at surface (p=pref)

# Write to *.time file as input for MicroHH.
timefile = open('patch.time','w')
timefile.write('{0:^20s} {1:^20s} {2:^20s}\n'.format('t','sbot[thl]', 'sbot[qt]'))
for t in range(time.size):
    timefile.write('{0:1.14E} {1:1.14E} {2:1.14E}\n'.format(time[t], thls[t], qts[t]))
timefile.close()

if (True):
    pl.close('all')

    pl.figure()
    pl.subplot(321)
    pl.title('Vertical grid', loc='left')
    pl.plot(z, zh[1:]-zh[:-1], '-x')
    pl.xlabel('z (m)')
    pl.ylabel('dz (m)')

    pl.subplot(322)
    pl.plot(th, z, '-x')
    pl.xlabel('th (K)')
    pl.ylabel('z (m)')

    pl.subplot(323)
    pl.plot(qt*1000, z, '-x')
    pl.xlabel('qt (g/kg)')
    pl.ylabel('z (m)')

    pl.subplot(324)
    pl.plot(acp, z, '-x')
    pl.xlabel('acp (-)')
    pl.ylabel('z (m)')

    #pl.subplot(325)
    #pl.plot(time/3600, wthl, '-x')
    #pl.xlabel('time (h)')
    #pl.ylabel('wthl (K/m/s)')

    #pl.subplot(326)
    #pl.plot(time/3600, wqt*1000., '-x')
    #pl.xlabel('time (h)')
    #pl.ylabel('wqt (g/kg/m/s)')
