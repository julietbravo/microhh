import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc

pl.close('all')

float_type = 'f8'

# Get number of vertical levels and size from .ini file,
with open('canopy.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

# Define vertical grid and fields,
dz = zsize / kmax
z = np.arange(0.5*dz, zsize, dz)
zh = np.arange(0, zsize+0.1, dz)

u = np.zeros(kmax)
pad = np.zeros(kmax)

# Vertical extent canopy.
# This is not very flexible with the area densities defined below.
zc = 10.

# One-sided plant area density (m2 m-3).
a0 = 5.32e-2

# Plant area index (m2 m-2).
pai = 0.479

# Vertical profiles PAD with different smoothness at top.
a_sharp = np.zeros(kmax+1)
a_sharp[:10] = a0

a_twostep = np.zeros(kmax+1)
a_twostep[:9] = a0
a_twostep[9] = 0.5*a0

a_fourstep = np.zeros(kmax+1)
a_fourstep[:8] = a0
a_fourstep[8] = 7/8*a0
a_fourstep[9] = 4/8*a0
a_fourstep[10] = 1/8*a0

"""
# Comparison with Ouwersloot (2017, Fig. 1)
pl.figure()
pl.plot(a_sharp/a0, zh/zc, label='sharp')
pl.plot(a_twostep/a0, zh/zc, label='two step')
pl.plot(a_fourstep/a0, zh/zc, label='four step')
pl.ylim(0, 1.3)
pl.legend()
"""

# Switch between cases.
pad_h = a_fourstep

# Scale profile to prescribed PAI value.
pad_h *= pai / np.sum(pad_h * dz)

# Interpolate to full levels
#pad[:] = 0.5*(pad_h[1:] + pad_h[:-1])

# Constant velocity of (6,0,0)
u[:] = 6.

# Write NetCDF input.
nc_file = nc.Dataset('canopy_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

nc_file.createDimension('z', kmax)
nc_file.createDimension('zh', kmax+1)
nc_z = nc_file.createVariable('z' , float_type, ('z'))
nc_zh = nc_file.createVariable('zh' , float_type, ('zh'))

nc_group_init = nc_file.createGroup('init');
nc_u = nc_group_init.createVariable('u' , float_type, ('z'))
nc_padh = nc_group_init.createVariable('padh', float_type, ('zh'))

nc_z[:] = z[:]
nc_zh[:] = zh[:]
nc_u[:] = u[:]
nc_padh[:] = pad_h[:]

nc_file.close()
