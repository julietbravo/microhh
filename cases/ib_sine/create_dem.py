import numpy as np

# microhh_tools.py is available in the microhh/python directory,
# you need to copy it to this directory for this script to work
import microhh_tools as mht

# Settings sinusoidal hills
amplitude  = 0.00254
wavelength = 0.0508
z_offset   = 0.002

# Read namelist
nl = mht.Read_namelist()

# Short cuts
itot  = nl['grid']['itot']
jtot  = nl['grid']['jtot']
xsize = nl['grid']['xsize']
ysize = nl['grid']['ysize']
zsize = nl['grid']['zsize']

# Horizontal grid (x,y location at grid center)
dx = xsize / itot
dy = ysize / jtot
x  = np.arange(0.5*dx, xsize, dx)
y  = np.arange(0.5*dy, ysize, dy)

# 1D IB function
zIB_x = z_offset + amplitude + amplitude * np.sin(2*np.pi*x/wavelength)

# Expand in y-direction to 2D field
zIB = np.zeros((jtot, itot), dtype=np.float)
zIB[:,:] = zIB_x[None,:]

# Save the MicroHH dem.0000000 input file
mht.write_restart_file('dem.0000000', zIB[np.newaxis,:,:], itot, jtot, 1)

# Plot for debugging
if (True):
    import matplotlib.pylab as pl
    pl.close('all')

    pl.figure()
    pl.pcolormesh(x, y, zIB, cmap=pl.cm.RdBu_r)
    pl.colorbar()
    pl.xlabel('x (m)')
    pl.ylabel('y (m)')

    pl.figure()
    pl.plot(x, zIB[0,:], 'k-',  label='zIB')
    pl.legend(frameon=False, loc='best')
    pl.ylim(0,zsize)
    pl.xlabel('x (m)')
    pl.ylabel('z (m)')
