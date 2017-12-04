import numpy as np

# microhh_tools.py is available in the microhh/python directory,
# you need to copy it to this directory for this script to work
import microhh_tools as mht

# Settings sinusoidal hills
amplitude  = 0.00254
wavelength = 0.0508
z_offset   = 0.002

# Example of adding roughness
add_cubes    = False
patch_size   = wavelength / 5.
block_size   = patch_size / 3.
block_height = 0.5*amplitude

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

# Add cubes
if (add_cubes):
    from scipy.special import erf

    # Create blocks on flat surface
    xmod     = np.fmod(x, patch_size);
    z_block  = (0.5 - 0.5*erf(2.*(np.abs(2.*xmod - patch_size) - block_size) / 1e-9)) * block_height

    # Add to sinusoidal hills
    zIB[:,:] += z_block[None,:]

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
    ax=pl.subplot(111, aspect='equal')
    pl.plot(x, zIB[0,:], 'k-',  label='zIB')
    pl.legend(frameon=False, loc='best')
    pl.ylim(0,zsize)
    pl.xlabel('x (m)')
    pl.ylabel('z (m)')
