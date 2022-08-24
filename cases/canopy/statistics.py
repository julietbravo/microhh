import matplotlib.pyplot as pl
import netCDF4 as nc4
import xarray as xr
import numpy as np

pl.close('all')

def xr_read_all(f, groups=['default'], decode_times=True):
    # Read all NetCDF groups into a single Dataset.
    dss = [xr.open_dataset(f, decode_times=decode_times)]
    for group in groups:
        dss.append(xr.open_dataset(f, group=group, decode_times=decode_times))
    return xr.merge(dss)


s1 = xr_read_all('canopy.default.0000000.nc', decode_times=False)
s1m = s1.sel(time=slice(1500, 2400)).mean(dim='time')

# Friction velocity at canopy top:
ustar_zc = ((s1m.u_flux[10]**2 + s1m.v_flux[10]**2))**(1/4)

# Canopy height:
zc = 10.

pl.figure()
pl.plot(s1m.u/ustar_zc, s1m.z/zc)
pl.xlabel(r'$u/u*$ (-)')
pl.ylabel(r'$z/z_c$ (-)')
pl.xlim(0, 8)
pl.ylim(0, 3)

pl.figure()
pl.subplot(221)
pl.plot(np.sqrt(s1m.u_2) / ustar_zc, s1m.z/zc)
pl.xlabel(r'$\sigma_u/u*$ (-)')
pl.ylabel(r'$z/z_c$ (-)')
pl.ylim(0, 6)
pl.xlim(0, 2.5)

pl.subplot(222)
pl.plot(np.sqrt(s1m.v_2) / ustar_zc, s1m.z/zc)
pl.xlabel(r'$\sigma_v/u*$ (-)')
pl.ylabel(r'$z/z_c$ (-)')
pl.ylim(0, 6)
pl.xlim(0, 2)

pl.subplot(223)
pl.plot(np.sqrt(s1m.w_2) / ustar_zc, s1m.zh/zc)
pl.xlabel(r'$\sigma_w/u*$ (-)')
pl.ylabel(r'$z/z_c$ (-)')
pl.ylim(0, 6)
pl.xlim(0, 1.5)
