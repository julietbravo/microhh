import matplotlib.pyplot as pl
import netCDF4 as nc4
import xarray as xr
import numpy as np
from cycler import cycler

pl.close('all')


def xr_read_all(f, groups=None):
    """
    Read all (or selection of) NetCDF groups from NetCDF
    file using xarray. If `groups=None`, all NetCDF groups
    are loaded.
    """
    nc = nc4.Dataset(f)
    if groups is None:
        groups = list(nc.groups)

    # Check of NetCDF file has meaningful time units.
    if nc.variables['time'].units == 'seconds since start':
        decode_times = False
    else:
        decode_times = True

    nc.close()

    # Read all groups into a single Dataset.
    dss = [xr.open_dataset(f, decode_times=decode_times)]
    for group in groups:
        dss.append(xr.open_dataset(f, group=group, decode_times=decode_times))
    return xr.merge(dss)


stat   = xr_read_all('weisman_klemp.default.0000000.nc', groups=['default', 'thermo'])
budget = xr_read_all('weisman_klemp.default.0000000.nc', groups=['sb06_budget'])


colors = ['tab:red', 'tab:blue', 'tab:green', 'tab:purple']
line_styles = (cycler(linestyle=['-', '--', ':', '-.']) * cycler(color=colors))

species = ['qv', 'qc', 'qi', 'qr', 'qs', 'qg', 'qh']

t = 60
pl.figure(figsize=(12,8))

ncol = 4
nrow = 2

sp = 1
for qx in species:
    ax=pl.subplot(nrow, ncol, sp); sp+=1
    ax.set_prop_cycle(line_styles)
    pl.title(qx, loc='left')
    for v in budget:
        if qx in v:
            label = v.replace(f'{qx}_', '')
            pl.plot(budget[v][t,:]*1e3, budget.z/1e3, label=label)
    pl.legend(fontsize=8)
    pl.xlabel('Tendency (g kg-1 s-1)')
    if sp not in [2,6]:
        ax.set_yticklabels([])
    else:
        pl.ylabel('z (km)')

pl.tight_layout()
