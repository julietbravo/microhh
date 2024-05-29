import matplotlib.pyplot as pl
import xarray as xr
import numpy as np
from cycler import cycler

pl.close('all')

nc = xr.open_dataset('weisman_klemp.timing.micro.0000000.nc', decode_times=False)

"""
Time series of all processes.
"""
pl.figure(figsize=(7,5))

colors = [pl.cm.Set1(i) for i in range(pl.cm.Set1.N)]
line_styles = (cycler(linestyle=['-', '--', ':', '-.']) * cycler(color=colors))
pl.gca().set_prop_cycle(line_styles)

for v in nc.variables:
    if v != 'time' and 'mean' in v and 'total' not in v:
        pl.plot(nc.time, nc[v], label=v)
pl.legend(ncol=3, fontsize=7)
pl.ylabel('Mean execution time (s)')
pl.xlabel('Simulation time (s)')
pl.tight_layout()


"""
Sort processes on mean execution time.
"""
means = {}
total = 0.
for v in nc.variables:
    if v != 'time' and 'mean' in v and 'total' not in v:
        means[v] = nc[v].values.mean()
        total += means[v]

means = {k: v for k, v in sorted(means.items(), key=lambda item: item[1])}

pl.figure()
ax=pl.subplot()
i = 0
for process,timing in means.items():
    process = process.replace('_mean', '')
    perc = timing/total*100
    pl.scatter(i, perc, marker='+', color='k')
    pl.text(i, perc+0.5, process, rotation=90, size=8, ha='center', va='bottom')
    i += 1
ax.set_xticklabels([])
pl.ylabel('Execution time (%)')
pl.tight_layout()
