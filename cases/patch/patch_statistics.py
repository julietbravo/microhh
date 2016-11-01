from microhh_tools import *

s  = Read_statistics('patch.default.0000000.nc')
sc = Read_statistics('patch.patch_high.0000000.nc')
sr = Read_statistics('patch.patch_low.0000000.nc')

t0 = 0
t1 = -1

pl.figure()
pl.subplot(221)
pl.plot(s. th[t0:t1,:].mean(axis=0), s. z, label='full')
pl.plot(sc.th[t0:t1,:].mean(axis=0), sc.z, label='city')
pl.plot(sr.th[t0:t1,:].mean(axis=0), sr.z, label='rural')
pl.legend(frameon=False)

pl.subplot(222)
pl.plot(s. u[t0:t1,:].mean(axis=0), s. z, label='full')
pl.plot(sc.u[t0:t1,:].mean(axis=0), sc.z, label='city')
pl.plot(sr.u[t0:t1,:].mean(axis=0), sr.z, label='rural')
pl.legend(frameon=False)

pl.subplot(223)
pl.plot(s. v[t0:t1,:].mean(axis=0), s. z, label='full')
pl.plot(sc.v[t0:t1,:].mean(axis=0), sc.z, label='city')
pl.plot(sr.v[t0:t1,:].mean(axis=0), sr.z, label='rural')
pl.legend(frameon=False)

pl.subplot(224)
pl.plot(s. w[t0:t1,:].mean(axis=0), s. zh, label='full')
pl.plot(sc.w[t0:t1,:].mean(axis=0), sc.zh, label='city')
pl.plot(sr.w[t0:t1,:].mean(axis=0), sr.zh, label='rural')
pl.legend(frameon=False)
