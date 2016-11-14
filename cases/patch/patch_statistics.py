import matplotlib.pylab as pl
import numpy as np

pl.ion()

from microhh_tools import *   # Available in the MICROHH_DIR/python directory

s  = Read_statistics('patch.default.0000000.nc')
sc = Read_statistics('patch.patch_high.0000000.nc')
sr = Read_statistics('patch.patch_low.0000000.nc')

t = -1

pl.figure()
pl.subplot(221)
pl.plot(s. thl[t,:], s. z, label='full')
pl.plot(sc.thl[t,:], sc.z, label='city')
pl.plot(sr.thl[t,:], sr.z, label='rural')
pl.legend(frameon=False)

pl.subplot(222)
pl.plot(s. qt[t,:]*1000, s. z, label='full')
pl.plot(sc.qt[t,:]*1000, sc.z, label='city')
pl.plot(sr.qt[t,:]*1000, sr.z, label='rural')
pl.legend(frameon=False)

pl.subplot(223)
pl.plot(s. u[t,:], s. z, label='full')
pl.plot(sc.u[t,:], sc.z, label='city')
pl.plot(sr.u[t,:], sr.z, label='rural')
pl.legend(frameon=False)

pl.subplot(224)
pl.plot(s. v[t,:], s. z, label='full')
pl.plot(sc.v[t,:], sc.z, label='city')
pl.plot(sr.v[t,:], sr.z, label='rural')
pl.legend(frameon=False)

pl.figure()
T = s.thl * exner(s.ph)
pl.plot(s.qt[0,:] / qsat(s.ph[0,:], T[0,:]), s.z)
