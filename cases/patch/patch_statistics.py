import matplotlib.pylab as pl
import numpy as np

pl.close('all')
pl.ion()

import microhh_tools as mhh    # Available in the MICROHH_DIR/python directory

s  = mhh.Read_statistics('patch.default.0000000.nc')

t = -1

# Vertical profiles
pl.figure()
pl.subplot(321)
pl.plot(s. thl[t,:], s. z)
pl.legend(frameon=False)
pl.subplot(322)
pl.plot(s. qt[t,:]*1000, s. z)
pl.subplot(323)
pl.plot(s. ql[t,:]*1000, s. z)
pl.subplot(324)
pl.plot(s. u[t,:], s. z)
pl.subplot(325)
pl.plot(s. v[t,:], s. z)

# Relative humidity
pl.figure()
pl.subplot(111)
T = s.thl * mhh.exner(s.ph)
pl.plot(s.qt[0, :] / mhh.qsat(s.ph[0, :], T[0, :]), s.z)
pl.plot(s.qt[-1,:] / mhh.qsat(s.ph[-1,:], T[-1,:]), s.z)
pl.ylabel('z (m)')
pl.xlabel('RH (-)')

# Surface fluxes
H  = s.thlflux[:,0] * s.rhorefh[0] * mhh.cp
LE = s.qtflux [:,0] * s.rhorefh[0] * mhh.Lv

pl.figure()
pl.subplot(131)
pl.plot(s.t/3600., H)
pl.xlabel('time (h)')
pl.ylabel('H (W/m2)')

pl.subplot(132)
pl.plot(s.t/3600., LE)
pl.xlabel('time (h)')
pl.ylabel('LE (W/m2)')

pl.subplot(133)
pl.plot(s.t/3600., LE/(H+LE))
pl.xlabel('time (h)')
pl.ylabel('EF (-)')
