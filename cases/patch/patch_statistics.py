import matplotlib.pylab as pl
import numpy as np

pl.close('all')
pl.ion()

import microhh_tools as mhh    # Available in the MICROHH_DIR/python directory

# Read the statistic files
s  = mhh.Read_statistics('patch.default.0000000.nc')
sr = mhh.Read_statistics('patch.patch_low.0000000.nc')
su = mhh.Read_statistics('patch.patch_high.0000000.nc')

# Calculate energetic surface fluxes for all statistics files
for case in [s, sr, su]:
    case.H  = case.thlflux[:,0] * case.rhorefh[0] * mhh.cp
    case.LE = case.qtflux [:,0] * case.rhorefh[0] * mhh.Lv

t = -1

# Vertical profiles
pl.figure()
pl.subplot(321)
pl.plot(s. thl[t,:], s. z, label='total')
pl.plot(sr.thl[t,:], sr.z, label='rural')
pl.plot(su.thl[t,:], su.z, label='urban')
pl.legend(frameon=False)
pl.xlabel('thl (K)')
pl.ylabel('z (m)')

pl.subplot(322)
pl.plot(s. qt[t,:]*1000, s. z)
pl.plot(sr.qt[t,:]*1000, sr.z)
pl.plot(su.qt[t,:]*1000, su.z)
pl.xlabel('qt (g/kg)')
pl.ylabel('z (m)')

pl.subplot(323)
pl.plot(s. ql[t,:]*1000, s. z)
pl.plot(sr.ql[t,:]*1000, sr.z)
pl.plot(su.ql[t,:]*1000, su.z)
pl.xlabel('ql (g/kg)')
pl.ylabel('z (m)')

pl.subplot(324)
pl.plot(s. u[t,:], s. z)
pl.plot(sr.u[t,:], sr.z)
pl.plot(su.u[t,:], su.z)
pl.xlabel('u (m/s)')
pl.ylabel('z (m)')

pl.subplot(325)
pl.plot(s. v[t,:], s. z)
pl.plot(sr.v[t,:], sr.z)
pl.plot(su.v[t,:], su.z)
pl.xlabel('v (m/s)')
pl.ylabel('z (m)')

# Surface fluxes
pl.figure()
pl.subplot(221)
pl.plot(s.t /3600., s.H,  label='total')
pl.plot(sr.t/3600., sr.H, label='rural')
pl.plot(su.t/3600., su.H, label='urban')
pl.xlabel('time (h)')
pl.ylabel('H (W/m2)')
pl.legend(frameon=False, loc='best')

pl.subplot(222)
pl.plot(s.t /3600., s.LE,  label='total')
pl.plot(sr.t/3600., sr.LE, label='rural')
pl.plot(su.t/3600., su.LE, label='urban')
pl.xlabel('time (h)')
pl.ylabel('LE (W/m2)')

pl.subplot(223)
pl.plot(s.t/3600., s.LE/(s.H+s.LE))
pl.plot(sr.t/3600., sr.LE/(sr.H+sr.LE))
pl.plot(su.t/3600., su.LE/(su.H+su.LE))
pl.xlabel('time (h)')
pl.ylabel('EF (-)')

pl.subplot(224)
pl.plot(s.t /3600., s.H +s.LE,  label='total')
pl.plot(sr.t/3600., sr.H+sr.LE, label='rural')
pl.plot(su.t/3600., su.H+su.LE, label='urban')
pl.xlabel('time (h)')
pl.ylabel('H+LE (W/m2)')
