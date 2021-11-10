# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 10:47:05 2019

@author: rory
"""

import sys
import copy
import numpy as np
import matplotlib.pyplot as plt

from obspy.core import read
from obspy.signal.invsim import paz_to_freq_resp
from obspy.core import UTCDateTime

tr = read(sys.argv[1])
tr.merge(method=0, fill_value=0, interpolation_samples=0)
stt=UTCDateTime(sys.argv[2])
dt=float(sys.argv[3])
tr.trim(starttime=stt,endtime=stt+dt)
tr.detrend(type='simple')
tr.taper(0.01, type='hann', max_length=None, side='both')
#tr.filter('bandpass',freqmin=5, freqmax=40, corners=6)
tr1 = copy.deepcopy(tr)
npts = tr[0].stats.npts
samprate = tr[0].stats.sampling_rate
t = np.arange(0, ((npts-0.5) / samprate), 1 / samprate)
plt.figure()
plt.title((tr[0].stats.station,tr[0].stats.channel,tr[0].stats.starttime))

pre_filt = (0.3, 0.5, 45.0, 50.0)

seedresp = {'filename': 'RESP.LV.L00X.00.HHX',  # RESP filename
            # when using Trace/Stream.simulate() the "date" parameter can
            # also be omitted, and the starttime of the trace is then used.
            #'date': date,
            # Units to return response in ('DIS', 'VEL' or ACC)
            # Use DIS for WA calculation
            'units': 'VEL'
            }


# Remove the T120 instrument response
tr.simulate(paz_remove=None, remove_sensitivity=True, seedresp=seedresp)
PGV = np.max(np.abs(tr[0].data))*1e3
plt.plot(t, tr[0].data*1000., 'b-')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (mm/s)')
plt.legend(['Ground Motion'])
plt.show()
"""
seedresp = {'filename': 'RESP.LV.L00X.00.HHX',  # RESP filename
            # when using Trace/Stream.simulate() the "date" parameter can
            # also be omitted, and the starttime of the trace is then used.
            #'date': date,
            # Units to return response in ('DIS', 'VEL' or ACC)
            # Use DIS for WA calculation
            'units': 'DIS'
            }


# Remove the T120 instrument response
tr1.simulate(paz_remove=None, pre_filt=pre_filt, remove_sensitivity=True, seedresp=seedresp)
PGD = np.max(np.abs(tr1[0].data))*1e3
plt.subplot(412)
plt.plot(t, tr1[0].data*1000., 'r-')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (mm)')
plt.legend(['Ground Data'])

poles =  [-5.49779+5.60886j,-5.49779-5.60886j]
zeros = [0j, 0j]
#gain = 2080.
#based on the M4 example seems like the WA gain is not applied
gain = 1.
sens = 1.0
wa = {'poles': poles,
      'zeros': zeros,
      'gain': gain,
      'sensitivity': sens 
     }

# Remove the Wood-Anderson Response
tr2 = copy.deepcopy(tr1)
tr2.simulate(paz_remove=None, paz_simulate=wa, simulate_sensitivity=True)
plt.subplot(413)
plt.plot(t, tr2[0].data*1000., 'g-')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (mm)')
plt.legend(['Wood-Anderson'])

plt.show()
"""