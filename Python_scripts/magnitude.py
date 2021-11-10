# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 17:04:25 2019

@author: rory
"""
def magnitude(data,start,dt,hypo_x):
	import copy
	import numpy as np
	from obspy.core import read
	from obspy.core import UTCDateTime
	from obspy.signal.invsim import paz_to_freq_resp


	tr = read(data)
	stt=(UTCDateTime(start))
	dt=float(dt)
	tr.trim(starttime=stt,endtime=stt+dt);
	tr.detrend(type='simple')
	tr.taper(0.01, type='hann', max_length=None, side='both')
	tr.filter('bandpass',freqmin=5, freqmax=40, corners=6)
	tr1 = copy.deepcopy(tr)
	pre_filt = (0.3, 0.5, 45.0, 50.0)
	seedresp = {'filename': 'RESP.LV.L00X.00.HHX',  # RESP filename
            # when using Trace/Stream.simulate() the "date" parameter can
            # also be omitted, and the starttime of the trace is then used.
            #'date': date,
            # Units to return response in ('DIS', 'VEL' or ACC)
            # Use DIS for WA calculation
            'units': 'VEL'
            }
	tr.simulate(paz_remove=None, pre_filt=pre_filt, remove_sensitivity=True, seedresp=seedresp)


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

# Check WA response
#h, f = paz_to_freq_resp(poles, zeros, gain*sens, 0.005, 16384, freq=True)
#plt.figure()
#plt.loglog(f,abs(h))
#plt.xlabel('Frequency (Hz)')
#plt.ylabel('Sensitivity')

# Remove the Wood-Anderson Response
	tr2 = copy.deepcopy(tr1)
	tr2.simulate(paz_remove=None, paz_simulate=wa, simulate_sensitivity=True)
	WA_peak = np.max(np.abs(tr2[0].data))*1e9
	R = float(hypo_x) # Hypocentral Distance (km)
	BGS_mag=np.log10(WA_peak) + 1.11*np.log10(R) + 0.00189*R - 1.16*np.exp(-0.2*R) - 2.09, np.log10(WA_peak) + 1.11*np.log10(2.*R) + 0.00189*2.*R - 1.16*np.exp(-0.2*2.*R) - 2.09, np.log10(WA_peak) + 1.11*np.log10(5.*R) + 0.00189*5.*R - 1.16*np.exp(-0.2*5.*R) - 2.09
	o_sarg_mag=np.log10(WA_peak) + 1.06*np.log10(R) + 0.00182*R - 1.98, np.log10(WA_peak) + 1.06*np.log10(2.*R) + 0.00182*2.*R - 1.98, np.log10(WA_peak) + 1.06*np.log10(5.*R) + 0.00182*5.*R - 1.98
	return BGS_mag, o_sarg_mag