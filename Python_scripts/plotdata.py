import sys
import copy
import numpy as np
import matplotlib.pyplot as plt

from obspy.core import read
from obspy.signal.invsim import paz_to_freq_resp
from obspy.core import UTCDateTime

tr = read(sys.argv[1])
stt=UTCDateTime(sys.argv[2])
dt=float(sys.argv[3])
tr.trim(starttime=stt,endtime=stt+dt);
tr.detrend(type='simple')
tr.taper(0.01, type='hann', max_length=None, side='both')
tr.filter('bandpass',freqmin=5, freqmax=40, corners=6)
tr1 = copy.deepcopy(tr)
npts = tr[0].stats.npts
samprate = tr[0].stats.sampling_rate
t = np.arange(0, ((npts-0.5) / samprate), 1 / samprate)
plt.figure()
plt.subplot(411)
plt.plot(t, tr[0].data, 'k-')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (counts)')
plt.legend(['Raw Data'])
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
tr.simulate(paz_remove=None, pre_filt=pre_filt, remove_sensitivity=True, seedresp=seedresp)
PGV = np.max(np.abs(tr[0].data))*1e3
plt.subplot(412)
plt.plot(t, tr[0].data*1000., 'b-')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (mm/s)')
plt.legend(['Ground Data'])

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
plt.subplot(413)
plt.plot(t, tr1[0].data*1000., 'r-')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (mm)')
plt.legend(['Ground Data'])
# Define Wood-Anderson Response
#poles =  [-6.28318+4.71239j, -6.28318-4.71239j]
#zeros = [0j, 0j]
#gain = 2800.
#sens = 1.0
#wa = {'poles': poles,
#      'zeros': zeros,
#      'gain': gain,
#      'sensitivity': sens 
#     }
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
plt.subplot(414)
plt.plot(t, tr2[0].data*1000., 'g-')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (mm)')
plt.legend(['Wood-Anderson'])
# Old WA
#POLES 2
#-6.28318 4.71239
#-6.28318 -4.71239
#ZEROS 2
#CONSTANT 2800.0

# New WA
#POLES 2
#-5.49779 5.60886
#-5.49779 -5.60886
#ZEROS 2
#CONSTANT 2080.0
tr.write(sys.argv[1]+'.wa.mseed')

#WA amplitude in nm
WA_peak = np.max(np.abs(tr2[0].data))*1e9
print ('Peak ground velocity = %.3f mm/s' % (PGV))
print ('Peak ground displacement = %.6f mm' % (PGD))
print ('Peak Wood-Anderson Amplitude (excl. x2080 gain) = %.2f nm' % (WA_peak))
R = float(sys.argv[4]) # Hypocentral Distance (km)
#Hutton & Boore (1987) are used, as shown to be appropriate by Booth (2007).  The equation used by the BGS is:
print ('ML (BGS) [Source at 1x, 2x, 5x %.1f km ] = %.2f %.2f %.2f' % (R, np.log10(WA_peak) + 1.11*np.log10(R) + 0.00189*R - 2.09, np.log10(WA_peak) + 1.11*np.log10(2.*R) + 0.00189*2.*R - 2.09, np.log10(WA_peak) + 1.11*np.log10(5.*R) + 0.00189*5.*R - 2.09))
#Ottemöller & Sargeant (2013) inverted for new values of A and B using UK data and found the relation [without site corrections - note different values if using site terms - see Butcher et al., 2017] (not used by BGS):
print ('ML (Ottemöller & Sargeant, 2013) [Source at 1x, 2x, 5x %.1f km ] = %.2f %.2f %.2f' % (R, np.log10(WA_peak) + 1.06*np.log10(R) + 0.00182*R - 1.98, np.log10(WA_peak) + 1.06*np.log10(2.*R) + 0.00182*2.*R - 1.98, np.log10(WA_peak) + 1.06*np.log10(5.*R) + 0.00182*5.*R - 1.98))
print ('ML (Butcher et al., 2017 - R < 5 km) [Source at 1x, 2x, 5x %.1f km ] =  %.2f %.2f %.2f' % (R, np.log10(WA_peak) + 1.17*np.log10(R) + 0.0514*R - 3.0, np.log10(WA_peak) + 1.17*np.log10(2.*R) + 0.0514*2.*R - 3.0,np.log10(WA_peak) + 1.17*np.log10(5.*R) + 0.0514*5.*R - 3.0))
print ('ML (BGS 2018) [Source at 1x, 2x, 5x %.1f km ] = %.2f %.2f %.2f' % (R, np.log10(WA_peak) + 1.11*np.log10(R) + 0.00189*R - 1.16*np.exp(-0.2*R) - 2.09, np.log10(WA_peak) + 1.11*np.log10(2.*R) + 0.00189*2.*R - 1.16*np.exp(-0.2*2.*R) - 2.09, np.log10(WA_peak) + 1.11*np.log10(5.*R) + 0.00189*5.*R - 1.16*np.exp(-0.2*5.*R) - 2.09))

plt.show()
