# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 21:30:01 2018

@author: rory
"""

from obspy.core import read
from obspy import UTCDateTime
from obspy.signal.trigger import plot_trigger, pk_baer

dt = UTCDateTime("2018-10-27T09:27:19.5")
tr = read("LV.L009..HHE.D.2018.300", starttime=dt-1, endtime=dt+2)[0]
df = tr.stats
df

# Filtering with a lowpass on a copy of the original Trace
tr_filt = tr.copy()
tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)



#p_pick, phase_info = pk_baer(tr_filt.data, df,20, 20, 60, 12.0, 100, 100)
#print(p_pick)
# Now let's plot the raw and filtered data...
t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
plt.subplot(211)
plt.plot(t, tr.data, 'k')
plt.ylabel('Raw Data')
plt.subplot(212)
plt.plot(t, tr_filt.data, 'k')
plt.ylabel('Highpassed Data')
plt.xlabel('Time [s]')
plt.suptitle(tr.stats.starttime)
plt.show()

#from obspy.signal.trigger import classic_sta_lta
#ft = classic_sta_lta(tr_filt.data, int(2 * df), int(10 * df))
#plot_trigger(tr_filt, cft, 4.5, 1.0)

