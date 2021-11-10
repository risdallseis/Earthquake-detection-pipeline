# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 21:30:01 2018

@author: rory
"""

from obspy.core import read
from obspy import UTCDateTime
import sys
import numpy as np
import matplotlib.pyplot as plt



start = UTCDateTime(sys.argv[1])
#dt=(sys.argv[3])
tr = read((sys.argv[2]))
tr.trim(starttime=start, endtime=start+10)
tr.merge(method=0, fill_value=0, interpolation_samples=0)



df = tr[0].stats


# Filtering with a lowpass on a copy of the original Trace
tr_filt = tr.copy()
tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)
l1='raw data',
l2='filtered data',


#p_pick, phase_info = pk_baer(tr_filt.data, df,20, 20, 60, 12.0, 100, 100)
#print(p_pick)
# Now let's plot the raw and filtered data...
t = np.arange(0, tr[0].stats.npts / tr[0].stats.sampling_rate, tr[0].stats.delta)
plt.subplot(211)
plt.plot(t, tr[0].data, 'k')
plt.tick_params(labelsize=12)
plt.ylabel('sensor counts',fontsize=16)
plt.legend(l1,fontsize=12)
plt.subplot(212)
plt.plot(t, tr_filt[0].data, 'k')
plt.tick_params(labelsize=12)
plt.ylabel('sensor counts',fontsize=16)
plt.xlabel('Time [s]',fontsize=16)
plt.legend(l2,fontsize=12)
plt.suptitle('Unfiltered and Filtered  Waveform',fontsize=18)
plt.show()

#from obspy.signal.trigger import classic_sta_lta
#ft = classic_sta_lta(tr_filt.data, int(2 * df), int(10 * df))
#plot_trigger(tr_filt, cft, 4.5, 1.0)

