# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 21:30:01 2018

@author: rory
"""

from obspy.core import read
from obspy import UTCDateTime
import sys
import numpy as np
from obspy.signal.trigger import recursive_sta_lta, trigger_onset, classic_sta_lta
import matplotlib.pyplot as plt
from conv2d import conv2d
start = UTCDateTime(sys.argv[1])
tr = read(sys.argv[2])
dt=UTCDateTime(sys.argv[3])
tr.trim(starttime=start, endtime=dt)
df = tr[0].stats.sampling_rate


# Filtering with a lowpass on a copy of the original Trace
tr_filt = tr.copy()
tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)
trace=tr_filt
trace=conv2d(trace)

df = trace[0].stats.sampling_rate
npts = trace[0].stats.npts
cft = recursive_sta_lta(trace[0].data*1000, int(0.5 * df), int(5 * df))
bft=classic_sta_lta(trace[0].data*1000,int(0.5*df),int(5*df))
trigger_onset(cft, 2.5, 1.5, 300, max_len_delete=False)

t = np.arange(npts, dtype=np.float32) / df
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(t, trace[0].data*1000, 'k')
ax1.set_ylabel("Displacement (mm)",fontsize=16)
ax2 = fig.add_subplot(212, sharex=ax1)
ax2.plot(t, cft, 'k')
on_off = np.array(trigger_onset(cft, 2.5, 1.5))
i, j = ax1.get_ylim()
try:
    ax1.vlines(on_off[:, 0] / df, i, j, color='r', lw=2,
               label="Trigger On")
    ax1.vlines(on_off[:, 1] / df, i, j, color='b', lw=2,
               label="Trigger Off")
    ax1.legend()
except IndexError:
    pass
ax2.axhline(2.5, color='red', lw=1, ls='--')
ax2.axhline(1.5, color='blue', lw=1, ls='--')
ax2.set_xlabel("Time (s)",fontsize=16)
ax2.set_ylabel("STA / LTA",fontsize=16)
fig.suptitle('Recursive STA/LTA Trigger',fontsize=18)
fig.canvas.draw()
plt.show()



#plot_trigger(tr_filt[0], cft, 2.5, 1.0)
#p_pick, phase_info = pk_baer(tr_filt.data, df,20, 20, 60, 12.0, 100, 100)
#print(p_pick)
# Now let's plot the raw and filtered data...
"""
t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
plt.subplot(211)
plt.plot(t, tr.data, 'k')
plt.ylabel('Raw Data')
plt.subplot(212)
plt.plot(t, tr_filt.data, 'k')
plt.ylabel('Highpassed Data')
plt.xlabel('Time [s]')
plt.suptitle((tr[0].stats.station,tr[0].stats.channel,tr[0].stats.starttime))
plt.show()

#from obspy.signal.trigger import classic_sta_lta
#ft = classic_sta_lta(tr_filt.data, int(2 * df), int(10 * df))
#plot_trigger(tr_filt, cft, 4.5, 1.0)

"""