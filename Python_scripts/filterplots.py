# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 18:11:35 2019

@author: rory
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 14:46:52 2019

@author: rory
"""
from obspy.core import read
from obspy import UTCDateTime
import obspy
import sys
import matplotlib.pyplot as plt
import numpy as np
import conv2d

t1 = obspy.UTCDateTime(sys.argv[1])

st = read((sys.argv[2]),starttime=t1-30, endtime=t1+20)
st2=read((sys.argv[2]),starttime=t1-30, endtime=t1+20)
st3=read((sys.argv[2]),starttime=t1-30, endtime=t1+20)
st.merge(method=0, fill_value=0, interpolation_samples=0)
st2.merge(method=0, fill_value=0, interpolation_samples=0)
st3.merge(method=0, fill_value=0, interpolation_samples=0)
#st.trim(starttime=t1-1, endtime=t1+3)
#st2.trim(starttime=t1-1, endtime=t1+4)
#st3.trim(starttime=t1-1, endtime=t1+4)
zchan=st.copy()
rawchan=st.copy()
zchan2=st2.copy()
zchan3=st3.copy()

print('red=z, blue=e, yellow=n')
zchan.detrend(type='simple')
zchan.taper(0.01, type='hann', max_length=None, side='both')
zchan2.detrend(type='simple')
zchan2.taper(0.01, type='hann', max_length=None, side='both')
zchan3.detrend(type='simple')
zchan3.taper(0.01, type='hann', max_length=None, side='both')
zchan.filt=zchan.filter('highpass',freq=5, corners=6)
zchan2.filt=zchan2.filter('lowpass',freq=45, corners=6)
zchan3.filt=zchan3.filter('bandpass',freqmin=12, freqmax=45, corners=6)
t = np.arange(0, st[0].stats.npts / st[0].stats.sampling_rate, st[0].stats.delta)
plt.figure()

plt.tick_params(labelsize=20)
plt.subplot(411)
plt.plot(t,rawchan[0].data*1000,'k-')
plt.title('Raw Data',fontsize=16)
plt.ylabel('Counts (-)',fontsize=15)
plt.subplot(412)
#plt.suptitle('Highpass, Lowpass and Bandpass Filters')
plt.plot(t,zchan[0].data,'k-')
plt.title('Highpass (5Hz)',fontsize=16)
plt.ylabel('Counts (-)',fontsize=15)

plt.subplot(413)
plt.plot(t,zchan2[0].data,'k-')
plt.title('Lowpass (45Hz)',fontsize=16)
plt.ylabel('Counts (-)',fontsize=15)

plt.subplot(414)
plt.plot(t,zchan3[0].data,'k-')
plt.title('Bandpass (5-45Hz)',fontsize=16)
plt.ylabel('Counts (-)',fontsize=15)
plt.xlabel('Time (s)',fontsize=15)
plt.tight_layout()
plt.show()
#zchan.filt.plot(color='red', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +2)
#zchan2.filt.plot(color='blue', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +2)
#zchan3.filt.plot(color='yellow', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +2)
"""
zchan.plot(color='red', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +60*60)
zchan2.plot(color='blue', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +60*60)
zchan3.plot(color='yellow', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +60*60)
"""


