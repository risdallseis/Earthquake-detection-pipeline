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
from conv2d import conv2d
t1 = obspy.UTCDateTime(sys.argv[1])

st = read((sys.argv[2]),starttime=t1-1, endtime=t1+3)
st2=read((sys.argv[3]),starttime=t1-1, endtime=t1+3)
st3=read((sys.argv[4]),starttime=t1-1, endtime=t1+3)
st.merge(method=0, fill_value=0, interpolation_samples=0)
st2.merge(method=0, fill_value=0, interpolation_samples=0)
st3.merge(method=0, fill_value=0, interpolation_samples=0)
#st.trim(starttime=t1-1, endtime=t1+3)
#st2.trim(starttime=t1-1, endtime=t1+4)
#st3.trim(starttime=t1-1, endtime=t1+4)
zchan=st.copy()
zchan2=st2.copy()
zchan3=st3.copy()
zchand=conv2d(zchan)
zchan2d=conv2d(zchan2)
zchan3d=conv2d(zchan3)
print('red=z, blue=e, yellow=n')

t = np.arange(0, st[0].stats.npts / st[0].stats.sampling_rate, st[0].stats.delta)
plt.figure()
plt.suptitle('Template: 1.1M earthquake, station 4: NON-BGS detection  -0.7M',fontsize=18)
plt.tick_params(labelsize=16)
plt.subplot(311)
plt.plot(t,zchand[0].data*1000,'b-')
plt.title('Vertical Component',fontsize=16)
plt.ylabel('Displacement (mm)',fontsize=15)
plt.legend([st[0].stats.starttime])
plt.xticks([])

plt.subplot(312)
plt.plot(t,zchan2d[0].data*1000,'b-')
plt.title('East Component',fontsize=16)
plt.ylabel('Displacement (mm)',fontsize=15)
plt.xticks([])

plt.subplot(313)
plt.plot(t,zchan3d[0].data*1000,'b-')
plt.title('North Component',fontsize=16)
plt.ylabel('Displacement (mm)',fontsize=15)
plt.xlabel('Time (s)',fontsize=15)
#plt.tight_layout()
plt.show()
#zchan.filt.plot(color='red', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +2)
#zchan2.filt.plot(color='blue', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +2)
#zchan3.filt.plot(color='yellow', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +2)
"""
zchan.plot(color='red', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +60*60)
zchan2.plot(color='blue', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +60*60)
zchan3.plot(color='yellow', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +60*60)
"""


