# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 13:43:03 2019

@author: rory
"""

from obspy.core import read
from obspy import UTCDateTime
import obspy
import sys
import numpy as np
import matplotlib.pyplot as plt
t1 = obspy.UTCDateTime("2018-10-29-T11:30:30")
st = read("LV.L004..HHZ.D.2018.302", starttime=t1, endtime=t1+60*20)
#st2=read("UR.AQ03..HHE.D.2018.302")
#st3=read("UR.AQ03..HHN.D.2018.302")

st.merge(method=0, fill_value=0, interpolation_samples=0)
#st2.merge(method=0, fill_value=0, interpolation_samples=0)
#st3.merge(method=0, fill_value=0, interpolation_samples=0)


data=st[0].data

t = np.arange(0, st[0].stats.npts / st[0].stats.sampling_rate, st[0].stats.delta)
#f = Fs*(0:(L/2))/L;

freq=np.multiply(100,np.arange(0,(st[0].stats.npts)/2,st[0].stats.delta*st[0].stats.sampling_rate))/len(data)
datafft=np.abs(np.fft.fft(data))
fre12=freq[:-1]
#P1(2:end-1) = 2*P1(2:end-1);
a=len(datafft)
datafft2=datafft[0:60000]
plt.figure()
plt.plot(fre12,datafft2)
plt.yscale("log")
plt.ylabel('Amplitude (counts) (log)')
plt.xlabel('Frequency (Hz)')
plt.title('One-Sided Fourier Amplitude Spectrum')
#plt.plot(data)
zchan=st.copy()
#zchan2=st2.copy()
#zchan3=st3.copy()
print('red=z, blue=e, yellow=n')
zchan.filt=zchan.filter('bandpass',freqmin=5, freqmax=45, corners=6)
#zchan2.filt=zchan2.filter('bandpass',freqmin=5, freqmax=45, corners=6)
#zchan3.filt=zchan3.filter('bandpass',freqmin=5, freqmax=45, corners=6)
zchan.filt.plot(color='red', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +2)
#zchan2.filt.plot(color='blue', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +2)
#zchan3.filt.plot(color='yellow', number_of_ticks=1, tick_rotations=1, tick_format='%I:%M %p', starttime=t1-1, endtime=t1 +2)
