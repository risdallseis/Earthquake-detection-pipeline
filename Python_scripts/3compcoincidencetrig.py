# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 14:04:25 2019

@author: rory
"""


from obspy.core import read
from obspy import UTCDateTime
from obspy.signal.trigger import coincidence_trigger
from pprint import pprint
dt = UTCDateTime("2018-10-29T11:02:00")
st=read()
files = ["LV.L004..HHZ.D.2018.302","LV.L004..HHE.D.2018.302","LV.L004..HHN.D.2018.302","LV.L002..HHZ.D.2018.302","LV.L002..HHE.D.2018.302","LV.L002..HHN.D.2018.302"]

for filename in files:
    st += read(filename,starttime=dt, endtime=dt+60*60*1)
    print(st)

#df = st.stats.sampling_rate

# Filtering with a lowpass on a copy of the original Trace
tr_filt = st.copy()
tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)



st2 = tr_filt.copy()
trig = coincidence_trigger("recstalta", 2.50, 1, st2, 3, sta=0.5, lta=2)

pprint(trig)