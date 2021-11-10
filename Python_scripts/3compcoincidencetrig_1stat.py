# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 14:04:25 2019

@author: rory
"""


from obspy.core import read
from obspy import UTCDateTime

dt = UTCDateTime("2018-10-29T11:02:00")
st=read()
files = ["LV.L004..HHZ.D.2018.302","LV.L004..HHE.D.2018.302","LV.L004..HHN.D.2018.302"]

for filename in files:
    st += read(filename,starttime=dt, endtime=dt+60*60*1)
    
   # st.merge()
#st = read("LV.L004..HHZ.D.2018.297",  starttime=dt, endtime=dt+60*60*4)[0]
#df = st.stats.sampling_rate

# Filtering with a lowpass on a copy of the original Trace
tr_filt = st.copy()
tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)


from obspy.signal.trigger import coincidence_trigger
st2 = tr_filt.copy()
trig = coincidence_trigger("recstalta", 3.5, 1, st2, 3, sta=0.5, lta=10)
from pprint import pprint
pprint(trig)