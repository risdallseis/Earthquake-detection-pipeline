# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 14:23:36 2019

@author: rory
"""

from obspy.core import read
from obspy import UTCDateTime
import copy
from tabulate import tabulate
import numpy as np
from obspy.signal.trigger import coincidence_trigger
from pprint import pprint
from obspy.signal.cross_correlation import correlate as xcorr
from obspy.signal.cross_correlation import xcorr_max 
from obspy.signal.filter import envelope

#Trace data- looking for earthquakes
starttr = UTCDateTime("2018-10-29T11:02:00")
end = starttr + 60*60*1
st=read()
#files to sta/lta
files = ["LV.L004..HHZ.D.2018.302","LV.L004..HHE.D.2018.302","LV.L004..HHN.D.2018.302","LV.L002..HHZ.D.2018.302","LV.L002..HHE.D.2018.302","LV.L002..HHN.D.2018.302"]
#files to xcorr
files2xcor=["LV.L002..HHZ.D.2018.302","LV.L002..HHE.D.2018.302","LV.L002..HHN.D.2018.302"]
filesxcor = ["LV.L004..HHZ.D.2018.302","LV.L004..HHE.D.2018.302","LV.L004..HHN.D.2018.302"]
#Zcomps
stationzcomps=["LV.L002..HHZ.D.2018.302","LV.L004..HHZ.D.2018.302"]


# termplate time/length
stt=UTCDateTime("2018-10-29T11:30:40")
dtt=2
xcorr_threshold=0.35

tr=read()

for filename in files:
	st += read(filename,starttime=starttr, endtime=starttr+60*60*1)
	#print(st)
	# Filtering with a lowpass on a copy of the original Trace
	tr_filt = st.copy()
	tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)
	st2 = tr_filt.copy()
	trig = coincidence_trigger("recstalta", 2.5, 1, st2, 3, sta=0.5, lta=2.0)

possquakes=[(d['coincidence_sum'],d['time'],d['trace_ids']) for d in trig]
output=[item for item in possquakes if item[0] == 3]
maybquakes=[output]
possqtimes=[d[1] for d in output]

xcorr_earthquakes=list([])

output2=[item for item in possquakes if item[0] >3]
stalta_earthquakes=[output2]


#stalta_eq_table=tabulate([[(d[stalta_earthquakes[0][0], stalta_earthquakes[0][1], stalta_earthquakes[0][2])], (d[stalta_earthquakes[1][1], stalta_earthquakes[1][1], stalta_earthquakes[1][2]])], headers=['Trigger Sum', 'Time', 'Stations'] for d in stalta_earthquakes])
print (stalta_earthquakes)#need to change the print to make it clearer

counter=0

#create command for if all 3 components trigger 
for d in possqtimes:
	counter=0
	for filen, filem in zip(files2xcor, filesxcor):
		st += read(filen,starttime=d-2, endtime=d+10)
		tr += read(filen,starttime=stt,endtime=stt+dtt)
		tr.trim(starttime=stt,endtime=stt+dtt)
		tr2 = copy.deepcopy(st)
		if len(tr2) > 0:
#tr.detend(type='simple')
			tr.taper(0.01, type='hann', max_length=None, side='both')
			tr.filter('bandpass',freqmin=5, freqmax=45, corners=6)
			tr_en = envelope(tr[0].data)
			tr2.merge(method=0, fill_value=0, interpolation_samples=0)
			tr2.trim(starttime=d,endtime=d+dtt)
			tr2.detrend(type='simple')
			tr2.taper(0.01, type='hann', max_length=None, side='both')
			tr2.filter('bandpass',freqmin=5, freqmax=45, corners=6)
			tr2_en = envelope(tr2[0].data)
	#		cc = xcorr(tr2[0].data/np.max(np.abs(tr2[0].data)), tr[0].data/np.max(np.abs(tr[0].data)), tr2[0].stats.npts, demean=True, normalize=True, domain='freq')
			cc = xcorr(tr2[0].data, tr[0].data, tr2[0].stats.npts, demean=True, normalize=True, domain='freq')
			shift, value = xcorr_max(cc, abs_max=True)
			cc_en = xcorr(tr2_en, tr_en, tr2[0].stats.npts, demean=True, normalize=True, domain='freq')
			shift_en, value_en = xcorr_max(cc_en, abs_max=True)
			times = tr2[0].times('utcdatetime')
			tbest=int(np.min([tr2[0].stats.npts-1,int((tr2[0].stats.npts/2.))+shift]))
			pick = times[tbest]
			tbest_en=int(np.min([tr2[0].stats.npts-1,int((tr2[0].stats.npts/2.))+shift_en]))
			pick_en = times[tbest_en]
			if np.abs(value) >= xcorr_threshold or np.abs(value_en) >= np.min([xcorr_threshold*1.5, 1]):
				print ('station:%s %s %s to %s - Best Xcorr: %s %.2f %.2f %s %.2f %.2f' % (filen,filem,d, d+dtt, pick, value, xcorr_threshold, pick_en, value_en, xcorr_threshold*1.5))
				counter=counter+1
				nocomp=counter
				print('Number of components triggered: %s' %(nocomp))
				if nocomp==3:
					print('EARTHQUAKE DETECTED at: %s'%(pick))
					xcorr_earthquakes.append(d)
                    
            
#Creating a list of all the trigger times from both methods
Tquakes=[d[1] for d in stalta_earthquakes[0]]
for d in xcorr_earthquakes:
    Tquakes.append(d)

fracking_loc=np.array([[-2.95],[53.79]])
