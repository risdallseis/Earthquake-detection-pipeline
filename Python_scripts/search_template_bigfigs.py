# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 16:14:47 2018

@author: Rory
"""

"""
Created on Thu Nov  8 19:26:03 2018

@author: Rory
"""

import obspy.core
import copy
import numpy as np
import matplotlib.pyplot as plt
from obspy import core
from obspy import read
from obspy.signal.invsim import paz_to_freq_resp
from obspy import UTCDateTime
from obspy.signal.cross_correlation import correlate as xcorr
from obspy.signal.cross_correlation import xcorr_max 
from obspy.signal.filter import envelope

# python search_template.py data/20181029/LV.L004..HHZ.D.2018.302 2018-10-29T11:30:38 10 data/20181029/LV.L004..HHZ.D.2018.302 2018-10-29T09:30:00 2018-10-29T11:59:59 0.5

# python search_template.py TEMPLATE_MS TEMPLATE_TIME TEMPLATE_LEN TARGET_MS TARGET_START TARGET_END XCORR_LOWER_LIMIT

# read the temlpate dat
tr = read("LV.L004..HHZ.D.2018.302")

# read the new data
trF=read()
stations = ["LV.L004..HHZ.D.2018.302"]
for filename in stations:
	trF += read(filename)
	print(tr)
	print(trF)
# fix gaps (merge to single trace)
#tr.merge()

# termplate time/length
	stt=UTCDateTime("2018-10-29T11:30:40")
	dt=4
	xcorr_threshold=0.4
# cut out template
	tr.trim(starttime=stt,endtime=stt+dt)
	tr.detrend(type='simple')
	tr.taper(0.01, type='hann', max_length=None, side='both')
	tr.filter('bandpass',freqmin=5, freqmax=45, corners=6)
	tr_en = envelope(tr[0].data)

	start=UTCDateTime("2018-10-29T11:30:40")
	end=start+60*20

#print ('Search Start: %s; Search End: %s' % (start, end))
	for offset in np.linspace(0, end-start-dt, 2.*int(1+(end-start-dt)/dt)):
		tr2 = copy.deepcopy(trF)
		if len(tr2) > 0:
			tr2.merge(method=0, fill_value=0, interpolation_samples=0)
			tr2.trim(starttime=start+offset,endtime=start+dt+offset)
			tr2.detrend(type='simple')
			tr2.taper(0.01, type='hann', max_length=None, side='both')
			tr2.filter('bandpass',freqmin=5, freqmax=45, corners=6)
			tr2_en = envelope(tr2[0].data)
		#	cc = xcorr(tr2[0].data/np.max(np.abs(tr2[0].data)), tr[0].data/np.max(np.abs(tr[0].data)), tr2[0].stats.npts, demean=True, normalize=True, domain='freq')
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
				print ('%s %s to %s - Best Xcorr: %s %.2f %.2f %s %.2f %.2f' % (filename, start+offset, start+dt+offset, pick, value, xcorr_threshold, pick_en, value_en, xcorr_threshold*1.5))
				npts = tr[0].stats.npts
				samprate = tr[0].stats.sampling_rate
				t = np.arange(0, ((npts-0.5) / samprate), 1 / samprate)
				plt.figure()
				plt.subplots_adjust(top=4,bottom=0,hspace=0)
				plt.subplot(811)
				plt.plot(t-(shift/samprate), tr2[0].data*1., 'g-')
				plt.plot(t-(shift/samprate), tr2_en*1., 'b:')
				plt.ylabel('Sensor Counts (-)')
				plt.xlim([0, dt])
				plt.legend(['Detected'],loc=9, bbox_to_anchor=(1.2, 1.0))
				plt.subplot(813)
				plt.plot(t, tr[0].data*1., 'k-')
				plt.plot(t, tr_en*1., 'b:')
				plt.xlim([0, dt])
				plt.ylabel('Sensor Counts (-)')
				plt.legend(['Template'],loc=9, bbox_to_anchor=(1.2, 1.0))
				plt.subplot(815)
				plt.plot(t, tr[0].data, 'k-')
				plt.plot(t-(shift/samprate), tr2[0].data, 'g-')
				plt.xlim([0, dt])
				plt.ylabel('Sensor Counts (-)')
				plt.legend(['Both (original)'],loc=9, bbox_to_anchor=(1.2, 1.0))
				plt.subplot(817)
				plt.plot(t, tr[0].data*1./np.max(tr[0].data), 'k-')
				plt.plot(t-(shift/samprate), tr2[0].data*1./np.max(tr2[0].data), 'g-')
				plt.xlim([0, dt])
				plt.xlabel('Time (s)')
				plt.ylabel('Sensor Counts (-)')
				plt.legend(['Both (scaled - equal amplitude)'],loc=9, bbox_to_anchor=(1.3, 1.0))
				plt.show()
