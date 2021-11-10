# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 12:38:49 2019

@author: rory
"""

#add more stations to STALTA

from obspy.core import read
from obspy import UTCDateTime
import copy
import random

import numpy as np
from obspy.signal.trigger import coincidence_trigger, pk_baer
from pprint import pprint
from obspy.signal.cross_correlation import correlate as xcorr
from obspy.signal.cross_correlation import xcorr_max 
from obspy.signal.filter import envelope
from obspy.signal.invsim import paz_to_freq_resp
from numpy.linalg import inv
from magnitude import magnitude
import sys
#Trace data- looking for earthquakes
starttr = UTCDateTime(sys.argv[1])
end = starttr + (60*60*12)
st=read()

#files to sta/lta (1)
files = [sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10],sys.argv[11],sys.argv[12],sys.argv[13]]
#less noisy stations to sta/lta
files1=[sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[11],sys.argv[12],sys.argv[13]]
#files to xcorr
xcorrtemplates=[sys.argv[15],sys.argv[16]]
files2xcorr=[sys.argv[6],sys.argv[7]]


#station locations MUST BE IN THE SAME ORDER AS Z COMPS
s1=np.array([[-2.84,53.84]])#L2
s2=np.array([[-2.94,53.76]])#L4
s3=np.array([[-2.88,53.69]])#L6
s4=np.array([[-3.04,53.78]])#L9
velocity=5#Seismic wave velocity
# termplate time/length
stt=UTCDateTime("2018-10-27T10:55:26")
dtt=4
xcorr_threshold=0.35
tr=read()
tr1=read()
#for filename in files:
#	st += read(filename,starttime=starttr, endtime=end)
#	#print(st)
#	# Filtering with a lowpass on a copy of the original Trace
#	tr_filt = st.copy()
#	tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)
#	st2 = tr_filt.copy()
#	trig = coincidence_trigger("recstalta", 1.9, 1, st2, 12, sta=0.5, lta=2.5)#4 stations must trigger in order to be recorded
    

for filename in files1:
	st+=read(filename,starttime=starttr, endtime=end)
	tr_filt=st.copy()
	tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)
	st3=tr_filt.copy()
	trig2=coincidence_trigger("recstalta",1.3,1,st3,6,sta=0.5,lta=2.5)
    #THIKABOUT THIS# def qUAKES ^^ TO XCORR
    
    

#stalta_earthquakes_info=[(d['coincidence_sum'],d['time'],d['trace_ids']) for d in trig]
#print(stalta_earthquakes_info)
#stalta_earthquakes=[d['time'] for d in trig]

possquakes=[d['time'] for d in trig2]#INDEX THE TIME SLO I.E [0][1]
#output=[item for item in possquakes if item[0] >=6]
#maybquakes=[output]
#possqtimes=[d[1] for d in output]


xcorr_earthquakes=[]
#print (stalta_earthquakes)#need to change the print to make it clearer
print('detections to be Xcorr')
print(possquakes)
counter=0
xcorr_val=[]
xcorr_en_val=[]
#create command for if all 3 components trigger 
for d in possquakes:
	counter=0
	total_coeff=[]
	total_coeff_en=[]
	st =read(sys.argv[5],starttime=d-1, endtime=d+4)#xcorr from same station
	tr =read(sys.argv[14],starttime=stt,endtime=stt+dtt)
	tr2 = copy.deepcopy(st)
	if len(tr2) > 0:
		tr.detrend(type='simple')
		tr.taper(0.01, type='hann', max_length=None, side='both')
		tr.filter('bandpass',freqmin=5, freqmax=45, corners=6)
		tr_en = envelope(tr[0].data)
		tr2.merge(method=0, fill_value=0, interpolation_samples=0)
		tr2.detrend(type='simple')
		tr2.taper(0.01, type='hann', max_length=None, side='both')
		tr2.filter('bandpass',freqmin=5, freqmax=45, corners=6)
		tr2_en = envelope(tr2[0].data)
		cc = xcorr(tr2[0].data, tr[0].data, tr2[0].stats.npts, demean=True, normalize=True, domain='freq')
		shift, value = xcorr_max(cc, abs_max=True)
		cc_en = xcorr(tr2_en, tr_en, tr2[0].stats.npts, demean=True, normalize=True, domain='freq')
		shift_en, value_en = xcorr_max(cc_en, abs_max=True)
		times = tr2[0].times('utcdatetime')
		tbest=int(np.min([tr2[0].stats.npts-1,int((tr2[0].stats.npts/2.))+shift]))
		pick = times[tbest]
		tbest_en=int(np.min([tr2[0].stats.npts-1,int((tr2[0].stats.npts/2.))+shift_en]))
		pick_en = times[tbest_en]
		print ('station: %s temp:%s %s %s to %s - Best Xcorr: %s %.2f %.2f %s %.2f %.2f' % (st[0].stats.station,st[0].stats.channel,tr[0].stats.channel,d, d+dtt, pick, value, xcorr_threshold, pick_en, value_en, xcorr_threshold*1.5))
		if np.abs(value) >= xcorr_threshold or np.abs(value_en) >= np.min([xcorr_threshold*1.5, 1]):
			counter=counter+1
			nocomp=counter
			total_coeff.append(np.abs(value))
			total_coeff_en.append(np.abs(value_en))
			print('Number of components triggered: %s' %(nocomp))
			for filem,filen in zip(xcorrtemplates,files2xcorr) :
				st = read(filen,starttime=d-1, endtime=d+4)#xcorr from same station
				tr1 = read(filem,starttime=stt,endtime=stt+dtt)
				tr2 = copy.deepcopy(st)
				if len(tr2) > 0:
					tr1.detrend(type='simple')
					tr1.taper(0.01, type='hann', max_length=None, side='both')
					tr1.filter('bandpass',freqmin=5, freqmax=45, corners=6)
					tr1_en = envelope(tr[0].data)
					tr2.merge(method=0, fill_value=0, interpolation_samples=0)
					tr2.detrend(type='simple')
					tr2.taper(0.01, type='hann', max_length=None, side='both')
					tr2.filter('bandpass',freqmin=5, freqmax=45, corners=6)
					tr2_en = envelope(tr2[0].data)
					cc = xcorr(tr2[0].data, tr1[0].data, tr2[0].stats.npts, demean=True, normalize=True, domain='freq')
					shift, value = xcorr_max(cc, abs_max=True)
					cc_en = xcorr(tr2_en, tr1_en, tr2[0].stats.npts, demean=True, normalize=True, domain='freq')
					shift_en, value_en = xcorr_max(cc_en, abs_max=True)
					times = tr2[0].times('utcdatetime')
					tbest=int(np.min([tr2[0].stats.npts-1,int((tr2[0].stats.npts/2.))+shift]))
					pick = times[tbest]
					tbest_en=int(np.min([tr2[0].stats.npts-1,int((tr2[0].stats.npts/2.))+shift_en]))
					pick_en = times[tbest_en]
					print ('station: %s temp:%s %s %s to %s - Best Xcorr: %s %.2f %.2f %s %.2f %.2f' % (st[0].stats.station,tr1[0].stats.channel,st[0].stats.channel,d, d+dtt, pick, value, xcorr_threshold, pick_en, value_en, xcorr_threshold*1.5))
					if np.abs(value) >= xcorr_threshold or np.abs(value_en) >= np.min([xcorr_threshold*1.5, 1]):
						counter=counter+1
						nocomp=counter
						total_coeff.append(np.abs(value))
						total_coeff_en.append(np.abs(value_en))
						print('Number of components triggered: %s' %(nocomp))
						if nocomp==3:
							xcorr_val.append((np.abs(total_coeff)))
							xcorr_en_val.append((np.abs(total_coeff_en)))
							print('EARTHQUAKE DETECTED at: %s'%(d))
							xcorr_earthquakes.append(d)
xcorr_earthquakes_singles=[]
xcorr_results=xcorr_earthquakes + xcorr_val +  xcorr_en_val
#removing duplicate times

for i in xcorr_earthquakes:
    if i not in xcorr_earthquakes_singles:
        xcorr_earthquakes_singles.append(i)
Tquakes=[]
#Creating a list of all the trigger times from both methods
#for i in stalta_earthquakes:
#    Tquakes.append(i)

for d in xcorr_earthquakes_singles:
    Tquakes.append(d)

fracking_loc=np.array([[-2.95,53.79]])#fracking site coords
s_locs=np.array([[s1],[s2],[s3],[s4]])#station locations in lat/lon
s_lon=np.array([[s1[0][0]],[s2[0][0]],[s3[0][0]],[s4[0][0]]])#station longitudes
s_lat=np.array([[s1[0][1]],[s2[0][1]],[s3[0][1]],[s4[0][1]]])#station latitudes
#lat/lon to cartesian with fracking station at x,y = (0,0)
X=(s_lon-fracking_loc[0][0])*40000*np.cos((s_lat+fracking_loc[0][1])*np.pi/360)/360
Y=(fracking_loc[0][1]-s_lat)*40000/360
s_z=np.array([[0],[0],[0],[0]])
s_clocs=np.column_stack((X,Y,s_z))#station locations (X,Y,Z)
hypo=np.array([[0],[0],[-2.5]])
station_x_hypo=(np.square(s_clocs[:,0]))+ (np.square((s_clocs[:,1])) + (np.square(s_clocs[:,2])))
station2hypo=np.power(station_x_hypo,0.5)


#Zcomps
ncomp4=sys.argv[7]
ncomp9=sys.argv[13]
zcomp9=sys.argv[11]
zcomp4=sys.argv[5]

eq_locs=[]
eq_mags_bgs=[]
eq_mags_sarg=[]
for i in Tquakes:
	BGS_mag, o_sarg_mag=magnitude(ncomp4,i,2,station2hypo[1])#using station9 to calc mag
	eq_mags_bgs.append(BGS_mag)
	eq_mags_sarg.append(o_sarg_mag)
	p_arrivals=[]
	zt = read(zcomp4,starttime=(UTCDateTime(i))-2, endtime=(UTCDateTime(i))+3) #TEST P-WAVE PICKS WITHOUT HAVING TO ADJUST TEMPLATE START DURATION FOR EACH STATION
	zt1=read(zcomp9,starttime=(UTCDateTime(i))-2, endtime=(UTCDateTime(i))+3)
	tr_filt = zt.copy()[0]
	tr_filt2 = zt1.copy()[0]
	df=tr_filt.stats.sampling_rate
	tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)
	st2 = tr_filt.copy()
	tr_filt2.filter('bandpass', freqmin=5, freqmax=45, corners=6)
	st3 = tr_filt2.copy()
	p_pick4, phase_info = pk_baer(st2.data, df,20, 20, 60, 12.0, 100, 100)
	if p_pick4<2:
		p_pick4=random.uniform(station2hypo[1]/0.055,station2hypo[1]/0.045)#if no p-wave if found then a random pick is generated
	p_arrivals.append(p_pick4/100)
	p_pick9, phase_info = pk_baer(st3.data, df,20, 20, 60, 12.0, 100, 100)
	if p_pick9<2:
		p_pick9=random.uniform(station2hypo[3]/0.055,station2hypo[3]/0.045)
	p_arrivals.append(p_pick9/100)
	D1=(p_pick4/100)+((station2hypo[0]-station2hypo[1])/velocity)
	D3=(p_pick4/100)+((station2hypo[2]-station2hypo[1])/velocity)
	p_arrivals.insert(0,D1)
	p_arrivals.insert(2,D3)
	Do=np.reshape(np.asarray(p_arrivals),(4,1))
#starting model
	Mi=np.array([[0],[0],[-2.5],[-1],[0]])
#X,Y,Z distances from starting model to stations
	R=np.subtract(s_clocs,np.transpose(Mi[0:3]))
#straight line ray paths from hypocentre
	Rh=(np.square(R[:,0]))+ (np.square((R[:,1])) + (np.square(R[:,2])))
	Rh=np.power(Rh,0.5)
	Rh=np.reshape(Rh,(4,1))
#Ray travel times calculated using Rh
	Tt=Rh/velocity
#Predicted P-wave arrival times
	Dp=np.transpose(Tt + Mi[3,:])
	Dp=np.reshape(Dp,(4,1))
#Inverting the residual of observed and calculated arrival times for the change in model (D'=GM')
	iterations=1
	G=np.zeros([len(s_clocs),4])
	misfit_array=[]
#Earthquake location iteration while loop
	while iterations<8:
    #Arrival time residual (data D')
		Dres=Do-Dp
    #Total data misfit
		misfit=np.sum([np.square(Dres)])
		misfit_array.append(misfit)
    #Connecting matrix G, partial derivatives of time with respect to (x,y,z,t)
		G1=-(np.reshape(R[:,0],(4,1))/velocity)*(np.power(Rh,-0.5))
		G2=-(np.reshape(R[:,1],(4,1))/velocity)*(np.power(Rh,-0.5))
		G3=-(np.reshape(R[:,2],(4,1))/velocity)*(np.power(Rh,-0.5))
		G4=np.ones([(len(s_clocs)),1])
		G=np.column_stack((G1,G2,G3,G4))
    #Lagrange multiplier to reduce weighting on Z
		h=0
		F=np.array([[0],[0],[1],[0]])
		H1=np.column_stack([np.dot(np.transpose(G),G), F])
		H2= np.append(np.transpose(F), F[0])
		H2=H2.astype(np.float64)
		H=np.vstack((H1,H2))
		Ginv=np.linalg.inv(H)
    #Previous model
		M_prev=Mi
    #calculating the change in model
    #constructing model change matrix
		M1=np.transpose(G)*Dres
		M2=np.append(M1,F[0])
		delta_model1=np.dot(np.transpose(G),Dres)
		delta_model1=np.append(delta_model1,F[0])
		delta_model=np.dot(Ginv,delta_model1)
    #new model
		Mi=M_prev+np.reshape(delta_model,(5,1))
    #calculating distances from station to the new model
		R=np.subtract(s_clocs,np.transpose(Mi[0:3]))
		Rh=(np.square(R[:,0]))+ (np.square((R[:,1])) + (np.square(R[:,2])))
		Rh=np.power(Rh,0.5)
		Rh=np.reshape(Rh,(4,1))
    #calcualting new travel times
		Tt=Rh/velocity
    #new calculated pick times
		Dp=np.transpose(Tt + Mi[3,:])
		Dp=np.reshape(Dp,(4,1))
		iterations=iterations+1
		if iterations>7:
			eq_locs.append(Mi)

#Cartesian to lon/lat
x_dist=[d[0] for d in eq_locs]
y_dist=[d[1] for d in eq_locs]
eq_times=[d[3] for d in eq_locs]
s_mlon=np.mean(s_lon)
s_mlat=np.mean(s_lat)
eq_lat=(np.add(y_dist,(1000/9)*(s_mlat)))/(1000/9)
eq_lon=s_mlon-(((np.multiply(360,x_dist))/40000)/np.cos((s_mlat+eq_lat)*np.pi/360))
eq_locs_lonlat=np.column_stack((eq_lon,eq_lat,eq_times))

#OUTPUTTING RESULTS
#with open("stalta_results.txt", "w") as output:
#        output.write("\n")
#        output.write(str(stalta_earthquakes))
        
with open("total_earthquakes.txt", "w") as output:
        output.write("\n")
        output.write(str(Tquakes))
        
with open("eq_locations.txt", "w") as output:
        output.write("\n")
        output.write(str(eq_locs_lonlat))
        
np.savetxt('eq_mags_bgs.out',eq_mags_bgs)




with open("xcorr_results.txt", "w") as output:
        output.write("\n")
        output.write(str(xcorr_results))
        
print(xcorr_earthquakes_singles)