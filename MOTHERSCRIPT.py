# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 12:38:49 2019

Author: Rory Tisdall
Project: The Detection of Low Amplitude Earthquake Signals

Script name Earthquake Detection, Location and Magnitude MOTHERSCRIPT
"""
#CONSOLE INPUTS: STARTTIME, L2_Z,L2_E,L2_N,L4_Z,L4_E,L4_N,L6_Z,L6_E,L6_N,L9_Z,L9_E,L9_N,TEMPLATE_Z,TEMPLATE_E,TEMPLATE_N
#SCRIPT INPUTS: ENDTIME, CROSS CORRELATION THRESHOLD, STL/LTA PARAMETERS, STATION LOCATIONS, EARTHQUAUKE LOCATION STARTING MODEL
#LOADING MODULES
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
#EARTHQUAKE DETECTION START TIME
starttr = UTCDateTime(sys.argv[1])
#EARTHQUAKE DETECTION END TIME
end = starttr + (60*60*12)
st=read()

#ALL STATIONS AND COMPONENTS IN ORDER Z,E,N 
files = [sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10],sys.argv[11],sys.argv[12],sys.argv[13]]
#CLOSEST 2 STATIONS TO EXCPECTED EARTHQUAKE EPICENTRE
files1=[sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[11],sys.argv[12],sys.argv[13]]
#CROSS-CORRELATION TEMPLATE IN ORDER E,N 
xcorrtemplates=[sys.argv[15],sys.argv[16]]
#TRACES TO BE CROSS-CORRELATED AGAINST THE TEMPLATE
files2xcorr=[sys.argv[6],sys.argv[7]]


#STATION LOCATIONS- MUST BE IN THE SAME ORDER AS INPUTED IN FILES (ALL STATIONS)
s1=np.array([[-2.84,53.84]])#L2
s2=np.array([[-2.94,53.76]])#L4
s3=np.array([[-2.88,53.69]])#L6
s4=np.array([[-3.04,53.78]])#L9
#SEISMIC VELOCITY USED FOR INVERSION
velocity=5
#TEMPLATE START TIME AND LENGTH
stt=UTCDateTime("2018-10-27T10:55:26")
dtt=4
#CROSS-CORRELATION COEFFICIENT TRIGGER THRESHOLD
xcorr_threshold=0.35
tr=read()
tr1=read()

#4-STATION STA/LTA NETWORK COINCIDENCE TRIGGER - SETTINGS CAN BE ADJUSTED AFTER COINCIDENCE_TRIGGER (STA/LTA TYPE, TRIGGER-ON, TRIGGER-OFF, DATA, COINCIDENCE SUM, STA WINDOW SIZE, LTA WINDOW SIZE)
for filename in files:
	st += read(filename,starttime=starttr, endtime=end)
	tr_filt = st.copy()
	tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)
	st2 = tr_filt.copy()
	trig = coincidence_trigger("recstalta", 1.9, 1, st2, 12, sta=0.5, lta=2.5)

#2-STATION STA/LTA NETWORK COINCIDENCE TRIGGER- DETECTIONS WILL BE CROSS-CORRELATED FOR CONFIRMATION OF TRIGGER
for filename in files1:
	st+=read(filename,starttime=starttr, endtime=end)
	tr_filt=st.copy()
	tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)
	st3=tr_filt.copy()
	trig2=coincidence_trigger("recstalta",1.3,1,st3,6,sta=0.5,lta=2.5)


#EARTHQUAKE DETECTIONS FROM 4-STATION COINCIDENCE STA/LTA TRIGGER
stalta_earthquakes=[d['time'] for d in trig]

#DETECTIONS TIMES FROM 2-STATION NETWORK COINCIDENCE TRIGGER TO BE CROSS-CORRELATED
possquakes=[d['time'] for d in trig2]#

#SETTING UP LIST TO HOLD EARTHQUAKES DETECTED BY CROSS-CORRELATION
xcorr_earthquakes=[]
print('detections to be Xcorr')
print(possquakes)
counter=0
#SETTING UP LIST TO HOLD CROSS-CORRELATION COEFFCIENTS PRODUCED ON TRIGGERED EVENTS
xcorr_val=[]
xcorr_en_val=[]

#CROSS-CORRELATION PROCESS
for d in possquakes:
	counter=0
	total_coeff=[]
	total_coeff_en=[]
    #FIRST THE Z-COMP IS CROSS-CORRELATED
	st =read(sys.argv[5],starttime=d-1, endtime=d+4)
	tr =read(sys.argv[14],starttime=stt,endtime=stt+dtt)
	tr2 = copy.deepcopy(st)
	if len(tr2) > 0:
        #DETRENDING, TAPERING AND FILTERING DATA AND TEMPLATE
		tr.detrend(type='simple')
		tr.taper(0.01, type='hann', max_length=None, side='both')
		tr.filter('bandpass',freqmin=5, freqmax=45, corners=6)
        #ENVELOPE OF TEMPLATE
		tr_en = envelope(tr[0].data)
		tr2.merge(method=0, fill_value=0, interpolation_samples=0)
		tr2.detrend(type='simple')
		tr2.taper(0.01, type='hann', max_length=None, side='both')
		tr2.filter('bandpass',freqmin=5, freqmax=45, corners=6)
        #ENVELOPE OF DATA
		tr2_en = envelope(tr2[0].data)
        #CROSS-CORRELATION FUNCTION
		cc = xcorr(tr2[0].data, tr[0].data, tr2[0].stats.npts, demean=True, normalize=True, domain='freq')
        #TIME SHIFT WHEN MAXIMUM CROSS-CORRELATION COEFFICIENT ACHIEVED
		shift, value = xcorr_max(cc, abs_max=True)
        #ENVELOPE CROSS-CORRELATION
		cc_en = xcorr(tr2_en, tr_en, tr2[0].stats.npts, demean=True, normalize=True, domain='freq')
        #TIME SHIFT WHEN MAXIMUM CROSS-CORRELATION COEFFCIENT ACHIEVED WITH ENVELOPES
		shift_en, value_en = xcorr_max(cc_en, abs_max=True)
        #UTC DATE TIME WHEN MAXIMUM CROSS-CORRELATION  COEFFICIENT IS ACHIEVED
		times = tr2[0].times('utcdatetime')
		tbest=int(np.min([tr2[0].stats.npts-1,int((tr2[0].stats.npts/2.))+shift]))
		pick = times[tbest]
		tbest_en=int(np.min([tr2[0].stats.npts-1,int((tr2[0].stats.npts/2.))+shift_en]))
		pick_en = times[tbest_en]
		print ('station: %s temp:%s %s %s to %s - Best Xcorr: %s %.2f %.2f %s %.2f %.2f' % (st[0].stats.station,st[0].stats.channel,tr[0].stats.channel,d, d+dtt, pick, value, xcorr_threshold, pick_en, value_en, xcorr_threshold*1.5))
        #IF STATEMENT TO CHECK WHETHER THE CROSS-CORRELATION TRIGGER THRESHOLD IS EXCEEDED WITH EITHER NORMAL OR ENVELOPED 
		if np.abs(value) >= xcorr_threshold or np.abs(value_en) >= np.min([xcorr_threshold*1.5, 1]):
            #IF EXCEEDED ON THE Z-COMPONENT THEN N AND E COMPS WILL ALSO BE CROSS-CORRELATED
			counter=counter+1
			nocomp=counter
			total_coeff.append(np.abs(value))
			total_coeff_en.append(np.abs(value_en))
			print('Number of components triggered: %s' %(nocomp))
            #LOOP THAT CROSS-CORRELATED N AND E COMPS
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
                            #IF ALL THREE COMPONENTS HAVE EXCEEDED THE TRIGGER THRESHOLD, THE TIME FROM THE 2-STATION COINCIDENCE TRIGGER WILL BE RECORDED
							xcorr_val.append((np.abs(total_coeff)))
							xcorr_en_val.append((np.abs(total_coeff_en)))
							print('EARTHQUAKE DETECTED at: %s'%(d))
							xcorr_earthquakes.append(d)
#CREATING ADDITIONAL CROSS CORRELATION LIST IN ORDER TO REMOVE DUPLICATE DETECTIONS
xcorr_earthquakes_singles=[]
#LIST TO HOLD DETECTION TIMES WITH THE CROSS CORRELATION COEFFICIENTS PRODUCED
xcorr_results=xcorr_earthquakes + xcorr_val +  xcorr_en_val

#REMOVING DUPLICATE CROSS CORRELATION TRIGGER TIMES
for i in xcorr_earthquakes:
    if i not in xcorr_earthquakes_singles:
        xcorr_earthquakes_singles.append(i)

#LIST OF TOTAL EARTHQUAKES DETECTED FROM BOTH METHODS
Tquakes=[]
#ADDING 4-STATION NETWORK COINCIDENCE DETECTIONS TO TOTAL DETECTIONS
for i in stalta_earthquakes:
    Tquakes.append(i)
#ADDING CROSS CORRELATION DETECTIONS TO TOTAL DETECTIONS
for d in xcorr_earthquakes_singles:
    Tquakes.append(d)

#FRACKING SITE COORDINATES
fracking_loc=np.array([[-2.95,53.79]])
#STATION LOCATIONS ARRAY
s_locs=np.array([[s1],[s2],[s3],[s4]])
#STATION LOCATION LONGITUDES
s_lon=np.array([[s1[0][0]],[s2[0][0]],[s3[0][0]],[s4[0][0]]])
#STATION LOCATION LATITUDES
s_lat=np.array([[s1[0][1]],[s2[0][1]],[s3[0][1]],[s4[0][1]]])
#LAT/LON TO CARTESIAN SYSTEM WITH THE FRACKING SITE COORDS AS X,Y=(0,0)
X=(s_lon-fracking_loc[0][0])*40000*np.cos((s_lat+fracking_loc[0][1])*np.pi/360)/360
Y=(fracking_loc[0][1]-s_lat)*40000/360
s_z=np.array([[0],[0],[0],[0]])
#STATION LOCATIONS IN CARTESIAN SYSTEM (X,Y,Z)
s_clocs=np.column_stack((X,Y,s_z))
#HYPOCENTRE STARTING MODEL
hypo=np.array([[0],[0],[-2.5]])
#HYPOCENTRAL DISTANCE TO EACH STATION
station_x_hypo=(np.square(s_clocs[:,0]))+ (np.square((s_clocs[:,1])) + (np.square(s_clocs[:,2])))
station2hypo=np.power(station_x_hypo,0.5)


#STATION 4 AND 9, N AND Z COMPS
ncomp4=sys.argv[7]
ncomp9=sys.argv[13]
zcomp9=sys.argv[11]
zcomp4=sys.argv[5]

#LISTS TO HOLD EARTHQUAKE LOCATIONS AND MAGNITUDES
eq_locs=[]
eq_mags_bgs=[]
eq_mags_sarg=[]

#ALL DETECTIONS TO BE QUANTIFIED
for i in Tquakes:
    #MAGNITUDE CALCULATION USING N COMP FROM STATION 4 
	BGS_mag, o_sarg_mag=magnitude(ncomp4,i,2,station2hypo[1])#using station9 to calc mag
	eq_mags_bgs.append(BGS_mag)
	eq_mags_sarg.append(o_sarg_mag)
    #LIST TO HOLD P-WAVE ARRIVAL TIES
	p_arrivals=[]
    #DATA USED TO PICK P-WAVE ARRIVALS , 1 SECOND PREVIOUS TO TRIGGER TIME TO ENSURE P-WAVE IS IN THE DATA WINDOW
	zt = read(zcomp4,starttime=(UTCDateTime(i))-1, endtime=(UTCDateTime(i))+3) 
	zt1=read(zcomp9,starttime=(UTCDateTime(i))-1, endtime=(UTCDateTime(i))+3)
	tr_filt = zt.copy()[0]
	tr_filt2 = zt1.copy()[0]
	df=tr_filt.stats.sampling_rate
    #FILTERING DATA
	tr_filt.filter('bandpass', freqmin=5, freqmax=45, corners=6)
	st2 = tr_filt.copy()
	tr_filt2.filter('bandpass', freqmin=5, freqmax=45, corners=6)
	st3 = tr_filt2.copy()
    #P-WAVE PICK USING BAER FUNCTION 
	p_pick4, phase_info = pk_baer(st2.data, df,20, 20, 60, 12.0, 100, 100)
    #IF PICK IS LESS THAN 2 SAMPLES, MEANING IT HAS NOT FOUND THE P-WAVE, THEN A PICK IS ESTIMATED USING THE SET VELOCITY AND HYPOCENTRAL DISTANCE
	if p_pick4<2:
		p_pick4=random.uniform(station2hypo[1]/0.055,station2hypo[1]/0.045)#if no p-wave if found then a random pick is generated
	p_arrivals.append(p_pick4/100)
	p_pick9, phase_info = pk_baer(st3.data, df,20, 20, 60, 12.0, 100, 100)
    #IF PICK IS LESS THAN 2 SAMPLES, MEANING IT HAS NOT FOUND THE P-WAVE, THEN A PICK IS ESTIMATED USING THE SET VELOCITY AND HYPOCENTRAL DISTANCE
	if p_pick9<2:
		p_pick9=random.uniform(station2hypo[3]/0.055,station2hypo[3]/0.045)
	p_arrivals.append(p_pick9/100)
    #ESTIMATING P-WAVE ARRIVAL TIES FOR NOISY STATIONS L2 AND L6, USING THE DIFFERENCE IN HYPOCENTRAL DISTANCE FROM THESE STATIONS AND STATION 4
	D1=(p_pick4/100)+((station2hypo[0]-station2hypo[1])/velocity)
	D3=(p_pick4/100)+((station2hypo[2]-station2hypo[1])/velocity)
	p_arrivals.insert(0,D1)
	p_arrivals.insert(2,D3)
	Do=np.reshape(np.asarray(p_arrivals),(4,1))
    
#STARTING MODEL
	Mi=np.array([[0],[0],[-2.5],[-1],[0]])
#X,Y,Z DISTANCES FROM STARTING MODEL TO STATION
	R=np.subtract(s_clocs,np.transpose(Mi[0:3]))
#STRAIGHT LINE RAY PATHS
	Rh=(np.square(R[:,0]))+ (np.square((R[:,1])) + (np.square(R[:,2])))
	Rh=np.power(Rh,0.5)
	Rh=np.reshape(Rh,(4,1))
#RAY TRAVEL TIMES USING STRAIGHT RAY PATHS
	Tt=Rh/velocity
#PREDICTED P-WAVE ARRIVALS
	Dp=np.transpose(Tt + Mi[3,:])
	Dp=np.reshape(Dp,(4,1))
#CREATING CONNECTION MATRIX SETTING UP ITERATION COUNTER
	iterations=1
	G=np.zeros([len(s_clocs),4])
	misfit_array=[]
#INVERSION TO ITERATE 8 TIMES
	while iterations<8:
    #ARRIVAL TIME RESIDUAL
		Dres=Do-Dp
    #RESIDUAL MISFIT
		misfit=np.sum([np.square(Dres)])
		misfit_array.append(misfit)
    #CONNECTING MATRIX MADE FROM PARTIAL DERIVATES OF TIME WITH RESPECT TO X,Y,Z,T
		G1=-(np.reshape(R[:,0],(4,1))/velocity)*(np.power(Rh,-0.5))
		G2=-(np.reshape(R[:,1],(4,1))/velocity)*(np.power(Rh,-0.5))
		G3=-(np.reshape(R[:,2],(4,1))/velocity)*(np.power(Rh,-0.5))
		G4=np.ones([(len(s_clocs)),1])
		G=np.column_stack((G1,G2,G3,G4))
    #LAGRANGE MULTIPLIERS TO LOCK DEPTH TO STARTING MODEL DEPTH
		h=0
		F=np.array([[0],[0],[1],[0]])
		H1=np.column_stack([np.dot(np.transpose(G),G), F])
		H2= np.append(np.transpose(F), F[0])
		H2=H2.astype(np.float64)
		H=np.vstack((H1,H2))
		Ginv=np.linalg.inv(H)
    #PREIVOUS MODEL
		M_prev=Mi
    #CHANGE IN MODEL MATRIX 
		delta_model1=np.dot(np.transpose(G),Dres)
		delta_model1=np.append(delta_model1,F[0])
		delta_model=np.dot(Ginv,delta_model1)
    #NEW MODEL BY ADDING PREVIOUS MODEL AND CHANGE IN MODEL
		Mi=M_prev+np.reshape(delta_model,(5,1))
    #NEW RAY PATHS USING DISTANCE FROM STATIONS TO THE NEW MODEL
		R=np.subtract(s_clocs,np.transpose(Mi[0:3]))
		Rh=(np.square(R[:,0]))+ (np.square((R[:,1])) + (np.square(R[:,2])))
		Rh=np.power(Rh,0.5)
		Rh=np.reshape(Rh,(4,1))
    #NEW TRAVEL TIMES
		Tt=Rh/velocity
    #NEW P WAVE ARRIVAL TIMES
		Dp=np.transpose(Tt + Mi[3,:])
		Dp=np.reshape(Dp,(4,1))
		iterations=iterations+1
		if iterations>7:
            #ADDING EARTHQUAKE LOCATION TO LIST
			eq_locs.append(Mi)

#CONVERTING CARTESIAN EARTHQUAKE LOCATIONS TO LON/LAT
x_dist=[d[0] for d in eq_locs]
y_dist=[d[1] for d in eq_locs]
eq_times=[d[3] for d in eq_locs]
s_mlon=np.mean(s_lon)
s_mlat=np.mean(s_lat)
eq_lat=(np.add(y_dist,(1000/9)*(s_mlat)))/(1000/9)
eq_lon=s_mlon-(((np.multiply(360,x_dist))/40000)/np.cos((s_mlat+eq_lat)*np.pi/360))
eq_locs_lonlat=np.column_stack((eq_lon,eq_lat,eq_times))

#OUTPUTTING RESULTS TO TEXT FILES
with open("stalta_results.txt", "w") as output:
        output.write("\n")
        output.write(str(stalta_earthquakes))
        
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
        
