# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 15:43:10 2019

@author: rory
"""
import numpy 

def lonlat2cart(s_lon,s_lat,fracking_loc):
       R = numpy.float64(6371000)  # in meters
       longitude = s_lon
       latitude = s_lat
       longitude_f=fracking_loc[0][0]
       latitude_f=fracking_loc[0][1]
       numpy.array([X]) = (R * numpy.cos(longitude) * numpy.sin(latitude))-(R * numpy.cos(longitude_f) * numpy.sin(latitude_f))
       numpy.array([Y]) = (R * numpy.sin(latitude) * numpy.sin(longitude))-(R * numpy.sin(latitude_f) * numpy.sin(longitude_f))
       return numpy.array([X,Y])


lat1=53.79#fracking
lon1=-2.95
lat2=53.84
lon2=-2.84

dx = (lon2-lon1)*40000*np.cos((lat1+lat2)*np.pi/360)/360
dy = (lat1-lat2)*40000/360