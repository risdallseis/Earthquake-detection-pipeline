# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 23:10:12 2019

@author: rory
"""

def conv2d(tr):
    import copy
    import numpy as np
    import matplotlib.pyplot as plt
    from obspy.core import read
    from obspy.signal.invsim import paz_to_freq_resp
    from obspy.core import UTCDateTime
    tr.detrend(type='simple')
    tr.taper(0.01, type='hann', max_length=None, side='both')
    tr.filter('bandpass',freqmin=5, freqmax=45, corners=6)
    pre_filt = (0.3, 0.5, 45.0, 50.0)
    seedresp = {'filename': 'RESP.LV.L00X.00.HHX',  # RESP filename
            # when using Trace/Stream.simulate() the "date" parameter can
            # also be omitted, and the starttime of the trace is then used.
            #'date': date,
            # Units to return response in ('DIS', 'VEL' or ACC)
            # Use DIS for WA calculation
            'units': 'DIS'}
    tr.simulate(paz_remove=None, pre_filt=pre_filt, remove_sensitivity=True, seedresp=seedresp)
    poles =  [-5.49779+5.60886j,-5.49779-5.60886j]
    zeros = [0j, 0j]
    gain = 1.
    sens = 1.0
    wa = {'poles': poles,
      'zeros': zeros,
      'gain': gain,
      'sensitivity': sens    }
    tr.simulate(paz_remove=None, paz_simulate=wa, simulate_sensitivity=True)
    trD=copy.deepcopy(tr)
    return trD