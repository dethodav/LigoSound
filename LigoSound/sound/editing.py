# -*- coding: utf-8 -*-
# Copyright (C) Derek Davis (2017)
#
# This file is part of the Ligo Sounds python package.
#
# Ligo Sounds is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ligo Sounds is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Ligo Sounds.  If not, see <http://www.gnu.org/licenses/>.


"""Methods and utilties for performing audio pipline scans
"""


# This is where we will put together the file that processes 
# and writes the audio file

# Format is sound(gps, channel, frametype, outdir, ASDloc, shift, stretch, Aweight, amp)

# Add each filtering option as own function, which is then plugged into sound()

from  gwpy.timeseries import TimeSeries
import scipy.io.wavfile as wav
import numpy as np
import scipy as sp
import scipy.signal as sig
import pickle

def fshift(timeseries,shift_size,method='push'):
	"""Frequency shift the spectrum of the Timeseries.
	Parameters
	----------
	shift_size:`float`  
		size and sign of frequency shift in Hz. 
	method:'string', optional
	       method to prefrom shift
	       default is push
	       other option is hilbert          
	"""

	data = timeseries.value
	samp_rate = timeseries.sample_rate.value
 
	if (method=='push'):
	    time_length = len(data)/float(samp_rate)
	    df = 1.0/time_length
	    nbins = int(shift_size/df)

	    freq_rep = np.fft.rfft(data)
	    shifted_freq = np.zeros(len(freq_rep),dtype=complex)
	    for i in range(0,len(freq_rep)-1):
		    if 0<(i-nbins)<len(freq_rep):
			   shifted_freq[i]=freq_rep[i-nbins]
	    output = np.fft.irfft(shifted_freq)
	    out_real = np.real(output)

	if (method=='hilbert'):
	    if (fshift < 0):
		timeseries_high = timeseries.highpass( (shift_size * -1.0) )
		data = timeseries_high.value

	    dt = 1.0/samp_rate
	    N = len(data)
	    t = numpy.arange(0, N)
	    out_real = (signal.hilbert(data)*numpy.exp(2j*numpy.pi*shift_size*dt*t)).real

	out = TimeSeries(out_real,sample_rate=samp_rate)
	out.__metadata_finalize__(timeseries)
	out._unit = timeseries.unit
	del out.times

	return out

def wavwrite(timeseries,file_name,rate=4096,amp=.1):
	"""Prepares the timeseries for audio and writes 
	to a .wav file.
	Parameters
	----------
	file_name: `str`
	    name of file to be written.
	
	rate: `float`, optional, default=4096
	    rate in Hz of the .wav file.
	amp: `float`, optional, default=.1
	    maximum amplitude of .wav file.
	See Also
	--------
	scipy.io.wavfile.write
	    for details on the write process. 
	"""

        if timeseries.sample_rate.value != rate:
            print "resampling"
	    timeseries_resamp = timeseries.resample(rate)
        else:
            timeseries_resamp = timeseries
	timeseries_normal  = amp * timeseries_resamp.value / (max(abs(timeseries_resamp.value)))

	wav.write(file_name,rate,timeseries_normal)

def time_expand_central_freq(timeseries,factor,central_freq=0.0):
	"""Changes the time length of  a timeseries while preserving
	   a specified frequency. Other frequencies will experience 
	   frequency modulation. 
	Parameters
	----------
	factor: 'float'
	    the factor that the timeseries will be 
	    lengthened in time.
	central_freq: 'float',optional, default=0.0
	    the specified frequency to preserve during
	    the change in time frame.
	See Also
	--------
	TimeSeries.fshift
	    for details on the frequency shifting  process. 
	"""    


	samp_rate = timeseries.sample_rate.value
	samp_rate_out = samp_rate * 1.0 / factor

	out = timeseries
	out.sample_rate = samp_rate_out
	del out.times
	if central_freq != 0.0:
	    shift_factor = central_freq * (1.0 - 1.0 / factor)
	    out = fshift(out,shift_factor)

	return out


def invA(f):
    #This designs the filter to be used in the filtering
    kA = 12194
    f1 = 20.6
    f2 = 107.7
    f3 = 737.9
    f4 = 12194
    return  ( (f**2+f1**2) * np.sqrt((f**2+f2**2) * (f**2+f3**2)) * (f**2+f4**2) ) / (kA**2*f**4) 

def invAWeight(TS,cutoff = 20):
    data = TS.value
    samp_rate = TS.sample_rate.value
    df = samp_rate / len(data)
    
    asd = np.fft.rfft(data)
    
    fi = 0

    for i in range(len(asd)):
        if (fi > cutoff):
            asd[i] = asd[i] * invA(fi*2*np.pi)
        else:
            asd[i] = 0.0
        fi += df
    time_out = np.fft.irfft(asd)
    out = TimeSeries(time_out, sample_rate = samp_rate)
    return out

def cutoff(timeseries, highfreq):
    timeseries_filtered = timeseries.lowpass(highfreq)
    timeseries_filtered = timeseries_filtered.resample(2*highfreq+100)
    return timeseries_filtered

def sound(gps, channel, duration, outdir, frame=None, ASDloc=None, 
          lpass=None, hpass=None, shift=0, stretch=1, 
          Aweight=False, maxamp=.1, whiten=True, center=False):
    print "begin"
    time = float(gps)

    if center == True:
        start = time - float(duration)/2.
        end = time + float(duration)/2.
    else: 
        start = time
        end = time + float(duration)

    if frame == None:
        data = TimeSeries.get(channel, start, end)
    else:
        data = TimeSeries.find(channel, start, end,
                               frametype=frame) 

    found = False
    base = 0
    close = 1e10
    if ASDloc != None:
	    asd_dict = pickle.load( open( ASDloc, "rb" ) )
		#finds closest asd to use
            for ts in asd_dict:
                k = time - ts
		if (close > k > 0):
                    base = ts
                    close = k
	            found = True
    if found == False:
        ASD = data.asd(4,2)
        print "own ASD"
    else:
        ASD = asd_dict[base]
        print base
    print "data found"
	#Now do filtering
    if whiten == True:
	   data = data.whiten(4,2,asd=ASD)

        #bandpass with sample rate reduction
    if lpass != None:
	    data = cutoff(data,float(lpass))
    if hpass != None:
        data = data.highpass(float(hpass))
	if (shift != 0):
	    data = fshift(data,float(shift))
    if Aweight == True:
	    data = invAWeight(data,cutoff=20)
    if stretch != 0:
        data = time_expand_central_freq(data,float(stretch))
    #set up file name
    path = outdir + '/' + gps  + '.wav'
    wavwrite(data,path,amp=maxamp,rate=data.sample_rate.value)
