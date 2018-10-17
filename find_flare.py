import numpy as np
import os
from os.path import expanduser
import time
import datetime
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
from astropy.table import Table, Column, join 
from scipy.signal import wiener
from scipy import signal
from astropy.io import fits
import matplotlib
import glob
matplotlib.rcParams.update({'font.size':18})
matplotlib.rcParams.update({'font.family':'serif'})
from scipy.signal import savgol_filter


def FINDflare(flux, error, mjd_array,  N1=3, N2=1, N3=3, avg_std = False, std_window=7): 
	
	'''
	The algorithm for local changes due to flares defined by
    	S. W. Chang et al. (2015), Eqn. 3a-d
   	http://arxiv.org/abs/1510.01005
	 
	Adpated code from jradavenport's appalossa for kepler flares. Adpated for DWF data with high 
	cadence and upperlimits for deep detections. 
	
	
	Parametres: 
	------------
	flux: 	numpy array
		data to search over
	error:	numpy array
		errors to corresponding data 
	N1: 	int, optional 
		Coefficient from original paper - default is 3
		How many times above thes stddev is requried.
	N2: 	int, optional 
		Coefficient from original paper - default is 1 
		How many times above the stddev and uncertanity is required 
		NOTE: May need to be tweaked to work with fainted magnitudes with higher uncerts. 
	N3: 	int, optional 
		coefficient form the original paper - default is 3 
		The number of consecutive points required to flare as a flare 
		NOTE: 3 is a reasonable number for DWF 
	avg_std: bool, optional 
		Should the sigma in this data be computed by the median of the rolling().std? (default is false)
		
	std_window: float, optional
		if avg_std = True, how big of a window should it use? 
		Default is 25 data points. 
		NOTE: DWF data may not have 25 data points for all detections. 
	
	'''
	
	#Find the Median value of the flux of the flare. IMPORTANT make sure that the flux is for ONE night ONLY for DWF 
	med_i = np.nanmedian(flux)
	#print(flux)
	print('The median flux over lightcurve =' + str(med_i))

	if avg_std is False: 
		sig_i = np.nanstd(flux)  
		print('The standard deviation over the flux = : ' + str(sig_i))#the staddev of the window 
	else: 
		sig_i = np.nanmedian(pd.Series(flux).rolling(std_window, centre=True).std())
	
	ca = flux - med_i 
	cb = np.abs(flux-med_i)/sig_i 
	cc = np.abs(flux - med_i - error) / sig_i 
	#print(ca, cb, cc)
	## Apply the cuts from eqn 3a-c
	
	#### because our flux's are
	ctmp = np.where((ca < 0) & (cb > N1) & (cc > N2))
	#print(ctmp) 
	cindx = np.zeros_like(flux)
	cindx[ctmp] = 1 
	#print(cindx) 
	
	
	# need to find cumulative number of points that pass ctmp 
	# Count in reverse 
	
	ConM = np.zeros_like(flux)
	# this requires a full pass thrus the data -> bottleneck 
	for k in range(2, len(flux)):
		ConM[-k]= cindx[-k]*(ConM[-(k-1)] + cindx[-k] )
	print('Flare event detected in array at points, descending count from start point to finish point: ' + str(ConM))
	
	#these only defined between dl[i] and dr[i] 
	# find flare start where values in ConM switch from 0 to >=N3 
	istart_i = np.where((ConM[1:] >= N3)& (ConM[0:-1] - ConM[1:] <0))[0] +1 
	print('Flare begins on point: ' + str(istart_i))
	#print( istart_i)
	#Use the value of ConM to determine how many points away stop is 
	
	istop_i = istart_i + (ConM[istart_i] -1)
	istart_i = np.array(istart_i, dtype = 'int') 
	istop_i = np.array(istop_i, dtype = 'int') 
	print('Flare stops on point: ' + str(istop_i))

	
	min_day = np.min(mjd_array)
	test = np.sum(ConM)
	print('sum of ConM : ' + str(test))
	#flare_time = mjd_array - min_day
	#total_days = np.sum(flare_time[int(istart_i):int(istop_i + 1)])
	#total_flare_time = total_days * 60 * 60 * 24 
	#print('Duration of Flare detected in seconds: ' + str(total_flare_time))
	#p= np.trapz(new_array[int(istart_i):int(istop_i + 1)], x=(time * 60 * 60 * 24))
	#print(p)
# Testing getting flux arrays 

flare_path = '/home/swebb/oz100/flare_project/testing_lc/3hr_151218_msystembis1_3-400_h400_e0/' 

for filename in os.listdir(flare_path):
	if filename.endswith('.lc'): 
		mjd = []
		mag = []
		emag = []
		ulmag = []
		file_dir = flare_path + filename 
		fp_i = open(file_dir)
		next(fp_i)
		for row in fp_i:
			row = row.strip()
			columns = row.split()
			mjd.append(float(columns[:][0]))
			mag.append(float(columns[:][1]))
			emag.append(float(columns[:][2]))
			ulmag.append(float(columns[:][3]))
			
		#print(mjd)
		FINDflare(mag, emag, mjd,  N1=2, N2=1, N3=3, avg_std = False, std_window=7)
		
		#flares_found = Table()
		#flares_found['candidate name'] = filename
		#flares_found['flare_starts'] = istart_i 
		#flares_found['flare_stops'] = istop_i
		#flares_
		#FINDtime(mjd, mag)
			
	
