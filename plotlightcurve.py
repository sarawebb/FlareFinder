import numpy as np 
import os 
import matplotlib.pyplot as plt


def plot_lc(path):
	
	for filename in os.listdir(path):
		if filename.endswith('.lc'):
			MJD = []
			MAG = []
			file_dir = path + filename
			mjd, mag , emag, ulmag = np.loadtxt(file_dir, unpack= True)
			mjd_min = np.min(mjd)
			mjd_mod = mjd - mjd_min
			#print(mjd_mod)
		
			plt.plot(mjd_mod ,mag, 'bo')
			plt.show()
			
	
	
testing = plot_lc('/home/swebb/oz100/flare_project/testing_lc/3hr_151218_msystembis1_3-400_h400_e0/')
