import sys
import os
import getopt
import shutil
import subprocess

from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd 


# Conversion factors

mp=1.6726e-24
b_fac_pluto=0.0458505
r_fac_pluto=2.149e+02
temp_fac_pluto=1
rho_fac_pluto=mp
time_fac_pluto=1.49597871e+08/3600 

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def read_cme_params(pluto_ini_file):
	""" Read pluto.ini file fo run and retrieve start time 
	of CME and duration. Return both in code units"""

	cme_params_dict = {'CME_START_TIME': 0.0, 'CME_DURATION': 0.0}

	with open(pluto_ini_file, 'r') as fp:
		line = fp.readline()
		while line:
			newline = line;
			for sub in cme_params_dict.keys():
				if (line.find(sub) != -1):
					line2 = line.strip().split(' ')
					last = len(line2) - 1
					cme_params_dict[sub] = line2[last]
			line = fp.readline()

	start = float(cme_params_dict['CME_START_TIME']) * time_fac_pluto
	dur   = float(cme_params_dict['CME_DURATION']) * time_fac_pluto

	return start, dur


# def read_obs_r0_dat(obs_file):
# 	""" read observer file at inner boundary and convert units to hours nT"""

# 	df = pd.read_table(obs_file,delim_whitespace=True, skiprows=1)
# 	df['time'] = df['time'] * time_fac_pluto
# 	df['Bp']   = df['Bp'] * b_fac_pluto

# 	return (df)

def read_obs_dat(obs_file):
	""" read observer file and convert units to hours nT"""
	df = pd.read_table(obs_file,delim_whitespace=True, skiprows=1)
	df['time'] = df['time'] * time_fac_pluto
	df['Bp']   = df['Bp'] * b_fac_pluto

	return (df)

def read_tracer_dat(tracer_file):
	""" read tracers file """
	df = pd.read_table(tracer_file,delim_whitespace=True, skiprows=0)
	df['time'] = df['time'] * time_fac_pluto
	df['tracer1'] = df['tracer1'] * r_fac_pluto
	df['tracer2'] = df['tracer2'] * r_fac_pluto

	return (df)

def plot_vars_at_r0(df, cme_start_time, cme_duration, filename):

	time = df['time']
	bpr0 = df['Bp']
	vrr0 = df['vr']
	npr0 = df['rho']
	tr0 = df['temp']

	cme_end_time = cme_start_time + cme_duration

	#########################################################################
	#
	# Find the time range we should use for the plot at the inner boundary
	#

	imax = np.argmax(vrr0)
	dt_h = time[imax] - time[(imax-1)] # estimate for the dt in hours

	time_max = time[imax]

	wwindow = 20.
	tmin = round(time_max) - wwindow
	tmax = round(time_max) + wwindow

	#########################################################################
	# 
	# Plot Vr, Bp, and np perturbations at the inner boundary 
	#

	f1 = figure(figsize=[15,8], num=1)

	ax1 = f1.add_subplot(221)
	time = np.array(time)
	ydata = np.array(vrr0) 
	ax1.plot(time, ydata, color = 'red')
	ax1.axvline(x = cme_start_time, color = 'grey', linestyle = '--')
	ax1.axvline(x = cme_end_time, color = 'grey', linestyle = '--')
	xlim([tmin,tmax])
	#xlabel(r'Time [hours]')
	ylabel(r'V [km/s]')
	title("Time Profile of Velocity Perturbation at Inner Boundary")
	ax1.grid(which='major',axis='y')

	ax1 = f1.add_subplot(222)
	ydata = np.array(bpr0) 
	ax1.plot(time, ydata, color = 'blue')
	ax1.axvline(x = cme_start_time, color = 'grey', linestyle = '--')
	ax1.axvline(x = cme_end_time, color = 'grey', linestyle = '--')
	xlim([tmin,tmax])
	xlabel(r'Time [hours]')
	ylabel(r'B$_{\rm \phi}$ [nT]')
	title("Time Profile of Magnetic Field Perturbation at Inner Boundary")
	ax1.grid(which='major',axis='y')

	ax1 = f1.add_subplot(223)
	ydata = np.array(npr0) 
	ax1.plot(time, ydata, color = 'green')
	ax1.axvline(x = cme_start_time, color = 'grey', linestyle = '--')
	ax1.axvline(x = cme_end_time, color = 'grey', linestyle = '--')
	xlim([tmin,tmax])
	xlabel(r'Time [hours]')
	ylabel(r'n [$cm^{-3}$]')
	title("Time Profile of Density Perturbation at Inner Boundary")
	ax1.grid(which='major',axis='y')
	plt.savefig(filename)
	#plt.show()
	plt.clf()
	plt.cla()
	plt.close()
	return


def plot_vars_at_1au(df, cme_start_time, cme_duration, filename):

	time  = df['time'] 
	bp1au = df['Bp']
	vr1au = df['vr']
	np1au = df['rho']
	t1au  = df['temp']

	cme_end_time = cme_start_time + cme_duration
	#########################################################################
	#
	# Find the time range we should use for the plot
	#

	imax = np.argmax(vr1au)
	
	time_max = time[imax] # This is the time that Vr has a maximum

	# find index of CME start time - start plotting a little before
	itmin, time_min = find_nearest(time, cme_start_time)

	itmin = itmin - 10

	itmax = len(time)-1   
	tmin = round(time[itmin])
	tmax = round(time[itmax])

	#########################################################################
	#
	# Plot the variables at 1AU
	#

	f1 = figure(figsize=[15,8], num=1)

	ax1 = f1.add_subplot(411)
	time_1au = np.array(time)
	ydata = np.array(vr1au) 
	ax1.plot(time_1au, ydata, color = 'red')
	ax1.axvline(x = cme_start_time, color = 'grey', linestyle = '--')
	ax1.axvline(x = cme_end_time, color = 'grey', linestyle = '--')
	ymin=min(ydata[itmin:itmax])
	ymax=max(ydata[itmin:itmax])     
	xlim([tmin,tmax])
	ylim([0.9*ymin, 1.1*ymax])
	xlabel(r'Time [hours]')
	ylabel(r'v$_{r}$ [km/s]')

    #--------------------------------------------------
    
    # a short digression to calculate the time of arrival (ToA) of 
    # the shock at 1 AU. This should probably be done somewhere 
    # else or, even better, wrapped up into a separate function...but

	# Ideally, this should probably include calculations of the 
	# jump conditions, but since it's a 1-D solution, we know that 
	# the solar wind is constant at 1 AU until the arrival of the disturbance. 
	# If it's a fast CME, then the first perturbation will be the arrival
	# of the forward shock. So, any reasonably sharp increase in 
	# solar wind speed should act to capture it. We'll call that 
	# threshold dv_crit. 

	dv_crit = 20.

	dv = ydata[1:itmax] - ydata[0:(itmax-1)]
	print("max(dv):",np.max(dv))
	ishock = np.min(np.where(dv > dv_crit))

	print("Shock time: ", time_1au[ishock])
	ax1.axvline(x = time_1au[ishock], color = 'red', linestyle = '--')
	ax1.text(time_1au[ishock],600,"Shock ToA: "+str(round(time_1au[ishock] - 200,2))+" hours")

	# And also calculate and print the ambient values 
	# of the solar wind

	iamb = np.min(np.where(time_1au > 210))
	print("V_amb: ", ydata[iamb])
	ax1.text(225,1.05*ydata[iamb],"$v_{amb}$: "+str(round(ydata[iamb],2))+" km/s")

    #--------------------------------------------------

	ax1 = f1.add_subplot(412)
	ydata = np.array(np1au) 
	ax1.plot(time_1au, ydata, color = 'green')
	ax1.axvline(x = cme_start_time, color = 'grey', linestyle = '--')
	ax1.axvline(x = cme_end_time, color = 'grey', linestyle = '--')
	yscale('log')
	ymin=min(ydata[itmin:itmax])
	ymax=max(ydata[itmin:itmax])        
	xlabel(r'Time [hours]')
	ylabel(r'n [$cm^{-3}$]')
	xlim([tmin,tmax])
	ylim([0.9*ymin, 1.1*ymax])
	ax1.axvline(x = time_1au[ishock], color = 'red', linestyle = '--')

	print("n_amb: ", ydata[iamb])
	ax1.text(225,1.1*ydata[iamb],"$n_{amb}$: "+str(round(ydata[iamb],2))+" cm$^{-3}$")

	ax1 = f1.add_subplot(413)
	ydata = np.array(bp1au) 
	ax1.plot(time_1au, ydata, color = 'blue')
	ax1.axvline(x = cme_start_time, color = 'grey', linestyle = '--')
	ax1.axvline(x = cme_end_time, color = 'grey', linestyle = '--')
	ymin=min(ydata[itmin:itmax])
	ymax=max(ydata[itmin:itmax])         
	xlabel(r'Time [hours]')
	ylabel(r'B$_{\rm \phi}$ [nT]')
	xlim([tmin,tmax])
	ylim([0.9*ymin, 1.1*ymax])
	ax1.axvline(x = time_1au[ishock], color = 'red', linestyle = '--')

	print("Bp_amb: ", ydata[iamb])
	ax1.text(225,1.1*ydata[iamb],"$Bp_{amb}$: "+str(round(ydata[iamb],2))+" nT")

	ax1 = f1.add_subplot(414)
	ydata = np.array(t1au) 
	ax1.plot(time_1au, ydata, color = 'magenta')
	ax1.axvline(x = cme_start_time, color = 'grey', linestyle = '--')
	ax1.axvline(x = cme_end_time, color = 'grey', linestyle = '--')
	yscale('log')
	ymin=min(ydata[itmin:itmax])
	ymax=max(ydata[itmin:itmax])     
	xlabel(r'Time [hours]')
	ylabel(r'T [K]')
	xlim([tmin,tmax])
	ylim([0.9*ymin, 1.1*ymax])
	ax1.axvline(x = time_1au[ishock], color = 'red', linestyle = '--')

	print("T_amb: ", ydata[iamb])
	ax1.text(225,1.1*ydata[iamb],"$T_{amb}$: "+str(round(ydata[iamb],2))+" K")

	plt.savefig(filename)
	#plt.show()
	plt.clf()
	plt.cla()
	plt.close()

	return

def plot_tracers(df, cme_start_time, cme_duration, filename):

	""" Plot tracers location as a function of time"""

	cme_end_time = cme_start_time + cme_duration

	time  = df['time'] 
	tracer1 = df['tracer1']
	tracer2 = df['tracer2']

	f1 = figure(figsize=[15,5], num=1)
	ax1 = f1.add_subplot(111)
	ax1.plot(time, tracer1, label = 'tracer1')
	ax1.plot(time, tracer2, label = 'tracer2')
	ax1.axvline(x = cme_start_time, color = 'grey', linestyle = '--')
	ax1.axvline(x = cme_end_time, color = 'grey', linestyle = '--')
	xlabel(r'Time [hours]')
	ylabel(r'R [Rs]')
	title('Tracers Position as a function of time')
	legend()
	plt.savefig(filename)
	#plt.show()
	plt.clf()
	plt.cla()
	plt.close()
