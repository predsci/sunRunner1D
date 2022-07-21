#!/Users/michal/opt/miniconda3/envs/psi/bin/python



import sys
import os
import getopt
#import glob
#import time
#import datetime
import shutil
import subprocess
from util import *
#
# Write the welcome banner                                               # 
#
print('################################################')
print('################################################')
print(' ')
print(' Welcome to sunRunner1D 0.00001 (04/07/21) ')
print(' This is a 1D MHD calculation with a simple CME ')
print(' You will need to provide a <run name>, e.g. run001 ')
print(' Number of grid points: ')
print(' 	low - 500, medium - 1000, large - 2000 ')
print(' Outer boundary (in Rs units) default is ~257 Rs')
print(' You will also need to provide the CME parameters at 30Rs ')
print(' Currently the CME has a sin^2 shape but the ')
print(' magnetic field pulse can be set to have a sin shape ')
print(' with half the period so that is has an alternating sign ')
print(' User provided pulse parameters: ')
print(' Velocity pulse, default 2000 km/s ')
print(' Magnetic Field (Bp) pulse, default 600 nT ')
print(' Density pulse, default is 0.0 cm^-3 ')
print(' Duration of pulse, default is 10 hours ')
print(' Profile shape for magnetic field pulse (integer): ')
print('    0 = sin^2 (default) ')
print('    1 = sin with half the period ')
print(' All output files are saved in: runs/<run_name>/output ')
print(' For help use: ')
print(' ./sunRunner1D.py -h')
print(' Usage: ')
print(' ./sunRunner1D.py --run <run_name> --grid <grid_size> --r1 <R1> --t0 <t0> --rho0 <rho0> --v0 <v0> --bp0 <bp0> --dv <del_V> --db <del_Bp> --dn <del_rho> --dur <duration> --prof <profile>')
print(' Example Usage: ')
print(' ./sunRunner1D.py --run run001 --grid medium --r1 260 --dv 4000 --db 500 --dn 0 --dur 12 --prof 0')
print('################################################')
print('################################################')

def main(argv):


	template_file='cfg/pluto_cfg.ini'

	run_file = 'cfg/pluto.ini'

	##########################################################################
	#
	# check to see if we already have this file if so rename with a time stamp
	#

	#file_list = glob.glob(run_file)
	#ts = time.time()
	#readable_ts = datetime.datetime.fromtimestamp(ts).isoformat()
	#if len(file_list) > 0:
	#	os.rename(run_file, run_file+str(readable_ts))

	# To convert from PLUTO time units to hours

	time_fac_pluto = 1.49597871e+08/3600

	# To convert from PLUTO magnetic field units to nT

	b_fac_pluto = 0.0458505

	# To convert from AU to Rs multiply by this
	r_fac_pluto=2.149e+02

	# velocity in PLUTO is already in km/s
	# density in PLUTO is already in cm^-3

	# default values 

	run_name='run001'
	
	grid='low'
	r1=1.2 * r_fac_pluto # out boundary ~ 257 Rs 

	del_V = 2000.0
	del_B = 600.0 # in nT will be converted below to PLUTO units
	del_RHO=0.0
	dur = 10.0
	prof = 0

	# default background values 
	T_0 = 1.0e6
	BP_0 = 400.0 * b_fac_pluto # this is the default value in nT will be converted below to PLUTO units
	V_0 = 300.0
	RHO_0 = 600.0 

	dict_grid={'low': 500, 'medium': 1000, 'large': 2000}
	ngrid=dict_grid[grid]
	# Need to come up with letters for rho0 v0 and bp0 can not be r,v,b,n,d
	# use rho0 - c
	# use v0 - s	
	# use bp0 - m

	
	try:
		opts, args = getopt.getopt(argv,"h:r:g:l:t:c:s:m:v:b:n:d:p:",["help=","run=","grid=","r1=","t0=","rho0=","v0=","bp0=","dv=","db=","dn=","dur=","profile="])
	except getopt.GetoptError:
		print("\nWrong key word please use:\n")		
		print('sunRunner1D.py --run <run_name> --grid <grid_size> --r1 <grind_extent> --t0 <t0> --rho0 <rho0> --v0 <v0> --bp0 <bp0> --dv <del_V> --db <del_Bp> --dn <del_rho> --dur <duration> --profile <profile>')
		sys.exit(2)


	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print('opt ', opt)
			#print('sunRunner1D.py -r <run_name> -g <grid_size> -t <t0> -c <rho0> -s <v0> -m <bp0> -v <del_V> -b <del_Bp> -n <del_rho> -d <duration> -p <profile>')
			print('sunRunner1D.py --run <run_name> --grid <grid_size> --t0 <t0> --rho0 <rho0> --v0 <v0> --bp0 <bp0> --dv <del_V> --db <del_Bp> --dn <del_rho> --dur <duration> --profile <profile>')	
			sys.exit()
		elif opt in ("-r","--run"):
			run_name = arg
		elif opt in("-g","--grid"):
			grid = arg
			ngrid = dict_grid[grid]
		elif opt in("-l","--r1"):
			r1 = float(arg)
		elif opt in ("-t", "--t0"):
			T_0 = float(arg)
		elif opt in ("-c", "--rho0"):
			RHO_0 = float(arg)
		elif opt in ("-s", "--v0"):
			V_0 = float(arg)
		elif opt in ("-m", "--bp0"):
			BP_0 = float(arg)			
		elif opt in ("-v", "--dv"):
			del_V = float(arg)
		elif opt in ("-b", "--db"):
			del_B = float(arg)
		elif opt in ("-n", "--dn"):
			del_RHO = float(arg)
		elif opt in ("-d", "--dur"):
			dur = float(arg)
		elif opt in ("-p", "--prof"):
			prof = int(arg)
		else:
			print("Wrong key word, please use ./sunRunner1D.py --help to see a list of keywords")
			sys.exit(2)

	BP_0 = BP_0 / b_fac_pluto
	del_B = del_B / b_fac_pluto
	dur = dur / time_fac_pluto

	r1 = r1 / r_fac_pluto

	##########################################################################
	#
	# Read template pluto_cfg.ini file and write pluto.ini with run parameters
	#

	dict_replace = {'T_0': T_0, 'RHO_0': RHO_0, 'V_0': V_0, 'BP_0': BP_0, 'RHO_PERT': del_RHO, 'X1-grid': ngrid, 
	'V_PERT': del_V, 'BP_PERT': del_B, 'CME_DURATION': dur, 'PROFILE': prof}

	out_file = open(run_file, 'w')

	# Reading template file line by line and writing out line by line 

	with open(template_file, 'r') as fp:

		line = fp.readline()
		while line:
			newline=line
			for sub in dict_replace.keys():
				if (line.find(sub) != -1):

					line2 = line.strip().split(' ')
					last = len(line2)-1
					line2[last] = dict_replace[sub]
					newline=""
					for eachword in line2:
						newline = newline + ' '+str(eachword)
						newline = newline.lstrip()
					newline = newline+'\n'
					if (sub == 'X1-grid'):
						newline='X1-grid    1  0.13 '+str(dict_replace[sub])+' u  '+str(r1)+'\n'

			out_file.writelines(newline)

			line = fp.readline()

	out_file.writelines('\n')		

	out_file.close()



	##########################################################################
	#
	# set up the working directory for the PLUTO run
	#

	wdir = run_name
	if not os.path.exists('runs'):
		os.mkdir('runs')
	if not os.path.exists('runs/'+wdir):
		os.mkdir('runs/'+wdir)

	#########################################################################
	#
	# set up the output directory for the PLUTO run
	#

	pluto_out = 'output'

	if not os.path.exists('runs/'+wdir+'/'+pluto_out):
		os.mkdir('runs/'+wdir+'/'+pluto_out)

	#########################################################################
	#
	# move pluto.ini and copy definitions.h and pluto executable 
	# to PLUTO run directory 
	#

	src_path = os.getcwd()+'/cfg/'
	dst_path = os.getcwd() + '/runs/'+wdir

	pluto_ini_org = src_path+'pluto.ini'
	pluto_def_org = src_path+'definitions.h'
	pluto_exe_org = os.getcwd()+'/bin/pluto'

	pluto_ini_dst = dst_path+'/pluto.ini'
	pluto_def_dst = dst_path+'/definitions.h'
	pluto_exe_dst = dst_path+'/pluto'	


	shutil.move(    pluto_ini_org, pluto_ini_dst)
	shutil.copyfile(pluto_def_org, pluto_def_dst)
	shutil.copyfile(pluto_exe_org, pluto_exe_dst)

	#########################################################################
	#
	# Run Pluto
	#

   
	os.chdir('runs/')
	os.system('chmod -R 777 *')
	os.chdir(wdir)

	print('\n\nRunning Pluto...For Progress see:\n\n')
	print(os.getcwd()+'/out.txt')


	with open('out.txt','w+') as fout:
		with open('err.txt','w+') as ferr:
			out=subprocess.call(["./pluto"],stdout=fout,stderr=ferr)
   
	print('\n\n Run Completed Successfully \n\n')
	print(' Output Files saved at: \n')
	print(os.getcwd()+'\n')

	#########################################################################
	#
	# make plots
	#

	cme_start_time, cme_duration = read_cme_params(pluto_ini_file = pluto_ini_dst)

	# plots at inner boundary
	df_r0 = read_obs_dat(obs_file = dst_path+'/'+pluto_out+'/obs_r0.dat')

	plot_file = dst_path+'/'+'vars_r0.png'

	print('\n Saving Plot of Variables at Inner boundary', plot_file)
	
	plot_vars_at_r0(df_r0, cme_start_time, cme_duration, plot_file)

	
	# plots at 1AU

	df_1au = read_obs_dat(obs_file = dst_path+'/'+pluto_out+'/obs_1au.dat')

	plot_file = dst_path+'/'+'vars_1au.png'

	print('\nSaving Plot of Variables at R=1AU', plot_file)

	plot_vars_at_1au(df_1au, cme_start_time, cme_duration, plot_file)

	# tracers plot

	tracer_file = dst_path+'/'+pluto_out+'/tracers.dat'

	df_tracers = read_tracer_dat(tracer_file)
	
	plot_file = dst_path+'/'+'tracers.png'

	print('\nSaving Tracers Plot: ', plot_file,'\n')

	plot_tracers(df_tracers, cme_start_time, cme_duration, plot_file)

	############################################################################

if __name__ == "__main__":
   main(sys.argv[1:])