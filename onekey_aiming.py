import numpy as np
import os
import re
import shutil
from Deviation_aiming_new3 import *
from scipy.optimize import curve_fit
from sys import path
from python_postprocessing import *
from Open_CSPERB import *
from Open_CSPERB_plots import *
from HC import *
from Tube_materials import Inconel740H
from Flux_reader import *
from cal_sun import *
from scipy import interpolate

class one_key_start:
	def __init__(self, folder,num_bundle,num_fp,r_height,r_diameter,bins,tower_h,phi,elevation,DNI,D0):
		self.folder=folder
		self.r_diameter=r_diameter # receiver diameter,m
		self.r_height=r_height # receiver height,m
		self.num_bundle=num_bundle # number of panels
		self.num_fp=num_fp # number of flow paths
		self.bins=bins # receiver vertical bin
		self.tower_h=tower_h # tower height,m
		self.phi=phi # azimuth angle, deg
		self.elevation=elevation # elevation angle, deg
		self.DNI=DNI # DNI,W/m2
		self.D0=D0 # pipe outer diameter, mm
		self.csv_aiming='%s/pos_and_aiming_new.csv' % self.folder
		self.csv_trimmed='%s/pos_and_aiming.csv'%self.folder
		self.latitude=34.85
		self.dis_delta=5 # discretising the declination angle
		self.dis_omega=25 # discretising the solar hour angle
		
	def run_SOLSTICE(self,dni,phi,elevation,att_factor,num_rays,csv): # the input is not a solstice style
		# transfer into SOLSTICE convention
		phi=270.-phi
		if phi > 360.:
			phi=phi-360.
		print phi, elevation,dni
		vtk_path='%s/vtk'%self.folder
		if os.path.exists(vtk_path):
			shutil.rmtree(vtk_path)
		# replace keywords in SOLSTICE.py
		file_path='%s/SOLSTICE.py' % self.folder
		old_file=file_path
		fopen=open(old_file,'r') 
		w_str=""
		for line in fopen:
			if re.search('dni_1=',line):
				line = 'dni_1=%s' % (dni) + '\n'
				w_str+=line
			elif re.search('azimuth_1=',line):
				line = 'azimuth_1=%s' % (phi) + '\n'
				w_str+=line
			elif re.search('zenith_1=',line):
				line = 'zenith_1=%s' % (elevation) + '\n'
				w_str+=line
			elif re.search('att_factor_1=',line):
				line = 'att_factor_1=%s' % (att_factor) + '\n'
				w_str+=line
			elif re.search('mainfolder_1=',line):
				line = "mainfolder_1='%s'" % (self.folder) + '\n'
				w_str+=line
			elif re.search('csv_1=',line):
				line = "csv_1='%s'" % (csv) + '\n'
				w_str+=line
			elif re.search('num_rays_1=',line):
				line = "num_rays_1=%s" % (num_rays) + '\n'
				w_str+=line
			elif re.search('r_cyl_1=',line):
				line = "r_cyl_1=%s" % (self.r_diameter/2.) + '\n'
				w_str+=line
			elif re.search('h_cyl_1=',line):
				line = "h_cyl_1=%s" % (self.r_height) + '\n'
				w_str+=line
			elif re.search('tower_h_1=',line):
				line = "tower_h_1=%s" % (self.tower_h) + '\n'
				w_str+=line
			elif re.search('num_bundle_1=',line):
				line = "num_bundle_1=%s" % (self.num_bundle) + '\n'
				w_str+=line
			else:
				w_str+=line
		wopen=open(old_file,'w')
		wopen.write(w_str)
		fopen.close()
		wopen.close()
		os.system('python2 %s/SOLSTICE.py ' % self.folder)
		
	def attenuation(self,csv): # calculate the attenuation factor
		hst_info=np.loadtxt(csv,delimiter=',', skiprows=2)
		foc=hst_info[:,3]
		# to get the attenuation factor
		def func(x, b):
			return np.exp(-b * x)
		def fun_two(x):
			return 0.99321-0.0001176*x+1.97e-8*x**2
		xdata = np.linspace(0, np.max(foc), np.max(foc)*100)
		y = fun_two(xdata)
		ydata = y
		popt, pcov = curve_fit(func, xdata, ydata)
		y2 = [func(i, popt[0]) for i in xdata]
		att_factor =popt[0]
		return att_factor
		
	def HT_model(self,T_amb,V_wind): # receiver thermal model
		rec = Cyl_receiver(radius=0.5*self.r_diameter, height=self.r_height, n_banks=self.num_bundle, n_elems=50, D_tubes_o=self.D0/1000., D_tubes_i=self.D0/1000.-2.*1.2e-3, 
		  abs_t=0.98, ems_t=0.91, k_coating=1.2, D_coating_o=self.D0/1000.+45e-6)
		Strt=rec.flow_path(option='cmvNib%s'%self.num_fp,fluxmap_file=self.folder+'/flux-table.csv')
		rec.balance(HC=Na(), material=Inconel740H(), T_in=520+273.15, T_out=740+273.15, T_amb=T_amb+273.15, h_conv_ext='SK', filesave=self.folder+'/flux-table',air_velocity=V_wind)
		flux_limits_file='%s/201015_N07740_thermoElasticPeakFlux_velocity/N07740_OD%s_WT1.20_peakFlux_vel.csv'%(self.folder,round(self.D0,2))
		results,aiming_results,vel_max=tower_receiver_plots(files=self.folder+'/flux-table', efficiency=False, maps_3D=False, flux_map=False, flow_paths=True,saveloc=None, billboard=False, flux_limits_file=flux_limits_file,C_aiming=self.C_aiming)
		return results,aiming_results,Strt

	def aiming_loop(self,C_aiming,Exp,A_f): # the aiming strategy loop: optical + thermal
		# the input for optical modelling
		self.C_aiming=C_aiming
		print C_aiming
		print Exp
		print A_f
		att_factor=self.attenuation(self.csv_trimmed)
		aiming(self.folder,self.r_height,self.r_diameter,C_aiming,self.csv_trimmed,self.tower_h,self.num_bundle,Exp,A_f) # change aiming points 
		self.run_SOLSTICE(dni=self.DNI,phi=self.phi,elevation=self.elevation,att_factor=att_factor,num_rays=5000000,csv=self.csv_aiming) # optical simulation
		eta,q_results,eta_exc_intec=proces_raw_results('%s/vtk/simul'% self.folder,'%s/vtk'% self.folder) # optical postprocessing
		eff_interception=eta/eta_exc_intec
		print 'Interception efficiency: ' + str(eff_interception)
		read_data(self.folder,self.r_height,self.r_diameter,self.num_bundle,self.bins,flux_file=True) # read flux map
		results,aiming_results,Strt=self.HT_model(20.,0.) # thermal simulation
		#print aiming_results[0]
		print aiming_results[1]
		
		return aiming_results,eff_interception,Strt
	
	def search_algorithm(self): # parametric study of aiming extent
		C_aiming=np.zeros(self.num_bundle) # aiming extent
		C_aiming[:]=0. # equatorial aiming
		Exp=np.zeros(self.num_bundle) # shape exponent
		Exp[:]=1.5 # initialising exponent
		A_f=np.zeros(self.num_bundle) # asymmetry factor
		A_f[:int(0.25*self.num_bundle)]=A_f[int(0.75*self.num_bundle):]=0.67 # initiliasing asymmetry factor
		A_f[int(0.25*self.num_bundle):int(0.75*self.num_bundle)]=0.33
		aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)
		
		# to output eff_interception at equatorial aiming
		savedir='%s/Equatorial_interception.csv' % self.folder
		f=open(savedir, 'w')
		f.write('%s'%eff_interception)
		f.write("\n")
		f.close()
		
		# search algorithm
		C_aiming[:]=0.5
		aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)
		while np.all(aiming_results[0])==False and np.all(C_aiming<1.):
			for i in range(self.num_bundle):
				if aiming_results[0][i]==False:
					C_aiming[Strt[i]]+=0.05
					if Strt[i]==self.num_bundle-1:
						C_aiming[0]+=0.05
					else:
						C_aiming[Strt[i]+1]+=0.05
					C_aiming[Strt[i]-1]+=0.05
			aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)
			print C_aiming
			print aiming_results[0]
		
		# to output the results of search algorithm
		savedir='%s/Search_results.csv' % self.folder
		f=open(savedir, 'w')
		f.write(",".join(map(str, C_aiming)))
		f.write("\n")
		f.close()
		
	def fit_algorithm(self): # optimisation of shape exponent and asymmetry factor
		if os.path.exists('%s/output'%self.folder):
			shutil.rmtree('%s/output'%self.folder)
		# run the results after search algorithm
		savedir='%s/Search_results.csv' % self.folder 
		C_aiming=np.loadtxt(savedir,delimiter=',', skiprows=0)
		Exp=np.zeros(self.num_bundle)
		Exp[:]=1.5
		A_f=np.zeros(self.num_bundle)
		A_f[:int(0.25*self.num_bundle)]=A_f[int(0.75*self.num_bundle):]=0.67
		A_f[int(0.25*self.num_bundle):int(0.75*self.num_bundle)]=0.33
		aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)

		#The optimisation uses Dakota from Sandia. Dakota is activated by running the dakota input file GA.in from os.system.
		#Dakota has an interface with python code, which is programmed in Dakota_interface.py.
		#There can be a series of optimisation in the fit algorithm, hence, another code Dakota_aiming.py is established. The current code
		#functions as changing the keywords in Dakota_aiming.py ready for different optimisations.
		
		# initialisation for the optimisation
		file_path='%s/Dakota_aiming.py' % self.folder
		old_file=file_path
		fopen=open(old_file,'r') 
		w_str=""
		for line in fopen:
			if re.search('num_bundle_2',line):
				line = "	num_bundle_2=num_bundle=%s\n" % self.num_bundle
				w_str+=line
			elif re.search('r_height_2',line):
				line = "	r_height_2=r_height=%s\n" % self.r_height
				w_str+=line
			elif re.search('r_diameter_2',line):
				line = "	r_diameter_2=r_diameter=%s\n" % self.r_diameter
				w_str+=line
			elif re.search('tower_h_2',line):
				line = "	tower_h_2=tower_h=%s\n" % self.tower_h
				w_str+=line
			elif re.search('num_fp_2',line):
				line = "	num_fp_2=num_fp=%s\n" % self.num_fp
				w_str+=line
			elif re.search('D0_2',line):
				line = "	D0_2=D0=%s\n" % self.D0
				w_str+=line
			elif re.search('phi_1',line):
				line = "	phi_1=phi=%s\n" % self.phi
				w_str+=line
			elif re.search('elevation_1',line):
				line = "	elevation_1=elevation=%s\n" % self.elevation
				w_str+=line
			elif re.search('DNI_1',line):
				line = "	DNI_1=DNI=%s\n" % self.DNI
				w_str+=line
			else:
				w_str+=line
		wopen=open(old_file,'w')
		wopen.write(w_str)
		fopen.close()
		wopen.close()
		
		# Strt comes from the receiver thermal model and shows the relationship between flow path and bank index.
		for i in range(self.num_fp):
			fopen=open(old_file,'r') 
			w_str=""
			for line in fopen:
				if re.search('C_aiming_%s_1'%(i+1),line):
					line = "	C_aiming_%s_1=C_aiming[%s]=%s\n" % (i+1,Strt[2*i],C_aiming[Strt[2*i]])
					w_str+=line
				elif re.search('Exp_%s_1'%(i+1),line):
					line = "	Exp_%s_1=Exp[%s]=%s\n" % (i+1,Strt[2*i],Exp[Strt[2*i]])
					w_str+=line
				elif re.search('A_f_%s_1'%(i+1),line):
					line = "	A_f_%s_1=A_f[%s]=%s\n" % (i+1,Strt[2*i],A_f[Strt[2*i]])
					w_str+=line
				elif re.search('C_aiming_%s_2'%(i+1),line):
					line = "	C_aiming_%s_2=C_aiming[%s]=%s\n" % (i+1,Strt[2*i+1],C_aiming[Strt[2*i+1]])
					w_str+=line
				elif re.search('Exp_%s_2'%(i+1),line):
					line = "	Exp_%s_2=Exp[%s]=%s\n" % (i+1,Strt[2*i+1],Exp[Strt[2*i+1]])
					w_str+=line
				elif re.search('A_f_%s_2'%(i+1),line):
					line = "	A_f_%s_2=A_f[%s]=%s\n" % (i+1,Strt[2*i+1],A_f[Strt[2*i+1]])
					w_str+=line
				else:
					w_str+=line
			wopen=open(old_file,'w')
			wopen.write(w_str)
			fopen.close()
			wopen.close()
		
		# This is for situation when the number of banks is larger than 16.
		for i in range(self.num_fp,12):
			fopen=open(old_file,'r') 
			w_str=""
			for line in fopen:
				if re.search('C_aiming_%s_1'%(i+1),line):
					line = "	C_aiming_%s_1=C_aiming[%s]=%s\n" % (i+1,0,C_aiming[0])
					w_str+=line
				elif re.search('Exp_%s_1'%(i+1),line):
					line = "	Exp_%s_1=Exp[%s]=%s\n" % (i+1,0,Exp[0])
					w_str+=line
				elif re.search('A_f_%s_1'%(i+1),line):
					line = "	A_f_%s_1=A_f[%s]=%s\n" % (i+1,0,A_f[0])
					w_str+=line
				elif re.search('C_aiming_%s_2'%(i+1),line):
					line = "	C_aiming_%s_2=C_aiming[%s]=%s\n" % (i+1,0,C_aiming[0])
					w_str+=line
				elif re.search('Exp_%s_2'%(i+1),line):
					line = "	Exp_%s_2=Exp[%s]=%s\n" % (i+1,0,Exp[0])
					w_str+=line
				elif re.search('A_f_%s_2'%(i+1),line):
					line = "	A_f_%s_2=A_f[%s]=%s\n" % (i+1,0,A_f[0])
					w_str+=line
				else:
					w_str+=line
			wopen=open(old_file,'w')
			wopen.write(w_str)
			fopen.close()
			wopen.close()
		
		# ready for optimisation
		for i in range(self.num_bundle):
			if aiming_results[2][i]==0.:
				continue
			print i,int(i/2)+1,i%2+1
			Tube=Strt
			
			file_path='%s/Dakota_aiming.py' % self.folder
			old_file=file_path
			fopen=open(old_file,'r') 
			w_str=""
			for line in fopen:
				if re.search('Exp_%s_%s' % (int(i/2)+1,i%2+1),line):
					line = "	Exp_%s_%s=Exp[%s]=%s\n" % (int(i/2)+1,i%2+1,Tube[i],'x[0]') # change the optimised variables to x[0],x[1]..
					w_str+=line
				elif re.search('A_f_%s_%s' % (int(i/2)+1,i%2+1),line):
					line = "	A_f_%s_%s=A_f[%s]=%s\n" % (int(i/2)+1,i%2+1,Tube[i],'x[1]')
					w_str+=line
				elif re.search('gx=',line):
					line = "	gx=aiming_results[2][%s]\n" % i # change the objective function to the crossover extent at the current tube bank
					w_str+=line
				else:
					w_str+=line
			wopen=open(old_file,'w')
			wopen.write(w_str)
			fopen.close()
			wopen.close()
			
			# run the optimisation
			if i%2==0:
				os.system('dakota -i opt_patternsearch1.in -o opt_patternsearch%s.out'% i) 
			else:
				os.system('dakota -i opt_patternsearch2.in -o opt_patternsearch%s.out'% i) 
			if not os.path.exists('%s/output'%self.folder):
				os.makedirs('%s/output'%self.folder)
			# move the optimisation results to the output folder
			if os.path.exists('%s/opt_patternsearch%s.out'%(self.folder,i)):
				shutil.copy('%s/opt_patternsearch%s.out' % (self.folder,i), '%s/output'%self.folder)
				os.remove('%s/opt_patternsearch%s.out'% (self.folder,i))
			
			# to read the output of Dakota
			out_file='%s/output/opt_patternsearch%s.out' % (self.folder,i)
			info=np.genfromtxt(out_file, delimiter=' ', dtype=str,usecols=(0))
			Exp[Tube[i]]=info[-10]
			A_f[Tube[i]]=info[-9]
			
			# update the optimisation results to Dakota_aiming.py
			file_path='%s/Dakota_aiming.py' % self.folder
			old_file=file_path
			fopen=open(old_file,'r') 
			w_str=""
			for line in fopen:
				if re.search('Exp_%s_%s' % (int(i/2)+1,i%2+1),line):
					line = "	Exp_%s_%s=Exp[%s]=%s\n" % (int(i/2)+1,i%2+1,Tube[i],info[-10])
					w_str+=line
				elif re.search('A_f_%s_%s' % (int(i/2)+1,i%2+1),line):
					line = "	A_f_%s_%s=A_f[%s]=%s\n" % (int(i/2)+1,i%2+1,Tube[i],info[-9])
					w_str+=line
				else:
					w_str+=line
			wopen=open(old_file,'w')
			wopen.write(w_str)
			fopen.close()
			wopen.close()
			
		# to output the results of fit algorithm
		savedir='%s/Fit_results.csv' % self.folder
		f=open(savedir, 'w')
		f.write(",".join(map(str, C_aiming)))
		f.write("\n")
		f.write(",".join(map(str, Exp)))
		f.write("\n")
		f.write(",".join(map(str, A_f)))
		f.write("\n")
		f.close()
		
	def adjust_algorithm(self,varying_DNI=False): # re-adjust the aiming extent if crossover still exists
		# to read the results of fit algorithm
		savedir='%s/Fit_results.csv' % self.folder
		if varying_DNI==True: # for DNI deviated from the clear-sky DNI
			savedir='%s/Interpolated_results.csv' % self.folder
		Results=np.loadtxt(savedir,delimiter=',', skiprows=0)
		C_aiming=Results[0]
		Exp=Results[1]
		A_f=Results[2]
		aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)
		print aiming_results[2][:16]
		print eff_interception
		# adjustment algorithm
		while (not np.all(aiming_results[2]<10.)):
			for i in range(self.num_bundle):
				if aiming_results[1][i]==False:
					C_aiming[Strt[i]]+=0.05
					if Strt[i]==self.num_bundle-1:
						C_aiming[0]+=0.05
					else:
						C_aiming[Strt[i]+1]+=0.05
					C_aiming[Strt[i]-1]+=0.05
			aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)
			print C_aiming
			print Exp
			print A_f
			print aiming_results[2]
			print eff_interception
			
		# to output the final results
		if varying_DNI==True:
			savedir='%s/Interpolated_final_results.csv' % self.folder
		else:
			savedir='%s/Final_results.csv' % self.folder
		f=open(savedir, 'w')
		f.write(",".join(map(str, C_aiming)))
		f.write("\n")
		f.write(",".join(map(str, Exp)))
		f.write("\n")
		f.write(",".join(map(str, A_f)))
		f.write("\n")
		f.close()
		print 'done'
		
		# defocusing if aiming extent > 1
		safety_aiming=[50.,50.,self.tower_h+self.r_height*0.5]
		if (not np.all(C_aiming<=1.)):
			title=np.array(['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z', 'm', 'm', 'm', 'm', 'm', 'm', 'm'])
			pos_and_aiming_new=np.array([])
			csv='%s/pos_and_aiming_new.csv' % path[0]
			hst_info=np.loadtxt(csv,delimiter=',', skiprows=2) 
			num_hst=hst_info.size/7
			
			for i in range(num_hst):
				if hst_info[i,6]>(self.tower_h+self.r_height) or hst_info[i,6]<self.tower_h:
					print hst_info[i,6],i
					hst_info[i,4]=safety_aiming[0]
					hst_info[i,5]=safety_aiming[1]
					hst_info[i,6]=safety_aiming[2]
					#print hst_info[i,:]
			
			pos_and_aiming_new=np.append(pos_and_aiming_new, hst_info)
			pos_and_aiming_new=pos_and_aiming_new.reshape(len(pos_and_aiming_new)/7, 7)
			np.savetxt(csv, pos_and_aiming_new, fmt='%s', delimiter=',')
			
			aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)
			print aiming_results[2][:16]
			print eff_interception
		
	def annual_aiming(self): # the get look-up tables after running MDBA optimisation at discretised points
		Delta=np.linspace(-23.45, 23.45, num=self.dis_delta) # decliation angle
		Omega=np.linspace(-180., 180., num=self.dis_omega) # solar hour angle
		sun=SunPosition() # import from cal_sun.py coded by Ye Wang
		for i in range(self.dis_delta):
			for j in range(self.dis_omega):
				# to get the azimuth angle, elevation angle and DNI
				daytime,sunrise=sun.solarhour(Delta[i], self.latitude)
				theta=sun.zenith(self.latitude, Delta[i], Omega[j])
				phi=sun.azimuth(self.latitude, theta, Delta[i], Omega[j])
				elevation=90.-theta
				DNI=self.get_I_Meinel(elevation) # clear-sky DNI
				print Delta[i],Omega[j],elevation,DNI
				if elevation<=15.:
					continue
				# create subfolders to store the results
				subfolder='%s/process/Aiming_%s_%s'%(path[0],round(Delta[i],2),round(Omega[j],2))
				if not os.path.exists(subfolder):
					os.makedirs(subfolder)
				else:
					continue
				self.phi=phi
				self.elevation=elevation
				self.DNI=DNI
				# MDBA optimisation
				self.search_algorithm()
				self.fit_algorithm()
				self.adjust_algorithm()
				# move the results to the subfolder
				shutil.copy('%s/Equatorial_interception.csv' % self.folder, subfolder)
				shutil.copy('%s/Final_results.csv' % self.folder, subfolder)
				shutil.copy('%s/Fit_results.csv' % self.folder, subfolder)
				shutil.copy('%s/flux-table_flux_fp.png' % self.folder, subfolder)
				shutil.copy('%s/Search_results.csv' % self.folder, subfolder)
				if os.path.exists('%s/output'%subfolder):
					shutil.rmtree('%s/output'%subfolder)
				shutil.copytree('%s/output' % self.folder, '%s/output'%subfolder)
		
		# generate the look-up tables
		E=np.arange((self.dis_delta+1)*(self.dis_omega+1),dtype=float).reshape(self.dis_delta+1,self.dis_omega+1)
		E[:,:]=0.
		Delta=np.linspace(-23.45, 23.45, num=self.dis_delta)
		Omega=np.linspace(-180., 180., num=self.dis_omega)
		E[0,1:]=Omega
		E[1:,0]=Delta
		
		# to read and output the aiming extent array
		for l in range(self.num_bundle):
			E[1:,1:]=0.5
			for i in range(self.dis_delta):
				for j in range(self.dis_omega):
					subfolder='%s/process/Aiming_%s_%s'%(path[0],round(Delta[i],2),round(Omega[j],2))
					if os.path.exists(subfolder):
						C_aiming=np.loadtxt('%s/Final_results.csv' % subfolder,delimiter=',', skiprows=0)
						E[i+1,j+1]=C_aiming[0][l]
					# to output 
					savedir='%s/process/Aiming_%s.csv' % (self.folder,l) # look-up tables for aming extents
					f=open(savedir, 'w')
					for k in range(self.dis_delta+1):
						f.write(",".join(map(str, E[k])))
						f.write("\n")
					f.close()
		
		# to read and output the shape exponent array
		for l in range(self.num_bundle):
			E[1:,1:]=1.5
			for i in range(self.dis_delta):
				for j in range(self.dis_omega):
					subfolder='%s/process/Aiming_%s_%s' % (path[0],round(Delta[i],2),round(Omega[j],2))
					if os.path.exists(subfolder):
						Exp=np.loadtxt('%s/Final_results.csv' % subfolder,delimiter=',', skiprows=0)
						E[i+1,j+1]=Exp[1][l]
					# to output 
					savedir='%s/process/Shape_%s.csv' % (self.folder,l)
					f=open(savedir, 'w')
					for k in range(self.dis_delta+1):
						f.write(",".join(map(str, E[k])))
						f.write("\n")
					f.close()
		
		# to read and output the asymmetry factor array
		for l in range(self.num_bundle):
			if l<4 or l>11:
				E[1:,1:]=0.33
			else:
				E[1:,1:]=0.67
			for i in range(self.dis_delta):
				for j in range(self.dis_omega):
					subfolder='%s/process/Aiming_%s_%s' % (path[0],round(Delta[i],2),round(Omega[j],2))
					if os.path.exists(subfolder):
						A_f=np.loadtxt('%s/Final_results.csv' % subfolder,delimiter=',', skiprows=0)
						E[i+1,j+1]=A_f[2][l]
					# to output 
					savedir='%s/process/Asymmetry_%s.csv' % (self.folder,l)
					f=open(savedir, 'w')
					for k in range(self.dis_delta+1):
						f.write(",".join(map(str, E[k])))
						f.write("\n")
					f.close()
	
	def interpolation(self,delta,omega): # to get the aiming variables by using interpolation from look-up tables
		Delta=np.linspace(-23.45, 23.45, num=self.dis_delta)
		Omega=np.linspace(-180., 180., num=self.dis_omega)
		C_aiming=np.zeros(self.num_bundle)
		Exp=np.zeros(self.num_bundle)
		A_f=np.zeros(self.num_bundle)
		# to interpolate the aiming extent
		for l in range(num_bundle):
			savedir='%s/process/Aiming_%s.csv' % (self.folder,l)
			z=np.loadtxt(savedir,delimiter=',', skiprows=0)
			z2=z[1:,1:].transpose()
			f = interpolate.interp2d(Delta,Omega,z2, kind='linear')
			znew = f(delta,omega)
			C_aiming[l]=znew[0]
		
		# to interpolate the shape exponent
		for l in range(num_bundle):
			savedir='%s/process/Shape_%s.csv' % (self.folder,l)
			z=np.loadtxt(savedir,delimiter=',', skiprows=0)
			z2=z[1:,1:].transpose()
			f = interpolate.interp2d(Delta,Omega,z2, kind='linear')
			znew = f(delta,omega)
			Exp[l]=znew[0]
		
		# to interpolate the asymmetry factor
		for l in range(num_bundle):
			savedir='%s/process/Asymmetry_%s.csv' % (self.folder,l)
			z=np.loadtxt(savedir,delimiter=',', skiprows=0)
			z2=z[1:,1:].transpose()
			f = interpolate.interp2d(Delta,Omega,z2, kind='linear')
			znew = f(delta,omega)
			A_f[l]=znew[0]
		
		# to output
		savedir='%s/Interpolated_results.csv' % self.folder
		f=open(savedir, 'w')
		f.write(",".join(map(str, C_aiming)))
		f.write("\n")
		f.write(",".join(map(str, Exp)))
		f.write("\n")
		f.write(",".join(map(str, A_f)))
		f.write("\n")
		f.close()
		print 'done'
		'''
		# do the full MDBA optimisation
		sun=SunPosition()
		daytime,sunrise=sun.solarhour(delta, self.latitude)
		theta=sun.zenith(self.latitude, delta, omega)
		phi=sun.azimuth(self.latitude, theta, delta, omega)
		elevation=90.-theta
		DNI=self.get_I_Meinel(elevation)
		print phi,elevation,DNI,round(delta,2), round(omega,2)
		self.phi=phi
		self.elevation=elevation
		self.DNI=DNI
		subfolder='%s/Aiming_%s_%s'%(path[0],round(delta,2),round(omega,2))
		if not os.path.exists(subfolder):
			os.makedirs(subfolder)
		# MDBA optimisation
		self.search_algorithm()
		self.fit_algorithm()
		self.adjust_algorithm()
		shutil.copy('%s/Equatorial_interception.csv' % self.folder, subfolder)
		shutil.copy('%s/Final_results.csv' % self.folder, subfolder)
		shutil.copy('%s/Fit_results.csv' % self.folder, subfolder)
		shutil.copy('%s/flux-table_flux_fp.png' % self.folder, subfolder)
		shutil.copy('%s/Search_results.csv' % self.folder, subfolder)
		shutil.copy('%s/Interpolated_results.csv' % self.folder, subfolder)
		if os.path.exists('%s/output'%subfolder):
			shutil.rmtree('%s/output'%subfolder)
			shutil.copytree('%s/output' % self.folder, '%s/output'%subfolder)
		'''
	
	def check(self): # to validate the interpolated results
		# choose 5*10 discretised points to validate the interpolated results
		Delta=np.linspace(-23., 23., num=5)
		Omega=np.linspace(-60., 60., num=10)
		Validation_results=np.full((5,10), False, dtype=bool) # boolean array for validation results
		sun=SunPosition()
		#A_over=np.array([0.,20.])
		for i in range(5):
			for j in range(10):
				daytime,sunrise=sun.solarhour(Delta[i], self.latitude)
				theta=sun.zenith(self.latitude, Delta[i], Omega[j])
				phi=sun.azimuth(self.latitude, theta, Delta[i], Omega[j])
				elevation=90.-theta
				DNI=self.get_I_Meinel(elevation)
				print i,j,phi,elevation,DNI,round(Delta[i],2), round(Omega[j],2)
				self.interpolation(delta=Delta[i],omega=Omega[j])
				self.phi=phi
				self.elevation=elevation
				self.DNI=DNI
				Results=np.loadtxt('%s/Interpolated_results.csv' % self.folder,delimiter=',', skiprows=0)
				C_aiming=Results[0]
				Exp=Results[1]
				A_f=Results[2]
				aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)
				if np.all(aiming_results[2]<10.)==True: # Very small deviation is allowed here.
					Validation_results[i,j]=True
		print Validation_results
	
	def varying_DNI(self,delta,omega,factor): # for cases when real DNI deviates from clear-sky values
		# get the interpolated results
		self.interpolation(delta,omega) 
		sun=SunPosition()
		daytime,sunrise=sun.solarhour(delta, self.latitude)
		theta=sun.zenith(self.latitude, delta, omega)
		phi=sun.azimuth(self.latitude, theta, delta, omega)
		elevation=90.-theta
		DNI=self.get_I_Meinel(elevation) # clear-sky DNI
		self.phi=phi
		self.elevation=elevation
		self.DNI=DNI*factor # the real DNI
		self.adjust_algorithm(varying_DNI=True) # use the adjustment algorithm to re-adjust the aiming extent
	
	def get_I_Meinel(self,elevation): # Meinel clear-sky model
		I0=1363.
		zenith=90.-elevation
		AM=1./np.cos(zenith/180.*np.pi)
		I=I0*0.7**(AM**0.678)
		return I
	
	def New_search_algorithm(self):  # the net one-key algorithm to get the optimised aiming points
		C_aiming=np.zeros(self.num_bundle) # aiming extent
		C_aiming[:]=0. # equatorial aiming
		Exp=np.zeros(self.num_bundle) # shape exponent
		Exp[:]=2.0 # initialising exponent
		A_f=np.zeros(self.num_bundle) # asymmetry factor
		if self.num_bundle/self.num_fp == 1:
			A_f[:]=0.75
		elif self.num_bundle/self.num_fp == 2:
			A_f[:int(0.25*self.num_bundle)]=A_f[int(0.75*self.num_bundle):]=0.67
			A_f[int(0.25*self.num_bundle):int(0.75*self.num_bundle)]=0.33
		
		# new search algorithm
		C_aiming[:]=0.5
		aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)
		gap=0.05
		while np.all(aiming_results[1])==False and np.all(C_aiming<1.):
			# for E
			C_aiming_old=np.ones(self.num_bundle)
			C_aiming_old[:]=C_aiming[:]
			for i in range(self.num_bundle):
				if aiming_results[1][i]==False:
					C_aiming[Strt[i]]+=gap
					if Strt[i]==self.num_bundle-1:
						C_aiming[0]+=gap
					else:
						C_aiming[Strt[i]+1]+=gap
					C_aiming[Strt[i]-1]+=gap
				# for A
				if A_f[Strt[i]]>0.5:
					if (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])<-0.1:
						A_f[Strt[i]]+=0.02
					elif (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])>0.1:
						A_f[Strt[i]]-=0.02
				else:
					if (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])<-0.1:
						A_f[Strt[i]]-=0.02
					elif (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])>0.1:
						A_f[Strt[i]]+=0.02
				# for S
				if aiming_results[5][i]>0.55:
					Exp[Strt[i]]-=0.2
				elif aiming_results[5][i]<0.45:
					Exp[Strt[i]]+=0.2
			C_aiming[C_aiming-C_aiming_old>gap]=C_aiming_old[C_aiming-C_aiming_old>gap]+gap
			aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)
	
if __name__=='__main__':
	from sys import path
	folder=path[0]
	num_bundle=16
	r_height=24.
	r_diameter=16.
	bins=50
	tower_h=175.
	phi=0.0
	elevation=55.15
	DNI=980.0
	num_fp=num_bundle/2
	D0=60.33
	Model=one_key_start(folder,num_bundle,num_fp,r_height,r_diameter,bins,tower_h,phi,elevation,DNI,D0)	
	#Model.search_algorithm()
	#Model.fit_algorithm()
	#Model.adjust_algorithm()
	#Model.annual_aiming()
	#Model.check()
	#Model.interpolation(delta=23.0,omega=20.0)
	#Model.varying_DNI(delta=23.0,omega=20.0,factor=1.2)
	Model.New_search_algorithm()
