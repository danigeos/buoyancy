#!/usr/bin/python
#D. Garcia-Castellanos,2019
print('\n\n########################## Analytical slab buoyancy ############################')
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import os, sys, math, time

import matplotlib as mpl
from random import randint

print('Number of arguments: %s' %(len(sys.argv)))	
print('Argument List: %s' %(str(sys.argv)))

g=9.8 #m s-2
alpha=1e-6 #m2 s-1  Diffusivity
Mytosec=1e6*365.25*24*3600


#------------------------------------------
#Compositional density contrast at the LAB (kg m-3):
#12 (100km tecton) 
#25 (200km proton) 
#33 (100km proton) 
#64 (100km archon)
#68 (200km archon)
Drho_comp=40. #Choose a value here.
vel=10. #Convergence rate [mm/yr]
drho_dT=.1175 #Density dependence with temperature, kg m-3 K-1
zm=35e3 #Moho depth (top pf plate), m
zl=200e3 #Base of the Lithosphere, LAB, m
Tm=600  #Moho temperature C
Ta=1300 #Asthenosphere temperature C
sigma=30. #Subduction angle in degrees
tmax=30 #Final time in Myr

force_plot_min=-20e12
force_plot_max=+20e12
#-------------------------------------

H=zl-zm 
Dl=H/4
bb=H*H/8./alpha

nit=101


if (len(sys.argv)>=2): 
	Drho_comp = float(sys.argv[1])
	print('\tDrho_composit= %.1f [kg/m3]' %Drho_comp)
if (len(sys.argv)>=3): 
	vel = float(sys.argv[2])
	print('\tVelocity = %.1f [km/Myr]' %Drho_comp)
	
vel *= 1e3/Mytosec 
sigma *= np.pi/180

vel_arr=		np.linspace(.1, 5, nit)
time_arr=		np.linspace(0., tmax, nit)
Fby_adv_temp=	np.linspace(0, 0, nit)
Fby_adv_comp=	np.linspace(0, 0, nit)
Fby_dif=		np.linspace(0, 0, nit)
Fby_total=		np.linspace(0, 0, nit)


for i in range(0, nit):
	print('\n________________________________________________________')
	#vel_arr[i]*1e3/Mytosec
	time = time_arr[i]*Mytosec #Dl/vel
	Dl=time*vel
	time1 = H/math.sin(sigma)/vel
	timeL = time
	if (time>time1): 
		timeL=time1
		print('\nTime longer than Stage 1 (question results for diffusion)')
	print("i=%d  vel=%.2f km/My    time=%.2f My  Dl=%.1f km" %(i, vel*Mytosec/1e3, time/Mytosec, Dl/1e3))
	#Critical Drho_comp that is overcome by advection:
	Drho_comp_crit = 1./2. * drho_dT * (Ta-Tm)
	dFsp_adv_temp_dt_ini = g * Drho_comp_crit * H * vel
	#dFsp_adv_temp_dt = + dFsp_adv_temp_dt_ini * math.pow((1 - time/time1), 2)
	#Average density difference between LM and asthenosphere at LAB
	#Drho_comp_avg = Drho_comp - Drho_comp_crit*math.sin(sigma)/H*vel*timeL
	#print('Drho_comp_avg = %.2e' %(Drho_comp_avg))
	#dFsp_adv_comp_dt = - Drho_comp_avg * g * H * vel
	#print('Buoyant_force_per_My= %+.1e + %+.1e = %+.1e [N/m/My]' %((dFsp_adv_temp_dt)*Mytosec, (dFsp_adv_comp_dt)*Mytosec, (dFsp_adv_temp_dt+dFsp_adv_comp_dt)*Mytosec))

	#Advection above zl:
	Fby_adv_temp[i] = - dFsp_adv_temp_dt_ini * time1 / 3. * (1-math.pow(1-timeL/time1,3))

	#Advection below zl:
	#Stage 1:
	Fby_adv_comp[i] = + g*H*Drho_comp*vel*timeL - 1./2.*g*Drho_comp_crit*math.sin(sigma) * math.pow(vel*timeL,2)
	#Stage 2:
	if (time>time1): 
		Fby_adv_comp[i] += g * (Drho_comp - Drho_comp_crit) * H * vel * (time-time1)

	print('\tBuoyancy by advection above zl = %+.2e [N/m] (blue)'  %(Fby_adv_temp[i]))
	print('\tBuoyancy by advection below zl = %+.2e [N/m] (green)' %(Fby_adv_comp[i]))

	#Diffusion 1st approach
	#Diffusion lowermost triangle (Stage 1+2):
	#ff = 4./3/np.pi**.5 * g * alpha**.5 * drho_dT * (Ta-Tm)/5 
	#Fby_dif[i] = ff / 5 / H * math.sin(sigma) * math.pow(vel,2.) * ((-2*time-3*timeL) * math.pow(time-timeL,3./2) + 2*math.pow(time,5./2))
	#Diffusion above triangle (Stage 2):
	#if (time>time1): 
	#	Fby_dif[i] += ff * math.pow(time-time1, 3./2) * vel

	#Diffusion below LAB - 2nd approach (simpler): 
	#ff = 4. * g * drho_dT * (Ta-Tm) * bb * alpha / H * vel / time1
	#expmint = math.exp(-time/bb)
	#Fby_dif[i] = ff * ( timeL*timeL/2. - bb*timeL - bb*bb*(expmint-1) ) 
	#if (time>time1): 
	#	expmint1 = math.exp(-time1/bb)
	#	explust1 = math.exp(+time1/bb)
	#	Fby_dif[i] += \
	#		+ ff * bb*explust1 * (expmint*(time+bb)-expmint1*(time1+bb)) \
	#		- ff * bb*bb*explust1*(expmint-expmint1) \
	#		- ff * time1 * ( time+bb*explust1*expmint - time1 - bb ) 

	#Diffusion below LAB - 3rd approach (simplest and most accurate!): 
	#ff = 4. * alpha * g * (Ta-Tm) / H * drho_dT * bb * vel
	ff = 1./2. * g * drho_dT * (Ta-Tm) * H * vel
	#print('\tFactor = %+.1e' %ff)
	expmint = math.exp(-time/bb)
	Fby_dif[i] = ff * ( time + bb*(expmint-1) ) 

	#Diffusion above LAB:
	#MISSING!!

	print('\tBuoyancy by diffusion below zl = %+.2e [N/m] (red)' %Fby_dif[i])


	Fby_total[i]=Fby_adv_temp[i]+Fby_adv_comp[i]+Fby_dif[i]
	print('\tTOTAL buoyant force (>0 =>upwards) = %+.1e [N/m] (GREY)' %Fby_total[i])

fig=plt.cla()
plt.hlines(0, 0, tmax, 'k', 'dotted', '')
#plt.hlines(+1e13, 0, tmax, 'k', 'dotted', '')
#plt.hlines(-1e13, 0, tmax, 'k', 'dotted', '')
plt.plot(time_arr,Fby_adv_temp, 'b', linewidth=1.5)
plt.plot(time_arr,Fby_adv_comp, 'g', linewidth=1.5)
plt.plot(time_arr,Fby_dif, 'r', linewidth=1.5)
plt.plot(time_arr,Fby_total, 'grey', linewidth=3)
plt.ylim([force_plot_min,force_plot_max])
#plt.xlabel("vel (km/My)")
plt.xlabel("time (My)")
plt.ylabel("buoyancy (N/m)")
line_1, = plt.plot([1,2,3], 'b', label='Line 2')
line_2, = plt.plot([3,2,1], 'g', label='Line 1')
line_3, = plt.plot([3,2,1], 'r', label='Line 1')
line_4, = plt.plot([3,2,1], 'gray', label='Line 1')
plt.legend([line_1, line_2, line_3, line_4], ['Advection_above_LAB', 'Advection_below_LAB', 'Diffusion', 'Total=Adv+Diff'], loc=2)
plt.show(block=False)
plt.savefig('analytic_buoyancy.png', dpi=1000)
#plt.savefig('analytic_buoyancy.ps', format='ps', dpi=1000)
print('Drho_comp_crit= %.1f [kg/m3]' %Drho_comp_crit)
print('timeSt1=%.2f My' %(time1/Mytosec))
input("Press Enter to continue...")
#plt.close()

