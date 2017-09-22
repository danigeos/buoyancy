#!/usr/bin/python
print('\n\n########################## Analytical slab buoyancy ############################')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os, sys, math, time

print('Number of arguments: %s' %(len(sys.argv)))	
print('Argument List: %s' %(str(sys.argv)))

g=9.8 #m s-2
alpha=1e-6 #m2 s-1  Diffusivity
Mytosec=1e6*365.25*24*3600


#-------------------------------------
#vel=3. #km/Myr
Drho_comp=25. #12, 33, 50 (archon), 64 kg m-3 Contrast between base lithosphere and top asthenosphere
vel=10 #mm/yr

drho_dT=.1175 #kg m-3 K-1
zm=25e3  #m
zl=200e3 #m
Tm=700  #C
Ta=1300 #C
sigma=45. #subduction angle in degrees

#time=5.*Mytosec #Myr
tmax=100
Dl=(zl-zm)/4
#-------------------------------------


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
Fsp_adv_temp=	np.linspace(0, 0, nit)
Fsp_adv_comp=	np.linspace(0, 0, nit)
Fsp_dif=		np.linspace(0, 0, nit)
Fsp_total=		np.linspace(0, 0, nit)


for i in range(0, nit):
	print('\n________________________________________________________')
	#vel_arr[i]*1e3/Mytosec
	time = time_arr[i]*Mytosec #Dl/vel
	Dl=time*vel
	time1 = (zl-zm)/math.sin(sigma)/vel
	timeL = time
	if (time>time1): 
		timeL=time1
		print('\nTime longer than Stage 1 (question results for diffusion)')
	print("i=%d  vel=%.2f km/My    time=%.2f My  timeSt1=%.2f My   Dl=%.1f km" %(i, vel*Mytosec/1e3, time/Mytosec, time1/Mytosec, Dl/1e3))
	Drho_comp_crit = 1./2 * drho_dT * (Ta-Tm)
	print('Drho_comp_crit= %.1f [kg/m3]' %Drho_comp_crit)
	dFsp_adv_temp_dt_ini = g * Drho_comp_crit * (zl-zm) * vel
	dFsp_adv_temp_dt = + dFsp_adv_temp_dt_ini * math.pow((1 - time/time1), 2)
	Drho_comp_avg = Drho_comp - 1./2 * drho_dT*(Ta-Tm)*math.sin(sigma)/(zl-zm)*vel*timeL
	print('Drho_comp_avg = %.2e  %.2e' %(Drho_comp_avg, 1./2 * drho_dT*(Ta-Tm) * math.sin(sigma)/(zl-zm)*vel*timeL))
	dFsp_adv_comp_dt = - Drho_comp_avg * g * (zl-zm) * vel
	print('Pull_force_per_My= %+.1e + %+.1e = %+.1e [N/m/My]' %((dFsp_adv_temp_dt)*Mytosec, (dFsp_adv_comp_dt)*Mytosec, (dFsp_adv_temp_dt+dFsp_adv_comp_dt)*Mytosec))
	#Advection Stage 1:
	Fsp_adv_temp[i] = dFsp_adv_temp_dt_ini * time1 / 3. * (1-math.pow(1-timeL/time1,3))
	Fsp_adv_comp[i] = - g*(zl-zm)*Drho_comp*vel*timeL + 1./4*g*drho_dT*(Ta-Tm)*math.sin(sigma) * math.pow(vel*timeL,2)
	#Advection Stage 2:
	if (time>time1): 
		Fsp_adv_comp[i] -= g * (Drho_comp - 1./2*drho_dT*(Ta-Tm)) * vel * (zl-zm) * (time-time1)

	print('\tPull by advection of temp. = %+.1e [N/m] (blue)'  %(Fsp_adv_temp[i]))
	print('\tPull by advection of comp. = %+.1e [N/m] (green)' %(Fsp_adv_comp[i]))



	#Diffusion lowermost triangle (Stage 1+2):
	ff = 4./3/np.pi**.5 * g * alpha**.5 * drho_dT * (Ta-Tm)/5 
	Fsp_dif[i] = - ff / 5 / (zl-zm) * math.sin(sigma) * math.pow(vel,2.) * ((-2*time-3*timeL) * math.pow(time-timeL,3./2) + 2*math.pow(time,5./2))
	#Diffusion above triangle (Stage 2):
	if (time>time1): 
		Fsp_dif[i] += - ff * math.pow(time-time1, 3./2) * vel


	#Diffusion below LAB - 2nd approach (simpler): 
	H=zl-zm
	bb=H*H/8./alpha
	ff = - 4. * g * drho_dT * (Ta-Tm) * bb * alpha / H * vel / time1
	expmint = math.exp(-time/bb)
	Fsp_dif[i] = ff * ( timeL*timeL/2. - bb*timeL - bb*bb*(expmint-1) ) 
	if (time>time1): 
		expmint1 = math.exp(-time1/bb)
		explust1 = math.exp(+time1/bb)
		Fsp_dif[i] += \
			- ff * bb*explust1 * (expmint*(time+bb)-expmint1*(time1+bb)) \
			+ ff * bb*bb*explust1*(expmint-expmint1) \
			+ ff * time1 * ( time+bb*explust1*expmint - time1 - bb ) 

	#Diffusion below LAB - 3rd approach (simplest and most accurate!): 
	ff = - 4. * alpha * g * (Ta-Tm) / H * drho_dT * bb * vel
	expmint = math.exp(-time/bb)
	Fsp_dif[i] = ff * ( time + bb*(expmint-1) ) 

	#Diffusion above LAB:
	#MISSING!!

	print('\tPull by diffusion of temp. = %+.1e [N/m] (red)' %Fsp_dif[i])


	Fsp_total[i]=Fsp_adv_temp[i]+Fsp_adv_comp[i]+Fsp_dif[i]
	print('\tTOTAL Slab Pull force      = %+.1e [N/m] (GREY)' %Fsp_total[i])

plt.cla()
plt.hlines(0, 0, tmax, 'k', 'dotted', '')
#plt.hlines(+1e13, 0, tmax, 'k', 'dotted', '')
#plt.hlines(-1e13, 0, tmax, 'k', 'dotted', '')
plt.plot(time_arr,Fsp_adv_temp, 'b', linewidth=1.5)
plt.plot(time_arr,Fsp_adv_comp, 'g', linewidth=1.5)
plt.plot(time_arr,Fsp_dif, 'r', linewidth=1.5)
plt.plot(time_arr,Fsp_total, 'grey', linewidth=3)
plt.ylim([-20e12,20e12])
#plt.xlabel("vel (km/My)")
plt.xlabel("time (My)")
plt.ylabel("Slab Pull (N/m)")
line_1, = plt.plot([1,2,3], 'b', label='Line 2')
line_2, = plt.plot([3,2,1], 'g', label='Line 1')
line_3, = plt.plot([3,2,1], 'r', label='Line 1')
line_4, = plt.plot([3,2,1], 'gray', label='Line 1')
plt.legend([line_1, line_2, line_3, line_4], ['Advection_above_LAB', 'Advection_below_LAB', 'Diffusion', 'Total=Adv+Diff'], loc=2)
plt.show(block=False)
plt.savefig('analytic_slab_pull.ps', format='ps', dpi=1000)
input("Press Enter to continue...")
#plt.close()
