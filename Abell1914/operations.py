from scipy.special import gamma
#from scipy.stats import chisquare
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import quad
from astropy.io import ascii
import math
import matplotlib.pyplot as plt
import time

#Reading OUtput file:


# Global constants[Only includes constants which are used multiple places]
alphae	=2.50
e	=4.803*1e-10
m_e	=9.109*1e-28
pi	=np.pi
c	=2.99*1e10
sigma	=5.670*1e-5




def integrand(p,nu,Bi,C0i,alphae,pstar0i):
# Input Parameter:
#nu	: Frequency of Observation
#p	: Integration variable
#Bi	: Magnetic Field at that phase
	const	=((2.00**(2.00/3.00))*((pi/3.00)**(3.00/2.00)))/gamma(11/6)
	nu_i	=(3.00*e*Bi*(p**2))/(4*pi*m_e*c)
	x	=nu/nu_i
	ftilde=1
	ans=ftilde*(const*(x**(1.00/3.00))*math.exp(-(11.00/8.00)*(x**(7.00/8.00))))*((C0i**((alphae+2)/3))*(p**(-alphae))*((1-(p/pstar0i))**(alphae-2)))
	return	ans

def formula(p_min,p_max,N,nu,Bi,C0i,alphae,pstar0i):
	h = (p_max-p_min)/N
	soln = h*(integrand(p_min,nu,Bi,C0i,alphae,pstar0i)+integrand(p_max,nu,Bi,C0i,alphae,pstar0i))/2.00
	for i in range(1,N):
		soln += h*integrand(p_min+(i*h),nu,Bi,C0i,alphae,pstar0i)
	return soln




def pstarji(cji,ti,uBj,uC,bi):
	a0	=0.325*1e-7
	if bi==0:
		ps=(1/a0)*(1/(uBj+uC)*(1/(3.154*1e16)))
	else:
		denm1	=a0*ti
		denm2	=(((cji**(-5*bi/3)-cji**(-1))/(1-(5*bi/3)))*uBj)+(((cji**(-bi/3)-cji**(-1))/(1-(bi/3)))*uC)
		ps=(cji**(-bi/3))/(denm1*denm2)
	return 	ps
	
		
def getdl(redshift):
#Input Parameter:
#redshift	:redshift of the source
#OUtput Paramter:
#Luminosity Distance	:Luminosity Distance of the source according to Lambda CDM model of the universe
	# Function to obtain Luminosity Distance of the source
	cosmo	=FlatLambdaCDM(H0=70,Om0=0.27)	
	dl	=cosmo.luminosity_distance(redshift)
	return dl.value	#result in Mpc
	


# Function to Calculate operational parameters at each phase
def phase_calc(ph,b,del_t,tau,B_src,V_src):
# Input Parameters:
#del_t	:Duration of phase		List	
#tau	:Timescale of expansion		List
#ph	:user input for phase		integer
#B_src	:Magnetic Field Of Source	float
#V_src	:Volume Of Source		float
#b	:compression index of phase	list

# Output Parameters:
#B_x	: Magnetic Field of source at each phase	List
#u_B_x	: Magnetic Energy Density at each phase		List
#V_x	: Volume of source at each phase		List
#C_x	: Compression Ratio ta each phase		List

# Defining Variables which we pass back to caller
	B_x	=[0,0,0,0,0]	# instantiating variable to store Magnetic field for all phase
	u_B_x	=[0,0,0,0,0]	# instantiating variable to store Magnetic field Energy Density for all phase
	V_x	=[0,0,0,0,0]	# instantiating variable to store Volume for all phase
	C_x	=[0,0,0,0,0]	# instantiating variable to store Compression Ratio for all phase
# Assigning the variables to be passed
	# Assigning Compresion Ratios for all phases
	C_x	=[((1+(del_t[i]/tau[i]))**(-b[i])) for i in range(0,5)]	#eq. 15
	# Assigning Volume for all phases
	V_x[ph]	=V_src	
	for i in range(0,ph):
		V_x[ph-i-1]	=V_x[ph-i]*C_x[ph-i]			#eq. 2

	# Assigning Magnetic Energy Density and Magnetic Field at each step
	B_x[ph]		=B_src
	u_B_x[ph]	=((B_src**2.00)/(8.00*pi))  
	u_B_x[0]	=u_B_x[ph]*(V_x[ph]/V_x[0])**(4.00/3.00)	#eq. 17
	for i in range(1,ph):
		u_B_x[i]=u_B_x[0]*((V_x[0]/V_x[i])**(4.00/3.00)) 	
	for i in range(0,ph):
		B_x[i]	=np.sqrt(8.00*pi*(u_B_x[i]))		#Units of [mu gauss]
	return B_x,u_B_x,V_x,C_x
	
#def pstarji(cji,ti,uBj,uC,b):
#def fp(C0i,alphae,pstar0i):

def chicalc(del_t,tau,ph,z,B_src,V_src,b,f_nuobs,f_err,flx_obs):
# Input Parameter #
#del_t	:Duration of phase		List	
#tau	:Timescale of expansion		List
#ph	:user input for phase		integer
#z	:redshift of the source		float
#B_src	:Magnetic Field Of Source	float
#V_src	:Volume Of Source		float
#b	:compression index of phase	list

# Output Parameters #
#Flux

	# Obtaining Luminosity Distance
	d	=getdl(z)*3.08*1e24	# will be used in calculation of Flux from Luminosity at each frequency observed 
	# Equivalent IC Magnetic Field and IC Energy Density
	u_c	=4.2*1e-13*(1+z)**4

	# Calculating energy density, Compression Index for each phase(x B_ix: i represents phase number and x represents previous)
	B,u_B,V,C=phase_calc(ph,b,del_t,tau,B_src,V_src)
	C0=[0,0,0,0,0]
	for i in range(0,ph+1):
		C0[i]=V[0]/V[i]
#******************************** print data for debugging ********************************
#	print('IC energy is',u_c/(1.6*1e-12))	
#	for i in range(0,len(B)):
#		print('magnetic field (uBi) in', i ,'phase is',B[i]/1e-6)
#		print('magnetic energy density in', i ,'phase is',u_B[i]/(1.6*1e-12))
#		print('volume (V) in', i ,'phase is',V[i]/(3.08*1e24)**3)
#		print('C0i in',i,'phase is',C[i])
#******************************************************************************************
#	time.sleep(10)

	
	ipstar0i=0
	p_min	=10
	p_max	=1e4
	# calculate pstar_0ph
	for i in range(1,ph+1):
		t=del_t[i]+tau[i]
		com=(1+(del_t[i]/tau[i]))
		mag	=u_B[i-1]
		ipstar0i+=(((C0[i-1])**(1.00/3.00))/pstarji(com,t,mag,u_c,b[i]))	#eq 14 & pstarji is eq 16.
	pstar0i=1/ipstar0i
	
#	
###	# Calculating Luminosity and Flux received by the observer
	flx_cal	=[] 
	c3	=math.sqrt(3)*(e**3)/(4*pi*m_e*c**2)	#constant in equation 19	
	Bi	=B[ph]
	C0i	=C0[ph]
	Vi	=V[ph]
	for nu in f_nuobs:
			nu=nu*1e9
			L_freq	=c3*Bi*Vi*formula(10,pstar0i,1000,nu,Bi,C0i,alphae,pstar0i)
#			L_freq	=c3*Bi*Vi*quad(integrand,10,pstar0i,args=(nu,Bi,C0i,alphae,pstar0i))[0]
			flx_cal.append(L_freq/(4*pi*d**2))

	norm	=sum(flx_obs)/sum(flx_cal)
##	#	normalizing calculated result
	flx_norm	=[norm*i for i in flx_cal]
	chi	=0
	for i in range (0,len(f_nuobs)):
		chi=chi+((flx_obs[i]-flx_norm[i])**2/flx_norm[i])
#	print(flx_norm)
	return chi,flx_norm	
		
	
	
	






## Instruction to print parameters
#magnetic field intensity	:	divide each term with 1e-6		[conversion to mu Gauss]
#magnetic field energy density	:	divide each term with 1.6*1e-12		[conversion to eV]
#Volume				:	divide each term with (3.08*1e24)**3	[conversion to Mpc3]









