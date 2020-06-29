# Modelling Enslin and GopalKrishna 2001
# Major updates 

from getpar import getparspace
from astropy.io import ascii
from operations import chicalc
import numpy as np
import matplotlib.pyplot as plt


####################################	

def readfile(myspec):
	dat = ascii.read(myspec)
	return dat 
myfile ='A1914_spec.dat'   #col1 =freq col2 = flux col3 = error
spec = readfile(myfile)
print(spec)
f_nuobs=np.array(spec['nughz'])
f_err=np.array(spec['fmyerr'])
flx_obs = np.array(spec['fmy'])
flx_obs	=[i*(1e-26) for i in flx_obs]
#f_nuobs=[0.1]
#flx_obs = [530]
#f_err=[0]



#####################################	Input Parameters	#####################################
# Extract the redshift and properties of the system from source.py
import source
#z	=0	#t1 source.z
#B_src	=2.7	#t1 source.B
#V_src	=0.12	#t1 source.V

z	=0.17120		#t1 source.z
B_src	=5*1e-6	#t1 source.B
V_src	=0.037*(3.08*1e24)**3		#t1 source.V


# To Take User Defined inputs regarding assumed scenario and phase for the system.
# Also checks for illegal entry and assigns default entry. 
scen=input("Enter the scenario[A, B or C](default: scenario B): ")
if scen not in ['A','B','C']:
	scen='B'


phase=input("Enter the phase[0,1,2,3 or 4](default: phase 3): ")
if phase not in ['0','1','2','3','4']:
	phase='3'
phase=int(phase)
print('Choice of Scenario is:',scen,'	','Choice of Phase is',phase)

# Compression Ratio Index
b=[1.8,1.2,0,2.0,0]



#####################################	Operational Parameters	#####################################
# To generate time-scale parameter space and save in parex.dat
getparspace(scen,phase)


# Reading data from table to phase-wise iterative solutions 
timeex= ascii.read('parex.dat')

print(timeex)
# DEFINE PARAMETERS TO STORE RESULT OF EACH SET OF TIMESCALES
chilist=[]
flx=[]
count=0	#counter
result	=open('myresult.dat','w')
for i in timeex:

	# instancing time scale for each set
	delt	=[0.0,0.0,0.0,0.0,0.0]	#delt0=0 for all sets
	tau	=[0.0,0.0,np.inf,0.0,np.inf]	#Time scale for tau_4 and tau_2 is prescribed infinity
	#setting up tau
	tau[0]	=i['tau0']
	tau[1]	=(2.0/3.0)*tau[0]
	tau[3]	=i['tau3']
	#setting up delt
	delt[0]	=0
	delt[1]	=i['delt1']
	delt[2]	=i['delt2']
	delt[3]	=i['delt3']
	delt[4]	=i['delt4']	
	
	delt	=[i*(3.154*1e16) for i in delt]
	tau	=[i*(3.154*1e16) for i in tau]
	
	print('Iteration number',count+1)
	chi,flux=chicalc(delt,tau,phase,z,B_src,V_src,b,f_nuobs,f_err,flx_obs)
	chilist.append(chi)
	flx.append(flux)
	myresult	=str(count)+'\t'+str(chi)+'\n'
	result.write(myresult)
	count=count+1
chi_min_ind	=np.nanargmin(chilist)
chi_min		=chilist[chi_min_ind]
flx_min		=flx[chi_min_ind]

# instancing time scale for each set
delt_m	=[0.0,0.0,0.0,0.0,0.0]	#delt0=0 for all sets
tau_m	=[0.0,0.0,np.inf,0.0,np.inf]	#Time scale for tau_4 and tau_2 is prescribed infinity
#setting up tau
tau_m[0]	=timeex['tau0'][chi_min_ind]
tau_m[1]	=(2.0/3.0)*tau_m[0]
tau_m[3]	=timeex['tau3'][chi_min_ind]
#setting up delt
delt_m[0]	=0
delt_m[1]	=timeex['delt1'][chi_min_ind]
delt_m[2]	=timeex['delt2'][chi_min_ind]
delt_m[3]	=timeex['delt3'][chi_min_ind]
delt_m[4]	=timeex['delt4'][chi_min_ind]	
print('Chi_min',chi_min)
print('tau_min',tau_m)
print('delt_m',delt_m)
fig,ax = plt.subplots()
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Flux density (mJy)')
flx_obs		=[i/(1e-26) for i in flx_obs]
ax.plot(f_nuobs,flx_obs,'bo',label='obs')
flx_min		=[i/(1e-26) for i in flx_min]
ax.plot(f_nuobs,flx_min,'r-')
ax.legend()
plt.show()





	
###	
###	 while writing the program for chi square calculation for a time step,
###	 we send only one value of b that coresponds to user input phase
###	c=c+1
















 
