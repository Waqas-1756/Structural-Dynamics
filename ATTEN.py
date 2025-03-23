############ Importing Libraries #################

import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt

####################### Input Data ############################

##############Campbell and Bozorgnia (1994) ####################


Mw=np.array([5,6,7,8])
R=np.linspace(1,100,100)
## Campbell
N=Mw.shape[0]
O=R.shape[0]
#l_R=np.log(R)
F=0   ## Normal Faulting 
Shr=1  ## for hard rock 
Ssr=0
PGA=np.zeros((N,O))
for i in range(N):
    for j in range(O):
        PGA[i,j]=-3.512+0.904*Mw[i]-1.328*np.log((np.sqrt((R[j]**2+(0.149*np.exp((0.647*Mw[i])))**2))))
        +(1.125-0.112*np.log(R[j])-0.0957*Mw[i])*F+(0.440-0.171*np.log(R[j]))*Ssr+(0.405-0.222*np.log(R[j]))*Shr
for i in range(N):
    plt.plot(np.flip(R),abs(PGA[i])/9.81,label=f'Mw={Mw[i]}')
plt.xlabel('R (Kms)')
plt.ylabel('ln(PGA) (g)')
plt.xlim(3,100)
plt.legend()


###################### Atkinson and Boore (1993)###################

for i in range(N):
    for j in range(O):
        PGA[i,j]=-1.841+0.686*(Mw[i]-6)-0.123*(Mw[i]-6)**2-np.log(R[j])-0.00311*R[j]
for i in range(N):
    plt.plot(np.flip(R),abs(PGA[i])/9.81,label=f'Mw={Mw[i]}')
plt.xlabel('R (Kms)')
plt.ylabel('ln(PGA) (g)')
plt.legend()

############################ Toro et al. (1994) ###############

for i in range(N):
    for j in range(O):
        Rm=np.sqrt(R[j]**2+9.3**2)
        PGA[i,j]=2.20+0.81*(Mw[i]-6)-1.27*np.log(Rm)+0.11*max(np.log(Rm/100),0)-0.0021*Rm
for i in range(N):
    plt.plot(np.flip(R),abs(PGA[i])/9.81,label=f'Mw={Mw[i]}')
plt.xlabel('R (Kms)')
plt.ylabel('ln(PGA) (g)')
plt.legend()

############################ Youngs et al. (1988)  ###############

Z=1 # For intraslab events
for i in range(N):
    for j in range(O):
        PGA[i,j]=19.16+1.045*(Mw[i])-4.738*np.log(R[j]+205.5*np.exp(0.0968*Mw[i]))+0.54*Z
for i in range(N):
    plt.plot(np.flip(R),abs(PGA[i])/9.81,label=f'Mw={Mw[i]}')
plt.xlabel('R (Kms)')
plt.ylabel('ln(PGA) (g)')
plt.legend()




