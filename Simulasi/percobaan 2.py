import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def k(T,Ea,Ao) :
    #T dalam Kelvin
    R = 8.314
    K = Ao*np.exp(-Ea/(R*T))
    return (K)

def r(C,t):
    #t dalam sekon
 
    Aox = np.array([5.57e12,1.34e21,4.77e13,5.08e32,1.08e32])
    Eax = np.array([175.4e3,250e3,190.9e3,387.7e3,377.2e3])
    Tx = 573
    #penentuan nilai k   
    k1 = k(Tx,Eax[0],Aox[0])
    k2 = k(Tx,Eax[1],Aox[1])
    k3 = k(Tx,Eax[2],Aox[2])
    k4 = k(Tx,Eax[3],Aox[3])
    k5 = k(Tx,Eax[4],Aox[4])
    #penentuan variabel konsentrasi
    CSA = C[0]
    CHEPD = C[1]
    COCTD = C[2]
    CPEND = C[3]
    CHEXD = C[4]
    COCTDL = C[5]
    
    dCSAdt = -k1*CSA
    dCHEPDdt = k2*COCTDL
    dCOCTDdt = k3*COCTDL
    dCPENDdt = k4*COCTDL
    dCHEXDdt = k5*COCTDL
    dCOCTDLdt = k1*CSA-(k2+k3+k4+k5)*COCTDL

    return (dCSAdt,dCHEPDdt,dCOCTDdt,dCPENDdt,dCHEXDdt,dCOCTDLdt)

#integrasikan
#buat array waktu
t = np.linspace(0,18000,18000)
Co = np.array([0.15,0,0,0,0,0])
C = odeint(r,Co,t)
ulang = 6
#plot
for i in range (ulang):
    plt.plot(t,C[:,i])

plt.xlabel('Time (s)')
plt.ylabel('Concentration (M)')
plt.legend(['CSA','CHEPD','COCTD','CPEND','CHEXD','COCTDL'])
plt.show()
