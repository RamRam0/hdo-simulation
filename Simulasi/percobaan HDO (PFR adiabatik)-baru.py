# if __name__ == '__main__':
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#konstanta
R = 8.314 # Konstanta gas ideal (J/mol K)
rhoLA = 880 #Massa jenis asam laurat (kg/m^3)

#Parameter kinetika
Aox = np.array([0.5887,8.696])#Nilai konstanta arrhenius (s^-1)
Eax = np.array([38887,49972]) #Nilai energi aktivasi (J/mol)

#variabel tetap
L = 3 #L = panjang reaktor (m)
Po = 3000000 #Tekanan awal (Pa)
Rasio = 150 #rasio h2/umpan (Nm3/m^3)

 #variabel bebas
To = 280+273 #Temperatur awal masuk reaktor (K)
tpfr = 3 #waktu tinggal di reaktor (jam)
RR = 2 #Rasio Recycle

#Data berat molekul (kg/kmol)
MrLA = 200.321
MrDD = 184.322
MrUD = 170.295
MrH2 = 2
MrH2O = 18.015
MrCO2 = 44.01

def k(T,Ea,Ao) :
    #T dalam Kelvin
    K = Ao*np.exp(-Ea/(R*T))
    return (K)

def r(C,z):

    #penentuan variabel integrasi
    CLA = C[0]
    CDD = C[1]
    CUD = C[2]
    CH2 = C[3]
    CH2O = C[4]
    CCO2 = C[5]
    T = C[6]
    
    #perhitungan nilai k (s^-1)
    k1 = k(T,Eax[0],Aox[0])
    k2 = k(T,Eax[1],Aox[1])
    
    #perhitungan Cp #J/mol K
    CpLA = (-4.295)+(1.2373*T)+((-8.2209e-4)*T**2)\
            +((2.768e-7)*T**3)+((-3.8871e-11)*T**4)
    CpDD = (71.498)+(0.72559*T)+((1.1553e-4)*T**2)\
            +((-4.12e-7)*T**3)+((1.4141e-10)*T**4)
    CpUD = (125.21)+(0.31401*T)+((7.9137e-4)*T**2)\
            +((-9.141e-7)*T**3)+((2.7568e-10)*T**4)
    CpH2 = ((3.249)+((0.422e-3)*T)+((0.083e5)*T**(-2)))*R
    CpH2O = ((3.47)+((1.45e-3)*T)+((0.121e5)*T**(-2)))*R
    CpCO2 = ((5.457)+((1.045e-3)*T)+((-1.157e5)*T**(-2)))*R
    
    #perhitungan dHf #J/mol
    dHfLA = ((-582.24)+(-0.23113*T)+((1.2546e-4)*T**2))*1000
    dHfDD = ((-225.66)+(-0.25979*T)+((1.3823e-4)*T**2))*1000
    dHfUD = ((-208.56)+(-0.24686*T)+((1.3203e-4)*T**2))*1000
    dHfH2 = 0+((3.249*(T-298))+(((0.422e-3)/2)*(T**2-298**2))\
            +(-(0.083e5)*((1/T)-(1/298)))*R)
    dHfH2O = -241818+((3.47*(T-298))+(((1.45e-3)/2)*(T**2-298**2))\
            +(-(0.121e5)*((1/T)-(1/298)))*R)
    dHfCO2 = -393509+((5.457*(T-298))+(((1.045e-3)/2)*(T**2-298**2))\
            +((1.157e5)*((1/T)-(1/298)))*R)

    #perhitungan variabel beta dan laju alir linear
    beta = T/To
    u = L/(tpfr*3600) #laju alir linear (m/s)
    
    #persamaan diferensial neraca massa
    dCLAdz = -(k1+k2)*CLA/(beta*u)
    dCDDdz = k1*CLA/(beta*u)
    dCUDdz = k2*CLA/(beta*u)
    dCH2dz = -(3*k1)*CLA/(beta*u)
    dCH2Odz = (2*k1)*CLA/(beta*u)
    dCCO2dz = k2*CLA/(beta*u)

    Ctot = CLA+CDD+CUD+CH2+CH2O+CCO2 #mol/L
    rho = ((CLA*MrLA)+(CDD*MrDD)+(CUD*MrUD)+(CH2*MrH2)\
            +(CH2O*MrH2O)+(CCO2*MrCO2)) #kg/m^3
    Cp = ((CpLA*CLA/MrLA)+(CpDD*CDD/MrDD)+(CpUD*CUD/MrUD)+(CpH2*CH2/MrH2)\
          +(CpH2O*CH2O/MrH2O)+(CpCO2*CCO2/MrCO2))*1000/Ctot #J/kg K
    dHHDO = ((dHfDD+(2*dHfH2O))-((3*dHfH2)+dHfLA)) #J/mol
    dHDCO2 = ((dHfDD+dHfCO2)-(dHfLA)) #J/mol

    #Persamaan diferensial neraca energi
    dTdz = ((-dHHDO*k1*CLA)+(-dHDCO2*k2*CLA))*1000/(rho*Cp*u*beta)

    return (dCLAdz,dCDDdz,dCUDdz,dCH2dz,dCH2Odz,dCCO2dz,dTdz)


#Konsentrasi umpan
CLAo = ((1*rhoLA/MrLA)*((Po/(R*To))/1000))/((Rasio*101325/\
        (1000*R*293))+(1*rhoLA/MrLA))#mol/L
CH2o = (Rasio*101325/(1000*R*293))*((Po/(R*To))/1000)/\
       ((Rasio*101325/(1000*R*293))+(1*rhoLA/MrLA))#mol/L

#konstanta parameter kinetika masuk reaktor
k1x = k(To,Eax[0],Aox[0])
k2x = k(To,Eax[1],Aox[1])

#Konsentrasi umpan masuk reaktor
#dengan asumsi konversi asam laurat di reaktor adalah 100%
CLAox = CLAo*(1/(RR+1))
CDDox = ((k1x*CLAo)/(k1x+k2x))*(RR/(RR+1))
CUDox = ((k2x*CLAo)/(k1x+k2x))*(RR/(RR+1))
CH2ox = CH2o
CH2Oox = 0
CCO2ox = 0

#syntax untuk integrasi
#buat array rentang integrasi
zspan = 100 #jumlah rentang integrasi
z = np.linspace(0,L,zspan)
#konsentrasi awal masuk reaktor
Co = np.array([CLAox,CDDox,CUDox,CH2ox,CH2Oox,CCO2ox,To])
#Integrasi di dalam reaktor (konsentrasi dan temperatur)
C = odeint(r,Co,z)

#plot temperatur di sepanjang reaktor
plt.plot(z,C[:,6])
plt.xlabel('z (m)')
plt.ylabel('Temperatur (K)')
plt.title('Profil Temperatur di Sepanjang Reaktor')
plt.show()

#Plot konsentrasi asam laurat, dodekana dan undekana di sepanjang reaktor
ulang = 3
for i in range (ulang):
    plt.plot(z,10**3*(C[:,i]))

plt.xlabel('z(m)')
plt.ylabel('Konsentrasi (x 10^(-3) M)')
plt.legend(['C Asam laurat','C Dodekana','C Undekana'])
plt.title('Profil Konsentrasi di Sepanjang Reaktor')
plt.show()

#Plot Konversi Asam Laurat di sepanjang reaktor
KonvLA = np.zeros(zspan-1)
for i in range (zspan-1):
    KonvLA[i] = -(C[i+1,0]-C[0,0])*100/C[0,0]
    
plt.plot(z[0:(zspan-1)],KonvLA)
plt.xlabel('z (m)')
plt.ylabel('% Konversi Asam Laurat')
plt.title('Profil % Konversi Asam Laurat di Sepanjang Reaktor')
plt.show()

#Plot selektivitas HDO di sepanjang reaktor
SHDO = np.zeros(zspan-1)
for i in range (zspan-1):
    SHDO[i] = -(C[i+1,1]-C[0,1])*100/(C[i+1,0]-C[0,0])
    
plt.plot(z[0:(zspan-1)],SHDO)
plt.xlabel('z (m)')
plt.ylabel('% Selektivitas HDO')
plt.title('Profil % Selektivitas HDO di Sepanjang Reaktor')
plt.show()
