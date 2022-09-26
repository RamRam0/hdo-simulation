# if __name__ == '__main__':
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#konstanta
R = 8.314 # Konstanta gas ideal (J/mol K)
rhoLA = 880 #Massa jenis asam laurat (kg/m^3)
rhoDD = 750 #Massa jenis dodekana (kg/m^3)
rhoUD = 740 #Massa jenis undekana (kg/m^3)

#Parameter kinetika
Aox = np.array([0.0268,0.708])#Nilai konstanta arrhenius (s^-1)
Eax = np.array([24242,38077]) #Nilai energi aktivasi (J/mol)

#variabel tetap
L = 3 #L = panjang reaktor (m)
D = 1.6 #D = diameter dalam reaktor (m)
vf = 0.4 #void fraction reaktor
Po = 3000000 #Tekanan awal (Pa)
Rasio = 300 #rasio h2/umpan (Nm3/m^3)

#Untuk kondisi H2 umpan (Nm3)
Pi = 101325 #Pi = tekanan normal (Pa)
Ti = 273 #Ti = suhu normal (K)

# variabel bebas
tpfr = 3 #waktu tinggal di reaktor (jam)
wtLA = 0.5 #fraksi massa
wtDD = (1-wtLA)/2 #fraksi massa
wtUD = wtDD #fraksi massa

# Massa jenis campuran
rhoX = (wtLA * rhoLA) + (wtDD * rhoDD) + (wtUD * rhoUD)  # Massa jenis campuran (kg/m^3)

#Data berat molekul (kg/kmol)
MrLA = 200.321
MrDD = 184.322
MrUD = 170.295
MrH2 = 2
MrH2O = 18.015
MrCO2 = 44.01

#variasi suhu umpan masuk reaktor
Tvariasi = np.array([280+273,300+273,320+273,340+273])

#array rentang integrasi
zspan = 100 #jumlah rentang integrasi
z = np.linspace(0,L,zspan)

#array temperatur akhir
SHDO = np.zeros((zspan-1,4))

for i in range (4):
    # Konversi Rasio dari Nm3/m3 ke m3/m3
    To = Tvariasi[i]
    Rasioakhir = Rasio*To*(Pi/1000000)/(Ti*(Po/1000000)) # m3/m3

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
        CpLA = (50.801)+(2.258*T)+((-4.966e-3)*T**2)+((4.3771e-6)*T**3)
        CpDD = (84.485)+(2.0358*T)+((-5.0981e-3)*T**2)+((5.2186e-6)*T**3)
        CpUD = (94.169)+(1.7806*T)+((-4.6303e-3)*T**2)+((4.9675e-6)*T**3)
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


    # Konsentrasi umpan masuk reaktor
    CLAo = wtLA * rhoX * (1 / (Rasioakhir + 1)) / MrLA  # mol/L
    CDDo = wtDD * rhoX * (1 / (Rasioakhir + 1)) / MrDD  # mol/L
    CUDo = wtUD * rhoX * (1 / (Rasioakhir + 1)) / MrUD  # mol/L
    CH2o = Po / (1000 * R * To) * Rasioakhir / (Rasioakhir + 1)  # mol/L
    CH2Oo = 0
    CCO2o = 0

    #syntax untuk integrasi

    #konsentrasi awal masuk reaktor
    Co = np.array([CLAo,CDDo,CUDo,CH2o,CH2Oo,CCO2o,To])
    #Integrasi di dalam reaktor (konsentrasi dan temperatur)
    C = odeint(r,Co,z)

    #array untuk selektivitas HDO
    for j in range(zspan - 1):
        SHDO[j,i] = -(C[j + 1, 1] - C[0, 1]) * 100 / (C[j + 1, 0] - C[0, 0])

# Plot selektivitas HDO di sepanjang reaktor
for i in range (4):
    plt.plot(z[0:(zspan - 1)], SHDO[:,i])
plt.xlabel('z (m)')
plt.ylabel('% Selektivitas HDO')
plt.legend(['T = 280 C','T = 300 C','T = 320 C','T = 340 C'])
plt.title('Profil % Selektivitas HDO di Sepanjang Reaktor')
plt.show()