import matplotlib.pyplot as plt
import numpy as np
import param as prm
import read_cesm_values as cesm
#-------------------------------------------------------------------------

# define constants
mma = 28.96   # molar mass of air (g/mole)
mmo = 48      # molar mass of ozone in (g/mole)
mmcl = 35.45  # molar mass of Cl (g/mole)
mmcn = 97.46  # molar mass of chlorine natrate ClONO2 (g/mole)
mms = 96.06   # molar mass of sulfate (g/mole)
an = 6.023e23 # Avogadros number of (molec/mole)
re = (6371)   # radius of earth (km)
pi = 3.14     # circum/diameter 
sin60 = 0.866
da20 = 0.088   # density of air at 20 km altitude, Kg/m^3
das = 1.293    # density of air at standard temp and pressure (STP) (surface air) in Kg/m^3
kb = 1.38e-23  # Boltzmann Constant
du = 2.6867e20  # molec/m^2  in one Dobson Unit
rad = 0.05     # approx radius (um) of Antacrtic spring NAT aerosol, Fiebig et al (2014)
ndN = 5e-4     # Antarctic (surface measured) NAT (number density cm^-3), Sept-Oct, Weimer et al, (2022) 
dna = 1.5      # approx density of nitric acid HNO3 g/cm^3
dH2O = 1.0     # density of H2O g/cm^3
gms = mms/an   # molecular mass of sulfate (g/molec)
#--------------------------------------------------------------------------
t0 = 1
# remaining Cl (as ratio) in atm at time (t) yrs in future from 2023
clrt = (prm.clp_noaa + (prm.dcl * t0))/1e12
# -------------------------------------------------------------------------

# calculate air density molecules/m^3
dac = das/da20   # air density conversion factor from 20 km tao surface air equivalent
gma = mma/an     # g/molec air
gmo = mmo/an     # g/molec ozone
kma = gma/1000   # kg/molec air
kmo = (gmo/1000) # kg/molec ozone
ma20 = da20/kma  # molec air / m^3 at 20 km
# -------------------------------------------------------------------------

# calculate density of Cl in g/cm^3
gmcl = mmcl/an   # g/molec Cl
mc20 = ma20 * prm.clr  # Cl molec/m^3 at 20 km
mc20c = mc20*10**-6     # Cl molec/cm^3 at 20km
dmc20 = ma20 * prm.dclr  # annual change in Cl molec/m^3 at 20km
dmc20c = (dmc20)*10**-6  # annual change in Cl molec/cm^3 at 20km
c20 = mc20c * gmcl       # Cl (g/cm^3) at 20 km
dc20 = dmc20c * gmcl     # annual change in Cl (g/cm^3) at 20km
print(c20)
c = c20
dc = dc20

#  clm = (mc20 * gmcl)  # Cl in g/m^3
#  clc = clm * 1e-6  # Cl in g/cm^3
# dclc = (dmc20c * gmcl)  # annual change in Cl in g/cm^3 
# ------------------------------------------------------------------------

# calculate volume of Antarctic mid stratosphere, from 18-25 km altitude
r25 = re+(25)  # radius of earth + 25 km altitude
r18 = re+(18)  # radius of earth + 18 km altitude
h25 = r25-(r25*sin60)
h18 = r18-(r18*sin60)
vlc = (2/3)*(pi*((r25)**2)*h25)  # volume of larger cone (60-90s), Re + 25 km
vsc = (2/3)*(pi*((r18)**2)*h18)  # volume of smaller cone (60-90s), Re + 18 km 
vams = vlc - vsc      # volume of Antarctic mid strat 18-25 km altitude (60-90s) (km^3)
vamsm = (vams)*10**9 # volume of Antarctic mid strat 18-25 km altitude (60-90s) (m^3)

#  calculate volume of global mid-stratosphere (18-25 km):
#  vol25 = (4.0/3.0)*pi*r25**3.0  # volume of sphere of radius, earth + 25 km altitude
#  vol18 = (4.0/3.0)*pi*r18**3.0  # volume of sphere of radius, earth + 18 km altitude
#  vmsk = vol25 - vol18  # volume of mid stratosphere (km^3) from 18 to 25 km altitude
#  vmsm = (vmsk)*10**9   # volume of mid stratosphere in m^3
# -------------------------------------------------------------------------

# calculate baseline aerosol mass density in air:
vNAT = (4/3)*pi*rad**3  # volume of representative NAT particle (um^3)
vdN = vNAT*ndN*1e-12   # volume density of NAT (cm^3/cm^3 air)
dNAT = ((dna + (dH2O*3))/4)  # liquid density of NAT (g/cm^3)
dNTa = vdN * dNAT      # density of NAT in air (g NAT/cm^3 air)
dNT20 = ((dNTa)*1e-3)/((da20)*1e-6)  #density of NAT in air (cm^3 NAT/cm^3 air) at 20 km altitude
blae = dNTa            # baseline aerosol (density of NAT in air)
#-------------------------------------------------------------------------

# number density of baseline aerosol (as NAT) in mid stratosphere:
sig = 1.6  # geometric standard deviation for Aitkin and Accumulation mode, from CESM2
rib = 0.6e-6  # avg particle radius (m) for size distribution (bin?)
xib = dNT20  # mixing ratio (kg NAT/kg air) at 20 km
ros = 1770  # density of sulfate kg/m^3
vib = (4/3)*(pi*(rib**3))*np.exp((9/2)*(np.log(sig))**2)
nib = xib/(ros*vib)  # number density of sulfate 
bNAT = 5e-4  # from Weimer et al (2023)  

tts = np.linspace(1,prm.year_e-prm.year_b+1,prm.year_e-prm.year_b+1)
sty = tts+prm.year_b-1  # potential start year
tse = np.array([0,2,7,12,17,22,27,32,37,42])  # yrs from 2023, five yr steps after 2025, extended to 2065

#-------------------------------------------------------------------------

calculate baseline O3 depletion without SAI: 
# def odb(x, y = bNAT):
# return (((c+(dc*x))/gmcl)*y)*k*2*dac
# a = tts
# for wt in (tts):
# dObt = odb(a)
#-------------------------------------------------------------------------

# calculate mass density of sulfate aerosol to add to Antarctic mid stratosphere for 1 K surface cooling, in g/cm^3 :

# dela = 2.0e9      # total mass (kg) of additional aerosol needed per yr for 1K cooling, 2 Tg (until 2045), from CESM2-WACCM
# dae = dela/vamsm  # mass density of added aerosol in Antarctic mid stratosphere, kg/m^3
# daec = dae*1e3*1e-6  # mass density of added aerosol in Antarctic mid stratosphere, g/cm^3

dela = 20  # mg/m^2 total column SO4 aerosol increase from CESM2-WACCM, fig. 2 Richter et al (2022), 60-90S, 2035-2054
daec = (dela/7000)*1e-3 * 1e-6  # mass density g/cm^3 of additional SO4 aerosol needed in Antarctic mid stratosphere 18-25km 
#--------------------------------------------------------------------------

gms = mms/an  # molecular mass of sulfate (g/molec)
mds = daec/gms  # molecular density of added sulfate (molec/cm^3) at 20km (aproxmated from 18-25km altitude)

# calculate mass density of SO4 from CESM values (ARISE), take year 2 as test
# dae_A2 = prm.dae_A[2]*da20 
# dae_A2 = dae_A2 * 1e-6 * 1e-6
# -------------------------------------------------------------------------

# calculate number density of sulfate aerosol to add to mid stratosphere for 1 K surface cooling:
sig = 1.6    # geometric standard deviation for Aitkin and Accumulation mode, from CESM2
ri = 0.6e-6  # avg particle radius (m) for size distribution (bin?)
xi = 0.2e-7  # mixing ratio (kg SO4/kg air)
vi = (4/3)*(pi*(ri**3))*np.exp((9/2)*(np.log(sig))**2)
ni = xi/(ros*vi)  # number density of additional mid-stratospheric sulfate needed for 1 K surface cooling 
# -------------------------------------------------------------------------

# calculate heterogeneous reaction rate for R1:
# k1 = 1000  # reaction rate for (R1); (rxns cm^-3 s^-1) at 200K in mid stratosphere, Borrmann et al, (1997)
ut = 0.2     # uptake coefficient (gamma)for ClONO2 + HCl on H2SO4/H2O (binary aerosol) at 198K, Peter, (1997)
te1 = 198    # initial temperature (K) (PSC formation)
mcn = (mmcn/an)/1000  # molecular mass of ClONO2 (kg/molec)
m = mcn
cgas = ((8*kb*te1)/(pi*m))**1/2  # mean molecular speed of gas phase ClONO2 (m/s)
cgsc = cgas*100              # mean molecular speed of gas phase ClONO2 in (cm/s)
print(cgsc)
sad = 8.6  # surface area density um^2 cm^-3 from Tilmes et al (2022) CESM2 data, multi-year average following initial 3 yr particle growth phase 
sadl = np.array([2.0,2.0,3.0,4.0,5.0,8.0,9.0,8.0,8.5,8.0,7.5,10.0,10.0,9.5,8.0,7.0,8.5,9.0,9.0,8.0,8.5,8.0,10.0]) # SAD data (SAI), Tilmes et al,(2022) for CESM2
sade = np.array([8.60,6.45,4.30,2.15,17.20,25.80,34.40,43.00,68.80,86.00])  # SAD error scenario values (CESM2 post growth avg x 1, .75, .5, .25, 2,3,4,5,8,10) um^2/cm^3

k  = 0.25*ut*cgas*sad  # Het reaction rate (rxns/cm^3 s) for (R1): ClONO2 + HCl --> Cl2 + HNO3; calculated using Wegner et al,(2012) Eq
kc = 0.25*ut*cgsc*sad  # Het reaction rate (rxns/cm^3 s) for (R1): ClONO2 + HCl --> Cl2 + HNO3; calculated using Wegner et al,(2012) Eq
print(k)

def ksal(x,y = cgsc):
    return (0.25*ut*y*x)
a = sadl
for sad in (sadl):
    ksl = ksal(a) # gives reaction rate (rxns cm^-3 s^-1) for each year into SAI deployment, based on successive temporal SAD values following initial deployment
ksl = np.array([ksl])

cgas = ((8*kb*prm.temp)/(pi*m))**1/2  # mean molecular speed of gas phase ClONO2 (m/s) over temp range
print(cgas)
#k = 0.25*ut*cgas*sad  # Het reaction rate (rxns/cm^3 s) for (R1): ClONO2 + HCl --> Cl2 + HNO3; calculated using Wegner et al,(2012) Eq

# sades = np.array([4.3,17.2,25.8,43.0,86.0])  # list of 'error margin scenario' SAD values (0.5,2,3,5,10 x SAD value from CESM2)  
sadh = 4.3  # half of CESM2 projected SAD
sad2 = 17.2  # 2 x CESM2 projected SAD
sad3 = 25.8
sad5 = 43.0
sad10 = 86.0
# Note: SAD values, for a 2020 start, (and for each year into SAI deployment) are applied here to later start dates, 
# assuming conditions (such as temperature) have not changed significantly since 2020 to affect SAD values.

k = np.zeros((np.size(prm.sadl),np.size(prm.gamma),np.size(cgas)), dtype = float)

for j,gamma_j in enumerate(prm.gamma):
    for i,sad_i in enumerate(prm.sadl):
        for l,cgas_l in enumerate(cgas):
            k[i,j,l] =   0.25*gamma_j*cgas_l*sad_i 
print(k[0])    
print('')     
# ---------------------------------------------------------------------------

# # calculate ozone depletion resulting from addition of new aerosol in SAI 1K cooling scenario, variable start times

dO23 = -(((c+(dc * t0))/gmcl) * ni)*k*2   # ozone depletion from (R1), 2023 SAI start (tts= 0) for 1K surface cooling, g cm^-3 s^-1
#dO35 = -(((c+(dc * t3))/gmcl) * ni)*k*2  # ozone depletion from (R1), 2035 SAI start (tts=12) for 1K surface cooling, g cm^-3 s^-1
#dO45 = -(((c+(dc * t5))/gmcl) * ni)*k*2  # ozone depletion from (R1), 2045 SAI start (tts=22) for 1K surface cooling, g cm^-3 s^-1

d23m = dO23   # ozone depletion from reaction R1 (2023 SAI start) in molec/cm^3 
#d35m = dO35  # ozone depletion from reaction R1 (2035 SAI start) in molec/cm^3 
#d45m = dO45  # ozone depletion from reaction R1 (2045 SAI start) in molec/cm^3 

d23s = d23m*dac  # ozone depletion (surface air equivalent) from R1 (2023 start) molec/cm^3, at STP
d23sr = np.around(d23s, 2)  # rounded to 2 decimal places
#d35s = d35m*dac  # ozone depletion (surface air equivalent) from R1 (2035 start) molec/cm^3, at STP
#d35sr = np.around(d35s, 2)
#d45s = d45m*dac  # ozone depletion (surface air equivalent) from R1 (2045 start) molec/cm^3, at STP
#d45sr = np.around(d45s, 2)

# d23l = -(((c+(dc * t0))/gmcl) * ni)*ksl*2*dac  # ozone depletion (surface air eq) from R1 (2023 start) molec/cm^3, at STP, variable SAD & k
# d25l = -(((c+(dc * t1))/gmcl) * ni)*ksl*2*dac  # ozone depletion (surface air eq) from R1 (2025 start) molec/cm^3, at STP, variable SAD & k
# d30l = -(((c+(dc * t2))/gmcl) * ni)*ksl*2*dac  # ozone depletion (surface air eq) from R1 (2030 start) molec/cm^3, at STP, variable SAD & k
# d35l = -(((c+(dc * t3))/gmcl) * ni)*ksl*2*dac  # ozone depletion (surface air eq) from R1 (2035 start) molec/cm^3, at STP, variable SAD & k
# d40l = -(((c+(dc * t4))/gmcl) * ni)*ksl*2*dac  # ozone depletion (surface air eq) from R1 (2040 start) molec/cm^3, at STP, variable SAD & k
# d45l = -(((c+(dc * t5))/gmcl) * ni)*ksl*2*dac  # ozone depletion (surface air eq) from R1 (2045 start) molec/cm^3, at STP, variable SAD & k
# --------------------------------------------------------------------------------------------------------

#  calculate stratospheric column ozone depletion in Dobson Units (DU)
d23du = (d23m*7000)/du  # O3 depletion from R1 (2023 start) in Dobson Units (using mid-strat column from 18-25 km)
d23rdu = np.around(d23du, 2)
#d35du = (d35m*7000)/du  # O3 depletion from R1 (2035 start) in Dobson Units   "       "        "           "
#d35rdu = np.around(d35du, 2)
#d45du = (d45m*7000)/du  # O3 depletion from R1 (2045 start) in Dobson Units   "       "        "           "
#d45rdu = np.around(d45du, 2)
# --------------------------------------------------------------------------------------------------------

# list of dO3 for 1 yr step successive future SAI start times (molec m^-3)

def odp(x, y = ni):
   return (((c+(dc*x))/gmcl)*y)*k*2*dac
a = tts
for wt in tts:
    dO3t = odp(a)

# print outputs with string ('words.......') labels

print('The value of k, (rxns/cm^3 s) calculated using Eq from Wegner et al, (2012) is:')
print(k[0])
print('')
print('O3 depletion in molec/cm^3 s (surface air equivelent) from (R1) for a 2023 SAI start, 1K surface cooling, is:')
print(d23sr)
print('')
print('O3 depletion, molec/cm^3 s (surface air equivaltent) from (R1) for a 2035 SAI start, 1K surface cooling, is:')
print(d35sr)
print('')
print('O3 depletion, molec/cm^3 s (surface air equivaltent) from (R1) for a 2045 SAI start, 1K surface cooling, is:')
print(d45sr)
print('')
print('The O3 depletion in DU, from (R1) for a 2023 SAI start, 1K surface cooling, is:')
print(d23rdu)
print('')
print('The O3 depletion in DU, from (R1) for a 2035 SAI start, 1K surface cooling, is:')
print(d35rdu)
print('')
print('The O3 depletion in DU, from (R1) for a 2045 SAI start, 1K surface cooling, is:')
print(d45rdu)

# -------------------------------------------------------------------------
 
#  plot O3 depletion results by yearly potential SAI start dates beginning with 2023 

fig=plt.figure(1,figsize=[24,12])
plt.plot(sty,dO3t,'o', markersize=10, label='1K cooling')
plt.plot(2025,cesm.clox_oct_150[0],'x', markersize=20) #baseline ClOx destruction
plt.plot(2040,cesm.clox_oct_150[4],'d', markersize=20) #baseline ClOx destruction
plt.xlabel('start date', fontsize=14)
plt.ylabel('dO3 (-molec cm^-3 s^-1)', fontsize=14)    
plt.title('SAI O3 depletion (Cl dependent) by start date', fontsize=18)
plt.savefig('ClOx_code.png')
plt.show()

"""
plt.plot()
#-------------------------------------------------------------------------

# Calculate O3 depletion vs Cl with SAD error scanarios:
def odps(x,y = ni):
    def ksae(w,z = cgsc):
        return (0.25*ut*z*w)   
    a = sade
    for sad in sade:
        kse = ksae(a) # reaction rate (rxns cm^-3 s^-1) based on error scenario values for SAD
    kse = np.array([kse]) 
    
    return (((c+(dc*x))/gmcl)*y)*kse*2*dac   
a = tts   
for ts in tts:
    docs = odps(a)    
docs = np.array([docs])   
#-------------------------------------------------------------------------

# temperature profile data from Wang and Lin (2007) fig 2b avg btw GPS and Radiosonde data 9/15/2006 Neumayer Sta Antarctica approx 70 S Lat
alt = np.array([2,4,6,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,27])  # altitude (km)
hPa = np.array([780,575,400,300,225,190,155,130,110,95,80,68,50,45,40,32,28,22,18,17,12])  # pressure altitude (hPa)
tpt = np.array([256,245,226,211,200,192,186,188,187,187,186,187,188,188,189,190,193,197,200,203,213]) # temperature (K) 
#------------------------------------------------------------------------  

# plot temperature profile
#------------------------------------------------------------------------
plt.close()
fig=plt.figure(figsize=[12,6])
plt.plot(tpt,alt,'o', markersize=1, label=('Neumayer Station 9/15/2006')
plt.xlabel('temperature (K)', fontsize=14)
plt.ylabel('altitude (km)', fontsize=14)    
plt.title('Antarctic spring 2006 temperature profile', fontsize=18)
plt.show()
plt.plot()
#-----------------------------------------------------------------------

"""


    

