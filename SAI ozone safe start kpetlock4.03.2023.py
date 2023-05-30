import matplotlib.pyplot as plt
import numpy as np
import param as prm
import read_cesm_values as cesm

# define initial chlorine (Cl) parameters

# remaining Cl (as ratio) in atm at time (t) yrs in future from 2023
clrt = (prm.clp_noaa + (prm.dcl * t0))/1e12
# -------------------------------------------------------------------------

# density of air

da20 = 0.088    # density of air at 20 km altitude, Kg/m^3
das = 1.293    # density of air at standard temp and pressure (STP) (surface air) in Kg/m^3
dac = das/da20  # air density conversion factor from 20 km tao surface air equivalent
# -------------------------------------------------------------------------

# calculate air density molecules/m^3

mma = 28.96      # molar mass of air in g/mole
mmo = 48         # molar mass of ozone in g/mole
an = 6.023e23    # Avogadros number of molec/mole
gma = mma/an     # g/molec air
gmo = mmo/an     # g/molec ozone
kma = gma/1000   # kg/molec air
kmo = (gmo/1000) # kg/molec ozone
ma20 = da20/kma  # molecules air / m^3 at 20 km
# -------------------------------------------------------------------------

# calculate density of Cl in g/cm^3

mc20 = ma20 * prm.clr  # Cl molec/m^3 at 20 km
mc20c = mc20*10**-6  # Cl molec/cm^3 at 20km
dmc20 = ma20 * prm.dclr  # annual change in Cl molec/m^3 at 20km
dmc20c = (dmc20)*10**-6  # annual change in Cl molec/cm^3 at 20km

mmcl = 35.45       # molar mass of Cl (g/mole)
gmcl = mmcl/an     # g/molec Cl
c20 = mc20c * gmcl # Cl (g/cm^3) at 20 km
dc20 = dmc20c * gmcl  # annual change in Cl (g/cm^3) at 20km
print(c20)
c = c20
dc = dc20

#  clm = (mc20 * gmcl)  # Cl in g/m^3
#  clgc = clm * 1e-6  # Cl in g/cm^3
# dclgc = (dmc20c * gmcl)  # annual change in Cl in g/cm^3 
# ------------------------------------------------------------------------

# calculate volume of Antarctic mid stratosphere, from 18-25 km altitude

re = (6371)  # radius of earth (km)
r25 = re+(25)  # radius of earth + 25 km altitude
r18 = re+(18)  # radius of earth + 18 km altitude
pi = 3.14
#  vol25 = (4.0/3.0)*pi*r25**3.0  # volume of sphere of radius, earth + 25 km altitude
#  vol18 = (4.0/3.0)*pi*r18**3.0  # volume of sphere of radius, earth + 18 km altitude
#  vmsk = vol25 - vol18  # volume of mid stratosphere (km^3) from 18 to 25 km altitude
#  vmsm = (vmsk)*10**9  # volume of mid stratosphere in m^3

sin60 = 0.866
h25 = r25-(r25*sin60)
h18 = r18-(r18*sin60)
vlc = (2/3)*(pi*((r25)**2)*h25)  # volume of larger cone (60-90s), Re + 25 km
vsc = (2/3)*(pi*((r18)**2)*h18)  # volume of smaller cone (60-90s), Re + 18 km 
vams = vlc - vsc # volume of Antarctic mid strat 18-25 km altitude (60-90s) (km^3)
vamsm = (vams)*10**9  # volume of Antarctic mid strat 18-25 km altitude (60-90s) (m^3)
# -------------------------------------------------------------------------

# calculate mass density of sulfate aerosol to add to mid stratosphere for 1 K surface cooling, in g/cm^3

# dela = 2.0e9 # total mass (kg) of additional aerosol needed per yr for 1K cooling, 2 Tg (until 2045), from CESM2-WACCM
# dae = dela/vamsm  # mass density of added aerosol in Antarctic mid stratosphere, kg/m^3
# daec = dae*1e3*1e-6  # mass density of added aerosol in Antarctic mid stratosphere, g/cm^3

dela = 20  # mg/m^2 total column SO4 aerosol increase from CESM2-WACCM, fig. 2 Richter et al (2022), 60-90S, 2035-2054
daec = (dela/7000)*1e-3 * 1e-6  # mass density g/cm^3 of added SO4 aerosol in Antarctic mid stratosphere 18-25km 

# calculate mass density from CESM values (ARISE), take year 2 as test

dae_A2 = prm.dae_A[2]*da20 
dae_A2 = dae_A2 * 1e-6 * 1e-6

# -------------------------------------------------------------------------

# calculate heterogeneous reaction rate for R1

# k1 = 1000  # reaction rate for (R1); (rxns cm^-3 s^-1) at 200K in mid stratosphere, Borrmann et al, (1997)
ut = 0.2  # uptake coefficient (gamma)for ClONO2 + HCl on H2SO4/H2O (binary aerosol) at 198K, Peter, (1997)
kb = 1.38e-23  # Boltzmann Constant
te1 = 198     # initial temperature (K)
mmcn = 97.46  # molar mass of chlorine natrate ClONO2 (g/mol)
mcn = (mmcn/an)/1000  # molecular mass of ClONO2 (kg/molec)
m = mcn
cgas = ((8*kb*te1)/(pi*m))**1/2  # mean molecular speed of gas phase ClONO2 (m/s)
cgsc = cgas*100  # mean molecular speed of gas phase ClONO2 in (cm/s)
print(cgsc)
sad = 8.6  # surface area density um^2 cm^-3 from Tilmes et al (2022) CESM2 data, multi-year average following initial 3 yr particle growth phase 
sadl = np.array([2.0,2.0,3.0,4.0,5.0,8.0,9.0,8.0,8.5,8.0,7.5,10.0,10.0,9.5,8.0,7.0,8.5,9.0,9.0,8.0,8.5,8.0,10.0]) # SAD data (SAI), Tilmes et al,(2022) for CESM2
sade = np.array([8.60,6.45,4.30,2.15,17.20,25.80,34.40,43.00,68.80,86.00])  # SAD error scenario values (CESM2 post growth avg x 1, .75, .5, .25, 2,3,4,5,8,10) um^2/cm^3

k  = 0.25*ut*cgas*sad  # Het reaction rate (rxns/cm^3 s) for (R1): ClONO2 + HCl --> Cl2 + HNO3; calculated using Wegner et al,(2012) Eq
kc = 0.25*ut*cgsc*sad  # Het reaction rate (rxns/cm^3 s) for (R1): ClONO2 + HCl --> Cl2 + HNO3; calculated using Wegner et al,(2012) Eq

def ksal(x,y = cgsc):
    return (0.25*ut*y*x)
a = sadl
for sad in (sadl):
    ksl = ksal(a) # gives reaction rate (rxns cm^-3 s^-1) for each year into SAI deployment, based on successive temporal SAD values following initial deployment
ksl = np.array([ksl])

print(k)    
cgas = ((8*kb*prm.temp)/(pi*m))**1/2  # mean molecular speed of gas phase ClONO2 (m/s)
print(cgas)
sad = 8.6  # surface area density um^2 cm^-3 from Tilmes et al (2022) CESM2 data, multi-year average following initial 5 yr particle growth phase 
#k = 0.25*ut*cgas*sad  # Het reaction rate (rxns/cm^3 s) for (R1): ClONO2 + HCl --> Cl2 + HNO3; calculated using Wegner et al,(2012) Eq

# sades = np.array([4.3,17.2,25.8,43.0,86.0])  # list of 'error margin scenario' SAD values (0.5,2,3,5,10 x SAD value from CESM2)  
sadh = 4.3  # half of CESM2 projected SAD
sad2 = 17.2  # 2 x CESM2 projected SAD
sad3 = 25.8
sad5 = 43.0
sad10 = 86.0
# Note: SAD values, from source research study, for a 2020 start, (and for each year into SAI deployment) are applied here to later start dates, 
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

tts = np.linspace(1,prm.year_e-prm.year_b+1,prm.year_e-prm.year_b+1)
sty = tts+prm.year_b-1
tse = np.array([0,2,7,12,17,22,27,32,37,42])   # yrs from 2023, five yr steps after 2025, extended to 2065

dO23 = -((c+(dc * t0)) * daec)*k*2  #a ozone depletion from (R1), 2023 SAI start (tts= 0) for 1K surface cooling, g cm^-3 s^-1
dO35 = -((c+(dc * t3)) * daec)*k*2  # ozone depletion from (R1), 2035 SAI start (tts=12) for 1K surface cooling, g cm^-3 s^-1
dO45 = -((c+(dc * t5)) * daec)*k*2  # ozone depletion from (R1), 2045 SAI start (tts=22) for 1K surface cooling, g cm^-3 s^-1

d23m = dO23/gmo  # ozone depletion from reaction R1 (2023 SAI start) in molec/cm^3 
d35m = dO35/gmo  # ozone depletion from reaction R1 (2035 SAI start) in molec/cm^3 
d45m = dO45/gmo  # ozone depletion from reaction R1 (2045 SAI start) in molec/cm^3 

d23s = d23m*dac  # ozone depletion (surface air equivalent) from R1 (2023 start) molec/cm^3, at STP
d23sr = np.around(d23s, 2)  # rounded to 2 decimal places
d35s = d35m*dac  # ozone depletion (surface air equivalent) from R1 (2035 start) molec/cm^3, at STP
d35sr = np.around(d35s, 2)
d45s = d45m*dac  # ozone depletion (surface air equivalent) from R1 (2045 start) molec/cm^3, at STP
d45sr = np.around(d45s, 2)

# d23l = -((c+(dc * t0)) * daec)*ksl*2*dac/gmo  # ozone depletion (surface air eq) from R1 (2023 start) molec/cm^3, at STP, variable SAD & k
# d25l = -((c+(dc * t1)) * daec)*ksl*2*dac/gmo  # ozone depletion (surface air eq) from R1 (2025 start) molec/cm^3, at STP, variable SAD & k
# d30l = -((c+(dc * t2)) * daec)*ksl*2*dac/gmo  # ozone depletion (surface air eq) from R1 (2030 start) molec/cm^3, at STP, variable SAD & k
# d35l = -((c+(dc * t3)) * daec)*ksl*2*dac/gmo  # ozone depletion (surface air eq) from R1 (2035 start) molec/cm^3, at STP, variable SAD & k
# d40l = -((c+(dc * t4)) * daec)*ksl*2*dac/gmo  # ozone depletion (surface air eq) from R1 (2040 start) molec/cm^3, at STP, variable SAD & k
# d45l = -((c+(dc * t5)) * daec)*ksl*2*dac/gmo  # ozone depletion (surface air eq) from R1 (2045 start) molec/cm^3, at STP, variable SAD & k
# --------------------------------------------------------------------------------------------------------

#  calculate stratospheric column ozone depletion in Dobson Units (DU)

du = 2.6867e20   # molec/m^2  in one Dobson Unit
d23du = (d23m*7000)/du  # O3 depletion from R1 (2023 start) in Dobson Units (using mid-strat column from 18-25 km)
d23rdu = np.around(d23du, 2)
d35du = (d35m*7000)/du  # O3 depletion from R1 (2035 start) in Dobson Units   "       "        "           "
d35rdu = np.around(d35du, 2)
d45du = (d45m*7000)/du  # O3 depletion from R1 (2045 start) in Dobson Units   "       "        "           "
d45rdu = np.around(d45du, 2)
# --------------------------------------------------------------------------------------------------------

# list of dO3 for 1 yr step successive future SAI start times (molec m^-3)

def odp(x, y = daec):
   return ((x+(dc*tts))*y)*k*2*dac/gmo
a = c
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

# plot ozone depletion vs Cl with SAD error scanarios:

def odps(x,y = daec2):
    def ksae(w,z = cgsc):
        return (0.25*ut*z*w)   
    a = sade
    for sad in sade:
        kse = ksae(a) # reaction rate (rxns cm^-3 s^-1) based on error scenario values for SAD
    kse = np.array([kse]) 
    
    return ((c+(dc*x))*y)*kse*2*dac/gmo   
a = tse   
for ts in tse:
    docs = odps(a)    
docs = np.array([docs])   

"""

