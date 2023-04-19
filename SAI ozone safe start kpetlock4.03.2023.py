import matplotlib.pyplot as plt
import numpy as np
import param as prm

# define initial chlorine (Cl) parameters

t0 = 0.0           # years in future from 2023 (initially set to 0 for 2023)
t1 = 2.0
t2 = 7.0            
t3 = 12.0
t4 = 17.0
t5 = 22.0

# remaining Cl (as ratio) in atm at time (t) yrs in future from 2023
clrt = (prm.clp + (prm.dcl * t0))/1e12
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

#  clkm = (mc20 * gmcl)  # Cl in g/m^3
#  clgc = clkm * 1e-6  # Cl in g/cm^3
# dclgc = (dmc20c * gmcl)  # annual change in Cl in g/cm^3 
# ------------------------------------------------------------------------

# calculate total volume of mid stratosphere, from 18-25 km altitude

re = (6371)  # radius of earth (km)
r25 = re+(25)  # radius of earth + 25 km altitude
r18 = re+(18)  # radius of earth + 18 km altitude
pi = 3.14
vol25 = (4.0/3.0)*pi*r25**3.0  # volume of sphere of radius, earth + 25 km altitude
vol18 = (4.0/3.0)*pi*r18**3.0  # volume of sphere of radius, earth + 18 km altitude
vmsk = vol25 - vol18  # volume of mid stratosphere (km^3) from 18 to 25 km altitude
vmsm = (vmsk)*10**9  # volume of mid stratosphere in m^3

sin60 = 0.866
h25 = r25-(r25*sin60)
h18 = r18-(r18*sin60)
vlc = (2/3)*(pi*((r25)**2)*h25)  # volume of larger cone (60-90s), Re + 25 km
vsc = (2/3)*(pi*((r18)**2)*h18)  # volume of smaller cone (60-90s), Re + 18 km 
vams = vlc - vsc # volume of Antarctic mid strat 18-25 km altitude (60-90s) (km^3)
vamsm = (vams)*10**9  # volume of Antarctic mid strat 18-25 km altitude (60-90s) (m^3)
# -------------------------------------------------------------------------

# calculate mass density of sulfate aerosol to add to mid stratosphere for 1 K surface cooling, in g/cm^3

dela = 2.0e9 # total mass (kg) of additional aerosol needed per yr for 1K cooling, 2 Tg (until 2045), from CESM2-WACCM
dae = dela/vamsm  # mass density of added aerosol in Antarctic mid stratosphere, kg/m^3
daec = (dae *1000) * 1e-6  # mass density of added aerosol in Antarctic mid stratosphere, g/cm^3

# calculate mass density from CESM values (ARISE), take year 2 as test

dae_A2 = prm.dae_A[2]*da20 
dae_A2 = dae_A2 * 1e-6 * 1e-6

# -------------------------------------------------------------------------

# calculate heterogeneous reaction rate for R1

# k1 = 1000  # reaction rate for (R1); (rxns cm^-3 s^-1) at 200K in mid stratosphere, Borrmann et al, (1997)
ut = 0.2  # uptake coefficient (gamma)for ClONO2 + HCl on H2SO4/H2O (binary aerosol) at 198K, Peter, (1997)
kb = 1.38e-23  # Boltzmann Constant
te1 = 198     # initial temperature (K)
mcn = (36.46/an)/1000  # molecular mass of HCl (kg/molec)
m = mcn
cgas = ((8*kb*te1)/(pi*m))**1/2  # mean molecular speed of gas (m/s?)
print(cgas)
sad = 8.6  # surface area density um^2 cm^-3 from Tilmes et al (2022) CESM2 data, multi-year average following initial 5 yr particle growth phase 
k = 0.25*ut*cgas*sad  # Het reaction rate (rxns/cm^3 s) for (R1): ClONO2 + HCl --> Cl2 + HNO3; calculated using Wegner et al,(2012) Eq
kh = 0.5*k
k2 = 2*k

# sadl = np.array([2.0,2.0,3.0,4.0,5.0,8.0,9.0,8.0,8.5,8.0,7.5,10.0,10.0,9.5,8.0,7.0,8.5,9.0,9.0,8.0,8.5,8.0,10.0]) # SAD data (SAI), Tilmes et al,(2022) for CESM2
# sades = np.array([4.3,17.2,25.8,43.0,86.0])
# sadh = 4.3
# sad2 = 17.2 
# sad3 = 25.8
# sad5 = 43.0
# sad10 = 86.0
# Note: SAD values, from source research study, for a 2020 start, (and for each year into SAI deployment) are applied here to later start dates, 
# assuming conditions (such as temperature) have not changed significantly since 2020 to affect SAD values.

# def ksal(x,y = cgas):
#    return (0.25*ut*y*x)
# a = sadl
# for sad in (sadl):
#    ksl = ksal(a) # gives reaction rate (rxns m^-3 s^-1) for each year into SAI deployment, based on successive temporal SAD values following initial deployment
# ksl = np.array([ksl])

# class A(param.Parameterized)
#     title = param.String(default = "rxnrate", doc = "Reaction rate (k) rxns/cm^3 sec")
# class B(A)
#     a = param.Temperature(192, bounds=(192,205), doc = "First temp")
#     b = parem.UptakeR1(0.001, bounds=(0.001, 10), doc = "Fist utc")
#     c = param.SAD(2.0, bounds=(2,10), doc = "First SAD")
#     d = param.tts(0, bounds=(0,22), doc = "Time from 2023 until start")
#     e = param.Cl(3.5e-13, bounds(0.00,0.001), doc = "Cl at 20km  g/cm^3")
    
print(k)    
print('')     
# ---------------------------------------------------------------------------

# # calculate ozone depletion resulting from addition of new aerosol in SAI 1K cooling scenario, variable start times

tts = np.linspace(1,prm.year_e-prm.year_b+1,prm.year_e-prm.year_b+1)
sty = tts+prm.year_b-1

dO23 = -((c+(dc * t0)) * daec)*k*2  #a ozone depletion from (R1), 2023 SAI start (tts= 0) for 1K surface cooling, g cm^-3 s^-1
dO35 = -((c+(dc * t3)) * daec)*k*2  # ozone depletion from (R1), 2035 SAI start (tts=12) for 1K surface cooling, g cm^-3 s^-1
dO45 = -((c+(dc * t5)) * daec)*k*2  # ozone depletion from (R1), 2045 SAI start (tts=22) for 1K surface cooling, g cm^-3 s^-1

d23m = dO23/gmo  # ozone depletion from reaction R1 (2023 SAI start) in molec/cm^3 
d35m = dO35/gmo  # ozone depletion from reaction R1 (2035 SAI start) in molec/cm^3 
d45m = dO45/gmo  # ozone depletion from reaction R1 (2045 SAI start) in molec/cm^3 

d23s = d23m*dac  # ozone depletion (surface air equivalent) from R1 (2023 start) molec/cm^3, at STP
d23sr = round(d23s, 2)  # rounded to 2 decimal places
d35s = d35m*dac  # ozone depletion (surface air equivalent) from R1 (2035 start) molec/cm^3, at STP
d35sr = round(d35s, 2)
d45s = d45m*dac  # ozone depletion (surface air equivalent) from R1 (2045 start) molec/cm^3, at STP
d45sr = round(d45s, 2)

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
d23rdu = round(d23du, 2)
d35du = (d35m*7000)/du  # O3 depletion from R1 (2035 start) in Dobson Units   "       "        "           "
d35rdu = round(d35du, 2)
d45du = (d45m*7000)/du  # O3 depletion from R1 (2045 start) in Dobson Units   "       "        "           "
d45rdu = round(d45du, 2)
# --------------------------------------------------------------------------------------------------------

 # [-(((c+(dc*wt))*daec)*k*2)*dac/gmo]  # list of dO3 for 1 yr step successive future SAI start times (molec m^-3)

def odp(x, y = daec):
   return ((x+(dc*tts))*y)*k*2*dac/gmo
a = c
for wt in (tts):
    dO3t = odp(a)
    
print(dO3t)
print('')

dO3y = (
    {2023:dO3t[0], 
    2024:dO3t[1], 
    2025:dO3t[2], 
    2026:dO3t[3], 
    2027:dO3t[4],
    2028:dO3t[5],
    2029:dO3t[6],
    2030:dO3t[7],
    2031:dO3t[8],
    2032:dO3t[9],
    2033:dO3t[10], 
    2034:dO3t[11],
    2035:dO3t[12],
    2036:dO3t[13],
    2037:dO3t[14], 
    2038:dO3t[15], 
    2039:dO3t[16],
    2040:dO3t[17],
    2041:dO3t[18],
    2042:dO3t[19],
    2043:dO3t[20], 
    2044:dO3t[21],
    2045:dO3t[22]}
    )
print(dO3y)


# print outputs with string ('words.......') labels

print('The value of k, (rxns/cm^3 s) calculated using Eq from Wegner et al, (2012) is:')
print(k)
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

plt.close()
fig=plt.figure(figsize=[12,6])
plt.plot(sty,dO3t,'o', markersize=1, label='1K cooling')
plt.xlabel('start date', fontsize=14)
plt.ylabel('dO3 (-molec cm^-3 s^-1)', fontsize=14)    
plt.title('SAI O3 depletion (Cl dependent) by start date', fontsize=18)
plt.show()

plt.plot()


# x = sty

# plt.plot(x, d25l)
# plt.plot(x, d30l)
# plt.plot(x, d35l)
# plt.plot(x, d40l)
# plt.plot(x, d45l)
# plt.show()
