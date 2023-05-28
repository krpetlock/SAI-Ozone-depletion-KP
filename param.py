import numpy as np
import xlrd

clp = 3275     # total chlorine (Cl) in atmosphere in ppt as of 2023

#### reading NOAA EECS
eecs_noaa = xlrd.open_workbook("NOAA_EECS.xls") #reading NOAA xls for EECS
sh = eecs_noaa.sheet_by_index(0)
ln_rec = 178 #length of record from NOAA
eecs_y = np.zeros((ln_rec,1),float)
eecs = np.zeros((ln_rec,1),float)
for rx in range(15,ln_rec+15):
    eecs_y[rx-15] = float(sh.cell_value(rowx=rx, colx=0)) #reading year column
    eecs[rx-15] = float(sh.cell_value(rowx=rx, colx=21))   # reading EECS column

clp_noaa = np.mean(eecs[(ln_rec-4):ln_rec,:]) # average over last year

dcl = -8.23     # avg annual change in Cl (dissipation) in ppt/yr
clr = clp/1e12   # total chlorine as ratio of molecules of Cl/molec air
dclr = dcl/1e12   # avg annual change in Cl (dissipation) molec Cl/molec air

year_b = 2023 # first year in which depletion is simulated
year_e = 2045 # last year in which depletion is simulated

temp = np.linspace(192,202,21) # defining a temperature profile, here between 192 and 202 simply
dae_A = [1.76, 3.768, 6.52, 8.54, 10.71, 12.30, 13.79, 15.55, 16.28, 17.56, 18.04, 18.50, 19.70, 20.75, 21.93, 22.68, 23.48, 24.19, 
         24.94, 26.74, 26.44, 26.40, 27.68, 28.26, 28.84, 29.60, 30.65, 30.63, 32.37, 31.70, 32.52, 32.99, 33.90, 34.71, 34.78] #so4 mass values (ug/kg)

sadl = np.array([2.0,2.0,3.0,4.0,5.0,8.0,9.0,8.0,8.5,8.0,7.5,10.0,10.0,9.5,8.0,7.0,8.5,9.0,9.0,8.0,8.5,8.0,10.0]) # SAD data (SAI), Tilmes et al,(2022) for CESM2

gamma = [0.2, 0.22, 0.24] #making these coefficients up as an example