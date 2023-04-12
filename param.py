import numpy as np

clp = 3275     # total chlorine (Cl) in atmosphere in ppt as of 2023
dcl = -8.23     # avg annual change in Cl (dissipation) in ppt/yr
clr = clp/1e12   # total chlorine as ratio of molecules of Cl/molec air
dclr = dcl/1e12   # avg annual change in Cl (dissipation) molec Cl/molec air

year_b = 2023 # first year in which depletion is simulated
year_e = 2045 # last year in which depletion is simulated

temp = np.linspace(192,202,21) # defining a temperature profile