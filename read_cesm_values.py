import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

fold = 'CESM_data'
fname = 'b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.SSP245-MA-GAUSS-DEFAULT.002.cam.h0zm.'
fname2 = '.203501-206912.nc'
var_clox = 'OddOx_CLOxBROx_Loss'
var_o3 = 'O3'
var_o3l = 'O3_Loss'

lat = xr.open_dataset(f"{fold}/{fname}{var_clox}{fname2}").lat
lev = xr.open_dataset(f"{fold}/{fname}{var_clox}{fname2}").lev
clox = xr.open_dataset(f"{fold}/{fname}{var_clox}{fname2}").OddOx_CLOxBROx_Loss

weights = np.cos(np.deg2rad(lat))  
weights.name = "weights"
lb = 0 #latitude bound 1, -90
lB = 31 #latitude bound 2, -60
clox_gm = clox[:,:,lb:lB].weighted(weights[lb:lB]).mean("lat").values #take latitude weighted average

clox_gm_mn = np.reshape(clox_gm,(35,12,70))
clox_oct = clox_gm_mn[:,9,:] # select month of october
clox_oct_150 = clox_oct[:,50] #select 150 hPa, max of ClOx depletion

years_m = np.linspace(2035,2069,35*12)
years = np.linspace(2035,2069,35)

fig = plt.figure(99,figsize=[24,12])
plt.plot(years,clox_oct_150,'x', markersize=10)
plt.xlabel('start date', fontsize=14)
plt.ylabel('dO3 (-molec cm^-3 s^-1)', fontsize=14)    
plt.title('ClOx driven O3 Loss in CESM2 at 150 hPa, October [60S-90S]', fontsize=18)
plt.savefig('ClOx_CESM.png')
plt.show()
