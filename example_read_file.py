#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 17:03:44 2023

@author: danielevisioni
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

#the three lines above import libraries in Python that are necessary for reading NetCDF files and do math.

fold = '/Users/danielevisioni/Documents/GitHub/' #this is my folder
fname = 'b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.SSP245-MA-GAUSS-DEFAULT.002.cam.h0zm.' #this is the case name 
fname2 = '.203501-206912.nc'
var_clox = 'OddOx_CLOxBROx_Loss' #various variables you will use: this is chlorine-driven loss
var_o3 = 'O3' #ozone concentrations (mol/mol)
var_o3l = 'O3_Loss' #overall ozone loss

lat = xr.open_dataset(f"{fold}/{fname}{var_o3}{fname2}").lat #this is loading the variable "latitude" from the file
lev = xr.open_dataset(f"{fold}/{fname}{var_o3}{fname2}").lev #this is loading the variable "pressure levels" from the file
o3 = xr.open_dataset(f"{fold}/{fname}{var_o3}{fname2}").O3 #this is loading the variable "Ozone concentrations" from the file