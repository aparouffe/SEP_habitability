#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 13:28:35 2021

@author: parouffe

computes the difference of oxygenation between 2006 and 2100
"""

from netCDF4 import Dataset
import glob
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colors import DivergingNorm

# depths = [0,100,200,300,500,800]
depths = [200]


path_fig = '/home/p/parouffe/Documents/figures/annexes/'
files = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/O2/*.nc')
f0 = Dataset(files[0])
lat = f0.variables['TLAT'][:]
lon = f0.variables['TLONG'][:]
zt = f0.variables['z_t'][:]
zt = zt/100
lon_ref = 360-110

lon2 = np.abs(lon[0,:]-lon_ref).argmin()

lat = lat[:,lon2:-3]
lon = lon[:,lon2:-3]

for depth in depths :
    conc_06 = []
    conc_21 = []
    depth_lvl = np.abs(zt-depth).argmin()
    
    for ff in files[:] :
        f = Dataset(ff)
        o2 = f.variables['O2'][0:12,depth_lvl,:,lon2:-3]
        conc_06.append(o2)
        o2 = f.variables['O2'][1128:,depth_lvl,:,lon2:-3]
        conc_21.append(o2)
    
    conc_06 = np.ma.array(conc_06)
    conc_2006 = np.ma.mean(conc_06,axis=0)  # yearly mean
    conc_2006 = np.ma.mean(conc_2006,axis=0) # ensemble mean
    
    conc_21 = np.ma.array(conc_21)
    conc_2100 = np.ma.mean(conc_21, axis=0)
    conc_2100 = np.ma. mean(conc_2100,axis=0)
    
    delta_O2 = conc_2100 - conc_2006
    print('delta max  ',np.max(delta_O2))





    vmin = -40; vmax=40
    levels = np.linspace(vmin,vmax,13)
    ticks =  np.arange(vmin,vmax+1,10)
    norm = DivergingNorm(vmin=vmin,vcenter=-0.001,vmax=vmax)
    
    print(depth_lvl)
    fig = plt.figure(figsize=(4.3,4))
    ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
    m = plt.contourf(lon,lat,delta_O2,\
                      vmin=vmin,vmax=vmax,levels=levels,
                      cmap='bwr',extend='both',norm=norm)
    cbar = plt.colorbar(m,ticks=ticks,shrink=0.8)
    cbar.ax.set_ylabel('mmol/m'+'\u00B3',fontsize=13)#,labelpad=5)
    # cbar.ax.tick_params(labelsize=12)
    ax.coastlines()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=1, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    plt.title('$\Delta$O2  (2100-2006) \n '+str(depth)+'m')
    plt.subplots_adjust(left=0.15,right=0.95)
    