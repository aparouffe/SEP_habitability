#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 11:20:24 2021

@author: parouffe

computes phi
"""
from netCDF4 import Dataset
import numpy as np             # version 1.19.2
import numpy.ma as ma

from scipy.stats import linregress
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from math import *
import glob
from matplotlib.colors import DivergingNorm

from metabolic_index_ST import *
from convers_O2_to_po2 import convers_O2_to_po2

# depths = [0,100,300,500,800]    # 1 to 60
depth = 500

lon_ref=110

path_fig = '/home/p/parouffe/Documents/figures/3.2./cv_phi/'

fO2 = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/O2/*.nc')
ftemp = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/TEMP/*.nc')
fsalt = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/SALT/*.nc')

fO2.sort()
ftemp.sort()
fsalt.sort()

n = len(fO2)

""" meatbolic index paramters """
a0 = 25
e0 = 0.4
A0 = a0*1.01325e-3

nc = Dataset(fO2[0])
zt = nc.variables["z_t"][:]
zt = zt/100
depth_lvl = np.abs(zt-depth).argmin()
lon = nc.variables['TLONG'][:,:]  # lon of T points        
lon2 = np.abs(lon[0,:]-(360-lon_ref)).argmin()
lon = lon[:,lon2:]
lat = nc.variables['TLAT'][:,lon2:]   # lat of T points
dxt = nc.variables["DXT"][:,lon2:]
dyt = nc.variables["DYT"][:,lon2:]    
dxt = dxt/10**5    # to km
dyt = dyt/10**5    # to km



phi_list_2006 = []
phi_list_2100 = []

list_o2 = []
list_temp = []
list_salt = []

for fi in range(0,n) :

    """   spatial variables + depth level"""
    print(fi)
    nc1 = Dataset(fO2[fi])
    nc2 = Dataset(ftemp[fi])
    nc3 = Dataset(fsalt[fi])
    # print(fO2[fi])
    # print(ftemp[fi])
    # print(fsalt[fi])
    
    o2_temp_2006 = nc1.variables['O2'][:,depth_lvl,:,lon2:]
    Temp_temp_2006 = nc2.variables['TEMP'][:,depth_lvl,:,lon2:]
    salt_2006 = nc3.variables['SALT'][:,depth_lvl,:,lon2:]                                               # dim [time,z,lat=116,lon=97]

    # o2_temp_2100 = nc1.variables['O2'][1128:,depth_lvl,:,lon2:]
    # Temp_temp_2100 = nc2.variables['TEMP'][1128:,depth_lvl,:,lon2:]
    # salt_2100 = nc3.variables['SALT'][1128:,depth_lvl,:,lon2:]  
    
    nc1.close()
    nc2.close()
    nc3.close()
    
    o2_temp_2006 = np.ma.mean(o2_temp_2006,axis=0)
    Temp_temp_2006 = np.ma.mean(Temp_temp_2006,axis=0)
    salt_2006 = np.ma.mean(salt_2006,axis=0)
    
    list_o2.append(o2_temp_2006)
    list_temp.append(Temp_temp_2006)
    list_salt.append(salt_2006)
    
    # o2_temp_2100 = np.ma.mean(o2_temp_2100,axis=0)
    # Temp_temp_2100 = np.ma.mean(Temp_temp_2100,axis=0)
    # salt_2100 = np.ma.mean(salt_2100,axis=0)
    
    
    """ PHI using metabolic_index function """
    # def meatabolic_index(O2,Temp,salt,z_t,depth_lvl)
    # phi_temp_2100 = metabolic_index_ST(o2_temp_2100,Temp_temp_2100,salt_2100,zt,depth_lvl,a0,e0)
    # phi_list_2100.append(phi_temp_2100)

    phi_temp_2006 = metabolic_index_ST(o2_temp_2006,Temp_temp_2006,salt_2006,zt,depth_lvl,a0,e0)
    phi_list_2006.append(phi_temp_2006)
    
""" ensemble mean + std """
phi_2006 = np.ma.array(phi_list_2006)
# phi_2006 = np.ma.mean(phi_2006,axis=1)
phi_std_2006 = np.nanstd(phi_2006,axis=0)
phi_2006 = np.ma.mean(phi_2006,axis=0)
phi_2006[phi_2006<0] = 0

# phi_2100 = np.ma.array(phi_list_2100)
# # phi_2100 = np.ma.mean(phi_2100,axis=1)
# phi_std_2100 = np.nanstd(phi_2100,axis=0)
# phi_2100 = np.ma.mean(phi_2100,axis=0)
# phi_2100[phi_2100<0] = 0

# std_o2 = np.nanstd(o2_list,axis=0)
# o2 = np.ma.mean(o2_list,axis=0) # temporal mean
""" conversion o2 to po2 """
list_po2 = []
for ii in range(len(list_o2)):
    po2_temp = convers_O2_to_po2(list_o2[ii],list_temp[ii],list_salt[ii],zt,depth_lvl)
    list_po2.append(po2_temp)

""" ensemble mean po2 et T """

po2 = np.ma.mean(list_po2,axis=0)
std_po2 = np.nanstd(list_po2,axis=0)

temp = np.ma.mean(list_temp,axis=0)
std_temp = np.nanstd(list_temp,axis=0)

salt = np.ma.mean(list_salt,axis=0)


""" computation of O2 and temperature contribution to phi """

#### dphi = dphidT + dphidpo2 ###
A0 = a0*1.01325e-3
E0 = e0  *1.60218e-19
kB = 1.38064852*1e-23
temp_K = temp+273.15
Tref = 15+273.15

f = np.exp(E0 / kB *(1/temp_K - 1/Tref))  # ~1 ok
kbt2 = kB*temp_K*temp_K

# dphi // T = A0 * po2 * (_e0/kbT2) * f
# dphi // po2 = A0 * f

dphi_dT = A0 * po2 * (-E0/kbt2) * f #* ((std_temp+273.15)/temp_K)
dphi_dpo2 = A0 * f #* (std_po2/po2)

""" plot """
fig, ax = plt.subplots(1,2, subplot_kw={'projection':ccrs.PlateCarree()},figsize =(7,3.5))
ax[0].set_title('Sensitivity to oxygen')
ax[1].set_title('Sensitivity to temperature')
cmap='bwr'
vmin=0.01; vmax=0.04
levels=np.linspace(vmin,vmax,9);print(levels)
ticks=np.linspace(vmin,vmax,5);print(ticks)
norm=DivergingNorm(vmin=vmin,vcenter=0.011,vmax=vmax)
m1 = ax[0].contourf(lon,lat,dphi_dpo2,cmap=cmap,extend='both',
                    vmin=vmin,vmax=vmax,levels=levels,norm=norm)
cbar1 = plt.colorbar(m1,ax=ax[0],shrink=0.7,ticks=ticks)
vmin=-0.3; vmax=0
levels=np.linspace(vmin,vmax,11);print(levels)
ticks=np.linspace(vmin,vmax,6)
norm=DivergingNorm(vmin=vmin,vcenter=-0.001,vmax=vmax)
m2 = ax[1].contourf(lon,lat,dphi_dT,cmap=cmap,extend='both',
                    vmin=vmin,vmax=vmax,levels=levels,norm=norm)
cbar2 = plt.colorbar(m2,ax=ax[1],shrink=0.7,ticks=ticks)

for i in range(2):
    gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              linewidth=1, color='gray', alpha=1, linestyle='--')
    ax[i].coastlines()
    gl.xlabels_top = False
    gl.ylabels_right = False
plt.subplots_adjust(right=0.95,left=0.1,bottom=0.2,wspace=0.4)
# cbar_ax = plt.gcf().add_axes([0.25, 0.15, 0.5, 0.04])
# cbar = plt.colorbar(m2,pad=0.12,ax=ax[1],ticks=ticks,orientation='horizontal',cax=cbar_ax)
# cbar.ax.set_xlabel('%',fontsize=13)#.set_rotation(0)
# cbar.ax.tick_params(labelsize=12)
# m1 = ax[0].contourf(lon, lat, phi_21, 60,
#               transform=ccrs.PlateCarree(), cmap='turbo',\
#                   vmin=vmin,vmax=vmax,levels=levels,extend=extend)
# m11 = ax[0].contour(m1,levels=[2.3],colors=cc,fontsize=12)
# m12 = ax[0].contour(lon,lat,o2_2100,levels=[90],colors='r',fontsize=12)
# ax[0].clabel(m12,fmt='%2.1f mmol/m \u00B3',colors='r')
# ax[0].clabel(m11,fmt='%2.1f $\Phi$ crit',colors='k')


