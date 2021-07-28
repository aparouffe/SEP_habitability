#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 14:37:54 2021

@author: parouffe

computes phi along a North-south transect
"""

from netCDF4 import Dataset
import numpy as np             # version 1.19.2
import numpy.ma as ma
import xarray as xr

from scipy.stats import linregress
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from math import *
import glob
from matplotlib.colors import DivergingNorm

from metabolic_index_Z import *

ref = 82
year = 2100

# import nc files

path_fig = '/home/p/parouffe/Documents/graphs/MOYENNES_ENS/PHI/coupes_vert/'

fO2 = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/O2/*.nc')
ftemp = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/TEMP/*.nc')
fsalt = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/SALT/*.nc')

fO2.sort()
ftemp.sort()
fsalt.sort()

n = len(fO2)
"""   spatial variables + lon idx """
nc = Dataset(fO2[0])
zt = nc.variables["z_t"][:]
zt = zt/100

time = nc.variables['time'][:]  # extract days since...
lat = nc.variables['TLAT'][:]   # lat of T points
lon = nc.variables['TLONG'][:]  # lon of T points

lon_ref = 360-ref
idx_lon = np.abs(lon[0,:]-lon_ref).argmin()

a0 = 25
e0 = 0.4

zmax = 40 # 1000m

""" O2, T salt at lon ref """
phi_list = []
o2_list = []

if year == 2006:
    d1 = 0; d2 = 12
elif year == 2100:
    d1 = 1128; d2 = 1140

for fi in range(0,n):
    
        nc1 = Dataset(fO2[fi])
        nc2 = Dataset(ftemp[fi])
        nc3 = Dataset(fsalt[fi])
        
        c = nc1.variables['O2'][d1:d2,:zmax,:,idx_lon]
        t = nc2.variables['TEMP'][d1:d2,:zmax,:,idx_lon]
        s = nc3.variables['SALT'][d1:d2,:zmax,:,idx_lon]         
    
        c = np.ma.mean(c,axis=0)    
        t = np.ma.mean(t,axis=0)    
        s = np.ma.mean(s,axis=0)    
    
        phi = metabolic_index_Z (c,t,s,zt[:zmax],25,0.4)
        phi[phi<0] == 0
        # phi = np.ma.mean(phi,axis=0)
        phi_list.append(phi)

        o2_list.append(c)


""" ensemble mean """

o2 = np.ma.mean(o2_list,axis=0)

phi_final = np.ma.mean(np.ma.array(phi_list),axis=0)
phi_std = np.ma.std(np.ma.array(phi_list),axis=0)

phi_min = phi_final-phi_std
phi_max = phi_final+phi_std

phi_final[phi_final<0]=0

# ax1.clabel(m3,fmt='%2.1f $\Phi$ $_{crit}$',colors='r',)
# ax1.clabel(m2,fmt='%2.1f $\Phi$ $_{crit}$',colors='k',)

vmin = 0
vmax = 7
#customize colorbar
levels = np.linspace(vmin,vmax,15)
cbar_ticks = np.linspace(0,vmax,8)
y, z = np.meshgrid(lat[:,0],-zt[:zmax])
prof = zt[::10]



fig = plt.figure(figsize=(4,3.5))
# plt.suptitle('PHI - vertical section at '+str(ref)+'W \n'+str(year)+'m',fontsize = 15)
ax1 = fig.add_subplot(111)
ax1.set_title('2100 - 82°W')
# plot climate veolocities
m = plt.contourf(y,z,phi_final, 60, \
              cmap = 'turbo' ,vmin=vmin,vmax=vmax,levels=levels)
    #levels=levels,
m2 = ax1.contour(m,levels=[3.5],colors='k')
ax1.clabel(m2,fmt='%2.1f $\Phi$$_{crit}$',colors='k',fontsize=13)
m3 = ax1.contour(y,z,o2,levels=[45],colors='w')
ax1.clabel(m3,fmt='%2.1f mmol/m\u00B3',colors='w',fontsize=13)


crs = ccrs.PlateCarree()
# colorbar
cbar = plt.colorbar(m,ticks = cbar_ticks,pad=0.05)
cbar.ax.set_ylabel('$\Phi$',fontsize=15,labelpad = 1.5)
plt.ylim([-1000,0])
plt.xlim([-50,-5])
ax1.set_yticks(np.linspace(-1000,0,6))
ax1.set_yticklabels(['1000','800','600','400','200','0'])#[::-1])
plt.subplots_adjust(left = 0.19,right=0.99,bottom=0.13, top = 0.89,hspace=0.9 )  
plt.ylabel('depth (m)',fontsize = 14)
plt.xlabel('latitude', fontsize = 14)
plt.grid()
# plt.set_aspect('auto')
# plt.imshow(fig)
# plt.savefig(path_fig+str(year)+'_'+str(ref)+'W.png')

# ax.set_xticks(np.linspace(-45,-5,5))
# ax.set_xticklabels(['45S','35S','25S','15S','5S'][::-1])
# fig = plt.figure(figsize=(5,5))
# # plt.suptitle('PHI - vertical section at '+str(ref)+'W \n'+str(year)+'m',fontsize = 15)
# ax1 = fig.add_subplot(111)
# # plot climate veolocities
# m = plt.contourf(y,z,phi[:zmax,:], 60, \
#               cmap = 'YlOrRd_r' ,vmin=vmin,vmax=vmax,levels=levels)
#     #levels=levels,
# m2 = ax1.contour(m,levels=[3.5],colors='k')
# ax1.clabel(m2,fmt='%2.1f $\Phi$crit',colors='k')
# crs = ccrs.PlateCarree()
# # colorbar
# cbar = plt.colorbar(m,ticks = cbar_ticks,pad=0.05)
# cbar.ax.set_ylabel('$\Phi$',fontsize=15,labelpad = 1.5)
# plt.ylim([-1000,0])
# plt.xlim([-50,-5])
# plt.subplots_adjust(left = 0.19,right=1,bottom=0.2, top = 0.8,hspace=0.9 )  
# plt.ylabel('profondeur (m)',fontsize = 14)
# plt.xlabel('latitude', fontsize = 14)
# plt.grid()