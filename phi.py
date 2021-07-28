#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 11:20:24 2021

@author: parouffe

computes phi + std phi
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
from scipy.stats import wilcoxon

from metabolic_index_ST import *

# depths = [0,100,300,500,800]    # 1 to 60
depth = 200

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

o2_list_2100 = []

for fi in range(n) :

    """   spatial variables + depth level"""
    print(fi)
    nc1 = Dataset(fO2[fi])
    nc2 = Dataset(ftemp[fi])
    nc3 = Dataset(fsalt[fi])
    # print(fO2[fi])
    # print(ftemp[fi])
    # print(fsalt[fi])
    
    o2_temp_2006 = nc1.variables['O2'][:12,depth_lvl,:,lon2:-3]
    Temp_temp_2006 = nc2.variables['TEMP'][:12,depth_lvl,:,lon2:-3]
    salt_2006 = nc3.variables['SALT'][:12,depth_lvl,:,lon2:-3]                                               # dim [time,z,lat=116,lon=97]

    o2_temp_2100 = nc1.variables['O2'][1128:,depth_lvl,:,lon2:-3]
    Temp_temp_2100 = nc2.variables['TEMP'][1128:,depth_lvl,:,lon2:-3]
    salt_2100 = nc3.variables['SALT'][1128:,depth_lvl,:,lon2:-3]  
    
    nc1.close()
    nc2.close()
    nc3.close()
    
    o2_temp_2006 = np.ma.mean(o2_temp_2006,axis=0)
    Temp_temp_2006 = np.ma.mean(Temp_temp_2006,axis=0)
    salt_2006 = np.ma.mean(salt_2006,axis=0)
    
    o2_temp_2100 = np.ma.mean(o2_temp_2100,axis=0)
    Temp_temp_2100 = np.ma.mean(Temp_temp_2100,axis=0)
    salt_2100 = np.ma.mean(salt_2100,axis=0)
    
    
    """ PHI using metabolic_index function """
    # def meatabolic_index(O2,Temp,salt,z_t,depth_lvl)
    phi_temp_2100 = metabolic_index_ST(o2_temp_2100,Temp_temp_2100,salt_2100,zt,depth_lvl,a0,e0)
    phi_list_2100.append(phi_temp_2100)

    phi_temp_2006 = metabolic_index_ST(o2_temp_2006,Temp_temp_2006,salt_2006,zt,depth_lvl,a0,e0)
    phi_list_2006.append(phi_temp_2006)
    
    o2_list_2100.append(o2_temp_2100)

""" list  to array """
o2_2100 = np.ma.mean(o2_list_2100,axis=0)

phi_2006 = np.ma.array(phi_list_2006)
phi_2100 = np.ma.array(phi_list_2100)    

""" wilcoxon test """
ns,nlat,nlon = np.shape(phi_2100)        
w_test = np.zeros((nlat,nlon))   
p_test = np.zeros((nlat,nlon))        

p_test[p_test==0]=np.nan # if the difference between the 2 series is equal to zero
                         # the test can't be applied
                         # the value in z,y is set to nan
w_test[w_test==0]=np.nan

d = phi_2006-phi_2100
dbis = np.zeros((n)) #to check if d different de zero
for ii in range(nlat):
    for jj in range(nlon):
        # d = nsup_06[:,zz,ii]-nsup_21[:,zz,ii]
               
        if (d[:,ii,jj] == dbis).all() == False :
            # print('yes',zz,ii)
            w,p = wilcoxon(phi_2006[:,ii,jj],phi_2100[:,ii,jj])#wilcoxon(d[:,zz,ii])
            w_test[ii,jj] = w
            p_test[ii,jj] = p   


""" ensemble mean + std """
# delta = phi_list_2100[0]-phi_list_2006[1]

phi_std_2006 = np.nanstd(phi_2006,axis=0)
phi_2006 = np.ma.mean(phi_2006,axis=0)
phi_2006[phi_2006<0] = 0

phi_std_2100 = np.nanstd(phi_2100,axis=0)
phi_2100 = np.ma.mean(phi_2100,axis=0)
phi_2100[phi_2100<0] = 0

phi_min = phi_2100-phi_std_2100
phi_max = phi_2100+phi_std_2100

"""" wilcoxon test """

lon = lon[:,:-3]
lat=lat[:,:-3]
i_lat,i_lon = np.where(p_test>=0.05)

# stat= np.zeros((nlat,nlon))
# stat[stat==0] = np.nan
# stat[i_lat,i_lon]=100

# x,y = np.meshgrid(lon[0,(i_lon)],lat[(i_lat),0])

# fig1, ax1 = plt.subplots(1,1, subplot_kw={'projection':ccrs.PlateCarree()},figsize=(4,4))
# ax1.set_title(str(depth)+'m')
# m1 = ax1.contourf(lon,lat,phi_2100,cmap='turbo')
# ax1.scatter(lon[0,i_lon],lat[i_lat,0],marker='.',color='k',s=4)
# ax1.coastlines()
# gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                   linewidth=1, color='gray', alpha=1, linestyle='--')
# gl.xlabels_top = False
# gl.ylabels_right = False
# cbar = plt.colorbar(m1,shrink=0.8)
# plt.subplots_adjust(bottom=0.05,top=0.95,right=0.98)


""" plot """
fig1, ax1 = plt.subplots(1,1, subplot_kw={'projection':ccrs.PlateCarree()},figsize=(3.5,3.5))
ax1.set_title('2100 - 200m')
levels = np.linspace(0,6,13)# np.arange(vmin,vmax+1,2.5)
ticks = np.linspace(0,6,7) #np.arange(vmin,vmax+1,5)
m1 = ax1.contourf(lon,lat,phi_2100,vmin=0,vmax=7,levels=levels,
                cmap='turbo')
m2 = ax1.contour(m1,levels=[3.5],colors='k')
m3 = ax1.contour(lon,lat,phi_2006,levels=[3.5],colors='r')
m4 = ax1.contour(lon,lat,o2_2100,levels=[45],colors='w')
m4 = ax1.contour(lon,lat,phi_min,levels=[3.5],colors='k',linestyles='dashed')
m4 = ax1.contour(lon,lat,phi_max,levels=[3.5],colors='k',linestyles='dotted')
ax1.scatter(281,-34,marker='*',color='red')
ax1.scatter(281,-25,marker='*',color='red')
ax1.scatter(251,-27,marker='*',color='red')
ax1.clabel(m3,fmt='%2.1f $\Phi$ $_{crit}$',colors='r',fontsize=13)
ax1.clabel(m2,fmt='%2.1f $\Phi$ $_{crit}$',colors='k',fontsize=13)

plt.subplots_adjust(bottom=0.02,top=0.98,right=0.95)

cbar1 = plt.colorbar(m1,ax=ax1,ticks=ticks,shrink=0.7,pad=0.065)
cbar1.ax.set_ylabel('$\Phi$',labelpad= 0.25,fontsize=13)
ax1.coastlines()
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=1, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False

plt.show()

""""""

# fig1, ax1 = plt.subplots(1,1, subplot_kw={'projection':ccrs.PlateCarree()})#,figsize=(3.5,3.5))
# ax1.set_title('2100 - 200m')
# vmin=-1;vmax=1
# levels = np.linspace(-1,1,11)# np.arange(vmin,vmax+1,2.5)
# ticks = np.linspace(-1,1,6) #np.arange(vmin,vmax+1,5)
# m1 = ax1.contourf(lon,lat,delta,vmin=vmin,vmax=vmax,levels=levels,
#                cmap='coolwarm')
# m2 = ax1.contour(lon,lat,phi_list_2100[0],levels=[3.5],colors='k')
# m3 = ax1.contour(lon,lat,phi_list_2100[1],levels=[3.5],colors='r')
# plt.colorbar(m1,ax=ax1,ticks=ticks)
# ax1.clabel(m3,fmt='$\Phi$ $_{crit}$ simu 1',colors='r',fontsize=13)
# ax1.clabel(m2,fmt='$\Phi$ $_{crit}$ simu 2',colors='k',fontsize=13)


# ax1.coastlines()
# gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                   linewidth=1, color='gray', alpha=1, linestyle='--')
# gl.xlabels_top = False
# gl.ylabels_right = False

# plt.show()


