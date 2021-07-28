#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 15:18:04 2021

@author: parouffe

This script computes the sensitivity of phi to O2 and Temperature along a vertical profile
sensitivity : coefficient of the derivative of phi
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

from metabolic_index_Z import *
from convers_O2_to_po2_Z import *



lon_ref= 75
lat_ref = 25
path_fig = '/home/p/parouffe/Documents/figures/3.2./cv_phi/'

fO2 = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/O2/*.nc')
ftemp = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/TEMP/*.nc')
fsalt = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/SALT/*.nc')

fO2.sort()
ftemp.sort()
fsalt.sort()


""" meatbolic index paramters """
a0 = 25
A0 = a0*1.013e-3     # mbar
e0 = 0.4         # eV
kB = 8.6e-5      # eV/K
Tref = 272.15+15 # K

nc = Dataset(fO2[0])
zt = nc.variables["z_t"][:]
zt = zt/100
lon = nc.variables['TLONG'][:]  # lon of T points        
lat = nc.variables['TLAT'][:]   # lat of T points

id_lon = np.abs(lon[0,:]-(360-lon_ref)).argmin()
id_lat = np.abs(lat[:,0]+lat_ref).argmin()
id_zmax = np.abs(zt-1200).argmin()
zt = zt[:id_zmax+1]


""" initialisation des variables """
list_phi = []
list_T =  []
list_po2 = []
list_a = []
list_b = []
list_dphi_approx = []
list_dphi_mdl = []
list_rms_dT = []


n = len(fO2)

for fi in range(0,n) :
    print(fi)
    """   spatial variables + depth level"""

    nc1 = Dataset(fO2[fi])
    nc2 = Dataset(ftemp[fi])
    nc3 = Dataset(fsalt[fi])

    o2 = nc1.variables['O2'][:,:id_zmax+1,id_lat,id_lon]
    Temp_temp = nc2.variables['TEMP'][:,:id_zmax+1,id_lat,id_lon]
    salt = nc3.variables['SALT'][:,:id_zmax+1,id_lat,id_lon]                                               # dim [time,z,lat=116,lon=97]

    # o2 = np.ma.mean(o2,axis=0)
    # Temp = np.ma.mean(Temp,axis=0)
    # salt = np.ma.mean(salt,axis=0)
    nc1.close()
    nc2.close()
    nc3.close()
    
        
    """ PHI  PO2 """
    
    po2_temp = convers_O2_to_po2_Z(o2,Temp_temp,salt,zt)   
    
    phi_ref = A0*po2_temp*np.exp(e0/kB*(1/(Temp_temp+272.15)-1/Tref))

    list_phi.append(phi_ref)
    list_po2.append(po2_temp)
    list_T.append(Temp_temp)

""" parametres initialisation """
phi = np.ma.array(list_phi)
Temp = np.ma.array(list_T)
po2 = np.ma.array(list_po2)

A0 = a0*1.01325e-3
E0 = e0  *1.60218e-19
kB = 1.38064852*1e-23
temp_K = Temp+273.15
Tref = 15+273.15

""" phi mean + std """
phi_mean = np.ma.mean(phi,axis=1) #moyenne temporelle 2006-2100
std_phi = np.nanstd(phi_mean,axis=0)   # ecart type sur l'ensemble
phi_mean = np.ma.mean(phi_mean,axis=0) # moyenne des simus

""" calcul dT et dpO2 : mean and for each sim at each depth """
ntime,nz = np.shape(Temp[0])
def temporal_trend(var):
    years = np.arange(0,ntime,12)
    nyears = len(years)
    var_years = np.zeros((nyears))
    i = 0
    for yy in years :
        var_years[i] = np.ma.mean(var[yy:yy+12],axis=0)
        i=i+1

    s,_,_,_,_ = linregress(np.arange(0,nyears),var_years)

    return s

trends_temp = np.zeros((n,nz)) # dim z
trends_po2 = np.zeros((n,nz))
trends_phi = np.zeros((n,nz))

for sim in range(n):
    # print(n)
    for zz in range(nz):
        # print(zz)
        trends_temp[sim,zz] = temporal_trend(Temp[sim,:,zz])
        trends_po2[sim,zz] = temporal_trend(po2[sim,:,zz])
        trends_phi[sim,zz] = temporal_trend(phi[sim,:,zz])

# tendances moyennes dim=z     
dT = np.ma.mean(trends_temp,axis=0)
dpo2 = np.ma.mean(trends_po2,axis=0)
dphi_mdl = np.ma.mean(trends_phi,axis=0)
std_dphi = np.nanstd(trends_phi,axis=0)
std_dT = np.nanstd(trends_temp,axis=0)
std_dpo2 = np.nanstd(trends_po2,axis=0)

# rmse à faire

""" calcul des coefficients a et b """
mean_Temp = np.ma.mean(Temp,axis=1) # axis = 1 : moyenne temporelle
mean_po2 = np.ma.mean(po2,axis=1)   # axis = 1 : moyenne temporelle
f = np.exp(E0 / kB *(1/(mean_Temp+273.15) - 1/Tref))  # ~1 ok
kbt2 = kB*((mean_Temp+273.15)**2)

# dphi // T = A0 * po2 * (_e0/kbT2) * f
# dphi // po2 = A0 * f

a = A0 * mean_po2 * (-E0/kbt2) * f  # dérivée partielle par rapport à T
b = A0 * f                     # dérivée partielle par rapport à pO2
    
std_a = np.nanstd(a,axis=0)    
std_b = np.nanstd(b,axis=0)

""" calcul de dPHI """
a_mean = np.ma.mean(a,axis=0)
b_mean = np.ma.mean(b,axis=0)
std_a = np.nanstd(a,axis=0)
std_b = np.nanstd(b,axis=0)

# dphi_approx = a_mean*dT + b_mean*dpo2
print(np.shape(a),np.shape(b))
print(np.shape(trends_temp),np.shape(trends_po2))

dphi_approx = np.ma.mean((a*trends_temp+b*trends_po2),axis=0)
std_dphi_approx = np.nanstd((a*trends_temp+b*trends_po2),axis=0)
# dphi_model = std_phi
print(np.shape(dphi_approx))
dpo2_max = dpo2.max()
dT_max = dT.max()

""" plot """
from matplotlib import gridspec

fig = plt.figure(figsize=(6.5,5.5))
spec = gridspec.GridSpec(ncols=2,nrows=1,
                         width_ratios = [2,1],wspace=0.1)
ax1 = fig.add_subplot(spec[0])
ax3 = fig.add_subplot(spec[1])#,sharey=ax1)
ax1.plot(a_mean,-zt,label=' a',color='green')
ax1.plot(b_mean,-zt,label=' b',color='blue')
ax2 = ax1.twiny()
ax2.plot(a_mean*dT,-zt,label=' a x dT',color='green',linestyle='dashed')
ax2.plot(b_mean*dpo2,-zt,label=' b x dpO2',color='blue',linestyle='dashed')
ax2.plot(dphi_approx,-zt,label='d$\Phi$ approx',color='k',linestyle='dashed')
ax2.plot(dphi_mdl,-zt,label='d$\Phi$ model',color='r',linestyle='dashed')
ax2.fill_betweenx(-zt,dphi_mdl+std_dphi,dphi_mdl-std_dphi,color='lightgrey')
ax1.legend(bbox_to_anchor=(0.7, 0))
ax1.set_yticks(np.linspace(-1000,0,6))
ax1.set_yticklabels(['1000','800','600','400','200','0'])
ax1.legend(ncol=2,bbox_to_anchor=(0.7, -0.06))
# ax2.legend()
ax2.legend(ncol=2,bbox_to_anchor=(0.88, 1.25))
# ax3 = plt.gca()
# ax3.invert_yaxis()
# ax3.xaxis.tick_top()
# ax3.xaxis.set_label_position('top')

# ax.fill_betweenx(-zt,dphi_mdl,v_change_min,color='lightgrey')
ax4 = ax3.twiny()
ax4.plot(dpo2,-zt,label='dpO2',color='blue')
ax4.fill_betweenx(-zt,dpo2+std_dpo2,dpo2-std_dpo2,color='lightblue')
# ax3.legend()
ax3.plot(dT,-zt,label='dT',color='green')
ax3.fill_betweenx(-zt,dT+std_dT,dT-std_dT,color='lightgreen')
ax3.legend(ncol=2,bbox_to_anchor=(0.8, -0.06))
ax4.legend(ncol=2,bbox_to_anchor=(0.8, 1.17))
# ax3.tick_params(levelleft=False)
ax1.set_ylim([-1000,0])
ax3.set_ylim([-1000,0])
ax1.grid(True)
ax1.set_ylabel('depth (m)',fontsize=13)
ax2.set_xticks( np.arange(-0.008,0.013,0.004))
ax2.set_xticklabels(np.round(np.arange(-0.008,0.013,0.004),3).astype(str).tolist())
ax3.grid(True)
plt.subplots_adjust(bottom=0.13,right=0.97,left=0.11,top=0.8)
ax3.set_yticklabels([])
plt.show()

fig,ax2 = plt.subplots()
ax2.plot(dphi_approx,-zt,label='d$\Phi$ approx',color='k',linestyle='dashed')
ax2.plot(dphi_mdl,-zt,label='d$\Phi$ model',color='r',linestyle='dashed')
ax2.fill_betweenx(-zt,dphi_mdl+std_dphi,dphi_mdl-std_dphi,color='lightgrey')
# ax2.fill_betweenx(-zt,dphi_approx+std_dphi_approx,dphi_approx-std_dphi_approx,color='lightcoral')

ax2.legend()
ax2.set_ylabel('depth (m)')
ax2.grid()