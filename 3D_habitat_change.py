#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 16:30:25 2021

@author: parouffe

computes  habitat change on a 2D surface at every depth layer
+
perform a wilcoxon test
"""

from netCDF4 import Dataset
import numpy as np             # version 1.19.2

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
# from math import *
import glob
from matplotlib.colors import DivergingNorm
import matplotlib as mpl
from fonctions import metabolic_index_Z_3D
import matplotlib as mpl
from scipy.stats import wilcoxon


lon_ref = 360-110
# lon_ref2 = 360-80
idx2_lon = 360-110

zone = 'coast'  # 'NE', 'E', 'full','large'
save = 'no' # 'yes' if yes

years = ['2006','2100']

path_fig = '/home/p/parouffe/Documents/graphs/habitat_loss/'
# path_fig = '/home/alexandra/Documents/graphs//MOYENNES_ENS_PHI/PHI_CONFIG/'

fO2 = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/O2/*.nc')
ftemp = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/TEMP/*.nc')
fsalt = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/SALT/*.nc')

fO2.sort()
ftemp.sort()
fsalt.sort()

n = len(fO2)


"""   spatial variables """
nc = Dataset(fO2[0])
zt = nc.variables["z_t"][:]
zt = zt/100
dxt = nc.variables["DXT"][:]
dyt = nc.variables["DYT"][:]    
dxt = dxt/10**5    # to km
dyt = dyt/10**5    # to km
time = nc.variables['time'][:]  # extract days since...
lat = nc.variables['TLAT'][:]   # lat of T points
lon = nc.variables['TLONG'][:]  # lon of T points

idx_lon = np.abs(lon[0,:]-lon_ref).argmin()
# idx2_lon = np.abs(lon[0,:]-lon_ref2).argmin()
zmax = np.abs(zt-1200).argmin()
lat = lat[:,idx_lon:]
lon = lon[:,idx_lon:]
zt = zt[:zmax]

""" physiological parameters """
a0 = 25
e0 = 0.4
phi_crit = 3.5

""" extract variables for years 2006 and 2100 """

# o2_06 = []; o2_21 = []; o2_36 = []
# temp_06 = []; temp_21 = []; temp_36 = []
# salt_06 = []; salt_21 = []; salt_36 = []

phi_2006 = []
phi_2100 = []
v_change = []
phi_2036 = []
phi_2070 = []
v3D_06 = 0
v3D_21 = 0


n = len(fO2)
""" load and store all simulations """
for fi in range(0,n) :
    print(fi)
    i0 = 0
    nc1 = Dataset(fO2[i0])
    nc2 = Dataset(ftemp[i0])
    nc3 = Dataset(fsalt[i0])

    o2_06 = nc1.variables['O2'][0:360,:zmax,:,idx_lon:]
    temp_06 = nc2.variables['TEMP'][0:360,:zmax,:,idx_lon:]
    salt_06 = nc3.variables['SALT'][0:360,:zmax,:,idx_lon:]                                               

    o2_21 = nc1.variables['O2'][780:,:zmax,:,idx_lon:]
    temp_21 = nc2.variables['TEMP'][780:,:zmax,:,idx_lon:]
    salt_21 = nc3.variables['SALT'][780:,:zmax,:,idx_lon:]

    nc1.close()
    nc2.close()
    nc3.close()        
    """ mean """
    o2_06 = np.mean(o2_06, axis=0)  # 2nd, yearly mean      
    temp_06 = np.mean(temp_06,axis=0)
    salt_06 = np.mean(salt_06,axis=0)
    
    o2_21 = np.mean(o2_21, axis=0)  # 1st, mean over simulation      
    temp_21 = np.mean(temp_21,axis=0)
    salt_21 = np.mean(salt_21,axis=0)


    """ meatbolic index """
    
    phi_06 = metabolic_index_Z_3D(o2_06,temp_06,salt_06,zt[:zmax],a0,e0)
    phi_21 = metabolic_index_Z_3D(o2_21,temp_21,salt_21,zt[:zmax],a0,e0)

    
    phi_2006.append(phi_06)
    phi_2100.append(phi_21)


    
""" mean + std """

# phi_2006 = np.nanmean(phi_2006,axis=0)
# phi_2100 = np.nanmean(phi_2100,axis=0)



# nz,nlat,nlon = np.shape(phi_2006)
# std_2006 = np.zeros((zmax))
# std_2100 = np.zeros((zmax))


# phi_06_max = np.zeros((nz,nlat,nlon))
# phi_21_max = np.zeros((nz,nlat,nlon))


# phi_06_min = np.zeros((nz,nlat,nlon))
# phi_21_min = np.zeros((nz,nlat,nlon))

# for zz in range(zmax):
#     std_2006[zz] = np.nanstd(phi_2006[zz,:,:])
#     std_2100[zz] = np.nanstd(phi_2100[zz,:,:])


#     phi_06_max[zz] = phi_2006[zz,:,:] + std_2006[zz]
#     phi_06_min[zz] = phi_2006[zz,:,:] - std_2006[zz]

#     phi_21_max[zz] = phi_2100[zz,:,:] + std_2100[zz]
#     phi_21_min[zz] = phi_2100[zz,:,:] - std_2100[zz]   

   


nlat,nlon = np.shape(lat)
nsup_06 = np.zeros((n,zmax,nlat))
nsup_21 = np.zeros((n,zmax,nlat))

for sim in range(n):
    # print(sim)
    for zz in range(zmax):
        for ii in range(nlat):
            nsup_06[sim,zz,ii] = (phi_2006[sim][zz,ii,:] > phi_crit).sum()
            nsup_21[sim,zz,ii] = (phi_2100[sim][zz,ii,:] > phi_crit).sum()
            # nsup_36[zz] = (phi_2036[zz,:,:] > phi_crit).sum()
            # nsup_70[zz] = (phi_2070[zz,:,:] > phi_crit).sum()
        
""" wilcoxon test """
        
w_test = np.zeros((zmax,nlat))   
p_test = np.zeros((zmax,nlat))        

p_test[p_test==0]=np.nan # if the difference between the 2 series is equal to zero
                         # the test can't be applied
                         # the value in z,y is set to nan
w_test[w_test==0]=np.nan

d = nsup_06-nsup_21
dbis = np.zeros((n)) #to check if d different de zero
for zz in range(zmax):
    for ii in range(nlat):
        # d = nsup_06[:,zz,ii]-nsup_21[:,zz,ii]
               
        if (d[:,zz,ii] == dbis).all() == False :
            # print('yes',zz,ii)
            w,p = wilcoxon(nsup_06[:,zz,ii],nsup_21[:,zz,ii])#wilcoxon(d[:,zz,ii])
            w_test[zz,ii] = w
            p_test[zz,ii] = p   

# ntot = np.ma.count(phi_2006)
# print('n3D_2006 : ',n3D_06/ntot)    
# ntot = np.ma.count(phi_2100)
# print('n3D_2100 : ',n3D_21/ntot)


""" volume change """
nsup_06_mean = np.ma.mean(nsup_06,axis=0)
# nsup_06_mean = np.ma.sum(nsup_06_mean,axis=2)
nsup_21_mean = np.ma.mean(nsup_21,axis=0)
# nsup_21_mean = np.ma.sum(nsup_21_mean,axis=2)

volume_change = (nsup_21_mean-nsup_06_mean)/nsup_06_mean *100

# print(volume_change)

""" plot """

prof = np.arange(0,1001,200)
X,Y = np.meshgrid(lat[:,0],-zt)

fig,ax = plt.subplots(1,1,figsize=(4,3.5))
# m = plt.contourf(X[:,::-1],Y[:,::-1],volume_change[:,::-1],cmap='YlOrRd_r')
m = plt.pcolormesh(X,Y,volume_change,cmap='YlOrRd_r')
ax.set_ylim([-1000,-5])
ax.set_yticks(np.arange(-1000,1,200))
ax.set_yticklabels(prof[::-1])
cbar = plt.colorbar(m)
ax.set_ylabel('depth (m)', fontsize=13)
ax.set_xlabel('latitude', fontsize=13)
cbar.ax.set_ylabel('% habitat change',fontsize=13)
ax.set_facecolor('grey')
plt.subplots_adjust(top=0.98,left=0.18,right=0.89,bottom=0.14)

plt.show()
# ax.fill_betweenx(zt[:zmax],v_change_max,v_change_min,color='lightgrey')
# ax.plot(volume_change,zt[:zmax],color='k')
# ax = plt.gca()
# ax.invert_yaxis()
# ax.xaxis.tick_top()
# ax.xaxis.set_label_position('top')
# ax.grid(True)
# ax.tick_params(axis='y',labelsize=9)
# ax.tick_params(axis='x',labelsize=9)
# ax.set_ylabel('profondeur (m)',fontsize=9)
# ax.set_title('2100 relatif a 2006 \n % variation habitat',fontsize=9)
# plt.subplots_adjust(bottom=0.02,top=0.75,left=0.25,right=0.99)

# fig,ax = plt.subplots(1,1,figsize=(2.5,3))
# volume_change = (nsup_36-nsup_06)/nsup_06 *100
# v_change_max = (nsup_36_max-nsup_06)/nsup_06 *100
# v_change_min = (nsup_36_min-nsup_06)/nsup_06 *100
# ax.fill_betweenx(zt[:zmax],v_change_max,v_change_min,color='lightgrey')
# ax.plot(volume_change,zt[:zmax],color='k')
# ax = plt.gca()
# ax.invert_yaxis()
# ax.xaxis.tick_top()
# ax.tick_params(axis='y',labelsize=9)
# ax.tick_params(axis='x',labelsize=9)
# ax.xaxis.set_label_position('top')
# ax.grid(True)
# ax.set_title('2036 relatif a 2006 \n % variation habitat',fontsize=9)
# plt.subplots_adjust(bottom=0.02,top=0.75,left=0.25,right=0.99)

# fig,ax = plt.subplots(1,1,figsize=(2.5,3))
# volume_change = (nsup_21-nsup_70)/nsup_70 *100
# v_change_max = (nsup_21_max-nsup_70)/nsup_70 *100
# v_change_min = (nsup_21_min-nsup_70)/nsup_70 *100
# ax.fill_betweenx(zt[:zmax],v_change_max,v_change_min,color='lightgrey')
# ax.plot(volume_change,zt[:zmax],color='k')
# ax = plt.gca()
# ax.invert_yaxis()
# ax.xaxis.tick_top()
# ax.xaxis.set_label_position('top')
# ax.tick_params(axis='y',labelsize=9)
# ax.tick_params(axis='x',labelsize=9)
# ax.grid(True)
# # ax[2].set_ylabel('depth (m)')
# ax.set_title('2100 relatif a 2070 \n % variation habitat',fontsize=9)
# plt.subplots_adjust(bottom=0.02,top=0.75,left=0.22,right=0.99)

#     # # fig1,ax1 = plt.subplots()
#     # ax1.plot(nsup_06/ntot*100,zt[:zmax],label='2006')
#     # ax1.plot(nsup_21/ntot*100,zt[:zmax],label='2100')
#     # ax1 = plt.gca()
#     # ax1.invert_yaxis()
#     # # ax1.xaxis.tick_top()
#     # ax1.set_xlabel('volume of habitat %')
#     # ax1.xaxis.set_label_position('top')
#     # # ax1.grid(True)
#     # # ax1.set_ylabel('depth (m)')
#     # plt.legend()
    
    