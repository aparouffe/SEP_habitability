#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 12:33:48 2021

@author: parouffe

computes std cv phi / cv phi
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
from matplotlib.colors import DivergingNorm
import matplotlib as mpl
from metabolic_index_ST import *

var_name = "PHI"
# depths = [0,100,200,300,500,800]    # 1 to 60
depth = 200


path_fig = '/home/p/parouffe/Documents/figures/3.2./'

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


nc = Dataset(fO2[0])
zt = nc.variables["z_t"][:]
zt = zt/100
depth_lvl = np.abs(zt-depth).argmin()        


time = nc.variables['time'][:]  # extract days since...
lat = nc.variables['TLAT'][:]   # lat of T points
lon = nc.variables['TLONG'][:]  # lon of T points

lon_ref = 360-110
idx_lon = np.abs(lon[0,:]-lon_ref).argmin()

lat = lat[:,idx_lon:]
lon = lon[:,idx_lon:]

dxt = nc.variables["DXT"][:,idx_lon:]
dyt = nc.variables["DYT"][:,idx_lon:]    
dxt = dxt/10**5    # to km
dyt = dyt/10**5    # to km


""" initialisation of variables """
    

cv = []
gradS = []; gradT = []
gradNS = []; gradWE = []
angles = []


for fi in range(0,n) :

    print(fi)
    """   spatial variables + depth level"""

    nc1 = Dataset(fO2[fi])
    nc2 = Dataset(ftemp[fi])
    nc3 = Dataset(fsalt[fi])
    # print(fO2[fi])
    # print(ftemp[fi])
    # print(fsalt[fi])
    
    o2 = nc1.variables['O2'][:360,depth_lvl,:,idx_lon:]
    Temp = nc2.variables['TEMP'][:360,depth_lvl,:,idx_lon:]
    salt = nc3.variables['SALT'][:360,depth_lvl,:,idx_lon:]                                               # dim [time,z,lat=116,lon=97]

    nc1.close()
    nc2.close()
    nc3.close()
    
    ntime, nlat, nlon = np.shape(o2)

    
    """ PHI using metabolic_index function """
    # def meatabolic_index(O2,Temp,salt,z_t,depth_lvl)
    var = metabolic_index_ST(o2,Temp,salt,zt,depth_lvl,a0,e0)

    """ Spatial Gradient """
    var_tmean = np.mean(var,axis=0)  # mean on time axis
    
    dx = dyt; dy = dxt
    gradS_df = np.zeros((nlat,nlon))
    gradS_WE = np.zeros((nlat,nlon))
    gradS_NS = np.zeros((nlat,nlon))
    ang = np.zeros((nlat,nlon))
    for ni in range(1,nlat-1):
        for nj in range(1,nlon-1):
            WE = (var_tmean[ni,nj+1]-var_tmean[ni,nj-1])/(2*dx[ni,nj])
            NS = (var_tmean[ni+1,nj]-var_tmean[ni-1,nj])/(2*dy[ni,nj])
            
            if np.array(WE) == 0:
                WE = (var_tmean[ni,nj]-var_tmean[ni,nj-1])/(dx[ni,nj-1])
            
            gradS_WE[ni,nj] = WE
            gradS_NS[ni,nj] = NS       
            gradS_df[ni,nj] = np.sqrt(np.multiply(WE,WE) + np.multiply(NS,NS))
            
            ang[ni,nj] = atan2(NS,WE)
   
    """  Temporal Gradient """ 
    years = np.arange(0,ntime,12)
    nyears = len(years)
    var_years = np.zeros((nyears,nlat,nlon))
    i = 0
    for yy in years :
        temp = var[yy:yy+12,:,:].mean(axis=0)
        var_years[i,:,:] = temp
        i=i+1
    grad_Temp = np.zeros((nlat,nlon))
    # print('gradT')

    for ni in range (0,nlat):
        for nj in range (0,nlon):
            s,_,_,_,_ = linregress(np.arange(0,nyears),var_years[:,ni,nj])
            grad_Temp[ni,nj] = s    
      
    climvel = np.divide(grad_Temp,gradS_df)
    grad_Temp[grad_Temp==0] = np.nan
    
    """ append variables to ensemble list """
    cv.append(climvel)
    gradS.append(gradS_df)
    gradT.append(grad_Temp)
    gradNS.append(gradS_NS)
    gradWE.append(gradS_WE)
    angles.append(ang)


""" ensemble mean + std """
for ii in range(len(cv)):
    cv[ii] = cv[ii][1:-1,1:-3]
cv_fin = np.ma.array(cv)
cv_std = np.nanstd(cv_fin,axis=0)
cv_mean = np.ma.mean(cv_fin,axis=0)


lat = lat[1:-1,1:-3]
lon = lon[1:-1,1:-3]

   

""" plot """
ratio = cv_std/cv_mean*100

# from matplotlib.colors import BoundaryNorm, ListedColormap

# cmap = mpl.cm.get_cmap('bwr',13)    # PiYG
# my_cmap = []
# for i in range(cmap.N):
#     my_cmap.append(cmap(i))
# cmap = mpl.colors.ListedColormap(my_cmap)
# levels = [-100,-50,-20,-10,-5,0,5,10,20,50,100]#np.linspace(vmin,vmax,9)
# ticks = [-100,-50,-20,-10,-5,-1,0,1,5,10,20,50,100]
# bounds = [-100,-50,-20,-10,-5,-1,0,1,5,10,20,50,100]
# norm =  BoundaryNorm(bounds, cmap.N)
    
    
vmin = 0; vmax=100
levels = np.linspace(vmin,vmax,11)
ticks = np.linspace(vmin,vmax,6)
cmap='OrRd'
fig,ax = plt.subplots(1,1,subplot_kw={'projection':ccrs.PlateCarree()})
plt.title('std cv phi / cv phi')
m = ax.contourf(lon,lat,np.abs(ratio),vmin=vmin,vmax=vmax,levels=levels,extend='max',cmap=cmap)#,norm=norm)
cbar = plt.colorbar(m,ax=ax,ticks=ticks,spacing='uniform')
ax.coastlines()
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=1, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
cbar.ax.set_ylabel('%',fontsize=14)        
plt.show()



my_cmap = []
cmap = mpl.cm.get_cmap('YlOrRd', 7)    # PiYG
for i in range(cmap.N):
    my_cmap.append(cmap(i))
cmap = mpl.colors.ListedColormap(my_cmap)
from matplotlib.colors import BoundaryNorm, ListedColormap
vmin = 0; vmax = 100

levels = [0,0.5,1,2,5,10,50,100]#np.linspace(vmin,vmax,9)
ticks = [0,0.5,1,2,5,10,50,100]#np.linspace(vmin,vmax,5)
bounds = [0,0.5,1,2,5,10,50,100]
norm =  BoundaryNorm(bounds, cmap.N)

# vmin = 0; vmax=100
# levels = np.linspace(vmin,vmax,11)
# ticks = np.linspace(vmin,vmax,6)
# cmap='hot'
fig,ax = plt.subplots(1,1,subplot_kw={'projection':ccrs.PlateCarree()})
plt.title('std cv phi')
m = ax.contourf(lon,lat,cv_std,vmin=vmin,vmax=vmax,levels=levels,extend='max',cmap=cmap,norm=norm)
cbar = plt.colorbar(m,ax=ax,ticks=ticks,spacing='uniform')
ax.coastlines()
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=1, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
cbar.ax.set_ylabel('km/yr',fontsize=14)        
plt.show()


# path = '/home/p/parouffe/Documents/article/fig_article/'
# plt.savefig(path+'std_cv_phi_sur_phi_200m.png')

# """ lpot std """

# stds = np.ma.array(stds)
# ll = [0,0,1,1,2,2]
# cc = [0,1,0,1,0,1]

# my_cmap = []
# cmap = mpl.cm.get_cmap('GnBu', 7)    # PiYG
# for i in range(cmap.N):
#     my_cmap.append(cmap(i))
# cmap = mpl.colors.ListedColormap(my_cmap)
# from matplotlib.colors import BoundaryNorm, ListedColormap

# fig1, ax1 = plt.subplots(3,2,figsize=(7,9),subplot_kw={'projection':ccrs.PlateCarree()})
# # cmap = mpl.colors.ListedColormap(['grey','red','yellow','orange','forestgreen','cyan','grey','brown'])
# # norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)  
# for ni in range(6):
#     vmin = 0; vmax = 100

#     levels = [0,0.5,1,2,5,10,50,100]#np.linspace(vmin,vmax,9)
#     ticks =  [0,0.5,1,2,5,10,50,100]#np.linspace(vmin,vmax,5)
#     bounds = [0,0.5,1,2,5,10,50,100]
#     norm =  BoundaryNorm(bounds, cmap.N)
    
#     ax1[ll[ni],cc[ni]].set_title(str(depths[ni])+'m')
#     m = ax1[ll[ni],cc[ni]].contourf(lon,lat,stds[ni],cmap=cmap,extend='max',
#                                     vmin=vmin,vmax=vmax,levels=levels,norm=norm)
#     ax1[ll[ni],cc[ni]].coastlines()
#     # cbar = plt.colorbar(m,ax=ax1[ll[ni],cc[ni]],shrink=0.8,ticks=ticks,spacing='uniform')
#     # if cc[ni]==1:
#     #     cbar.ax.set_ylabel('Ecart-type',fontsize=12,labelpad = 10)
#     gl = ax1[ll[ni],cc[ni]].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                   linewidth=1, color='gray', alpha=1, linestyle='--')
#     if ll[ni] != 2 and cc[ni] ==1:
#         gl.xlabels_top = False
#         gl.xlabels_bottom = False
#         gl.ylabels_right = False
#         gl.ylabels_left = False
#     if ll[ni]==2 and cc[ni] ==1 :
#         gl.xlabels_top = False
#         gl.ylabels_right = False
#         gl.ylabels_left = False
#     if ll[ni]==2 and cc[ni]==0:
#         gl.xlabels_top = False
#         gl.ylabels_right = False
#     if ll[ni]!=2 and cc[ni]==0:
#         gl.xlabels_top = False
#         gl.xlabels_bottom = False
#         gl.ylabels_right = False

# cb_ax = fig1.add_axes([0.85, 0.25, 0.02, 0.5])
# cbar = plt.colorbar(m,cax=cb_ax,shrink=0.8,ticks=ticks,spacing='uniform')
# cbar.ax.set_ylabel('Ecart-type',fontsize=12,labelpad = 10)


# fig1.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.85,
#                     wspace=0.06, hspace=0.15)  
# path = '/home/p/parouffe/Documents/figures/finales/'
# plt.savefig(path+'stds_PHI.png')


