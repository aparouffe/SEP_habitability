#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 09:37:04 2021

@author: parouffe
"""

from netCDF4 import Dataset
import numpy as np             # version 1.19.2
from scipy.stats import linregress
import datetime
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from math import *
import glob
from matplotlib.colors import DivergingNorm


var_name = "O2"
# depths = [0,100, 300,500,800]    # 1 to 60
depths = [200]    # 1 to 60
lon_ref = 110
zone = 'omz'  # 'NE', 'E', 'full','large'
save = 'no' # 'yes' if yes

path_fig = '/home/p/parouffe/Documents/figures/3.1./cv_O2_bis/'

files = glob.glob('/data/sysco2/travail_en_cours/parouffe/data/HUMBOLT/rcp/O2/*.nc')
files.sort()

""" config of graphs displaying all depths """
# fig, ax = plt.subplots(6,3,subplot_kw={'projection':ccrs.PlateCarree()})

stds = []  # dim[1:-1,1:-3]    si = 0 #subplot counter
cvs = []
gS = []
gT = []

angss = []
ulats =[]; ulons = []
uns = []; vns = []

si = 0 #subplot counter

for depth in depths:
    
    """ initialisation of variables """
    
    cv = []
    gradS = []; gradT = []
    gradNS = []; gradWE = []
    angles = []
    

    for filename in files[:] :
        print(filename)
        nc = Dataset(filename)
        
        # depth level to extract
        zt = nc.variables["z_t"][:]
        zt = zt/100
        depth_lvl = np.abs(zt-depth).argmin()
        
        # variable extraction  
        dxt = nc.variables["DXT"][:]
        dyt = nc.variables["DYT"][:]    
        
        dxt = dxt/10**5    # to km
        dyt = dyt/10**5    # to km
        
        # extraction of variables
        time = nc.variables['time'][:]  # extract days since...
        lat = nc.variables['TLAT'][:]   # lat of T points
        lon = nc.variables['TLONG'][:]  # lon of T points
        var = nc.variables[var_name][:,depth_lvl,:,:]  # T,O2 at T-points
                                               # dim [time,z,lat=116,lon=97]
        ntime, nlat, nlon = np.shape(var)
    
        """ Spatial Gradient """
        var_tmean = np.ma.mean(var,axis=0)  # mean on time axis
        
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
          
        grad_Temp[grad_Temp==0] = np.nan
    
        """ CV """
        
        climvel = np.divide(grad_Temp,gradS_df)
        
        """ append variables to ensemble list """
        cv.append(climvel)
        gradS.append(gradS_df)
        gradT.append(grad_Temp)
        gradNS.append(gradS_NS)
        gradWE.append(gradS_WE)
        angles.append(ang)
    
    
    """ ensemble mean + std """
    cv = np.array(cv)
    cv_std = np.std(cv,axis=0)
    gradS = np.array(gradS); gradS_mean = gradS.mean(axis=0)
    gradNS = np.array(gradNS); gradNS_mean = gradNS.mean(axis=0)
    gradWE = np.array(gradWE); gradWE_mean = gradWE.mean(axis=0)
    gradT = np.array(gradT); gradT_mean = gradT.mean(axis=0)
    angles = np.array(angles); angles = angles.mean(axis=0)
    cv_mean = cv.mean(axis=0)
    
    """ cut boarders """
    cv_mean = cv_mean[1:-1,1:-3]
    cv_std = cv_std[1:-1,1:-3]
    gradS_mean = gradS_mean[1:-1,1:-3]
    gradNS_mean = gradNS_mean[1:-1,1:-3]
    gradWE_mean = gradWE_mean[1:-1,1:-3]
    gradT_mean = gradT_mean[1:-1,1:-3]
    angles = angles[1:-1,1:-3]
    lat = lat[1:-1,1:-3]
    lon = lon[1:-1,1:-3]
    
    """ configuration of geographical area"""
    
    lon2 = np.abs(lon[0,:]-(360-lon_ref)).argmin()
  
    # selection zone
    lat = lat[:,lon2:]
    lon = lon[:,lon2:]
    u = gradWE_mean[:,lon2:].reshape(-1)
    v = gradNS_mean[:,lon2:].reshape(-1)
    u_norm = np.divide(u, np.sqrt(u**2. + v**2.))  # normalize for same arrow size
    v_norm = np.divide(v, np.sqrt(u**2. + v**2.))
    ulon = lon.reshape(-1)
    ulat = lat.reshape(-1) 
    
    liste_angles = angles[:,lon2:-1].reshape(-1)
    gradS_mean = gradS_mean[:,lon2:]
    gradT_mean = gradT_mean[:,lon2:]
    cv_mean = cv_mean[:,lon2:]
    cv_std = cv_std[:,lon2:]

    
    cvs.append(cv_mean)
    gS.append(gradS_mean)
    gT.append(gradT_mean)
    stds.append(cv_std) 

    angss.append(liste_angles)
    ulats.append(ulat)
    ulons.append(ulon)
    uns.append(u_norm)
    vns.append(v_norm)
    
    """ max min cv et std """    
    print(depth)
    print(np.nanmax(cv_mean))
    print(np.nanmin(cv_mean))
        

        
""" PLOT CV + gradS + gradT a une profondeur"""
fig, ax = plt.subplots(1,1,figsize=(4,4),subplot_kw={'projection':ccrs.PlateCarree()})
plt.subplots_adjust(top=0.95,left=0.04,right=0.99,bottom=0.05)

# plt.suptitle(str(depths[ni])+'m')
# CV
vmin = -10; vmax = 10
levels = np.linspace(vmin,vmax,11)# np.arange(vmin,vmax+1,2.5)
ticks = np.linspace(vmin,vmax,6) #np.arange(vmin,vmax+1,5)
norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
cmap = 'RdYlBu_r'
m = ax.contourf(lon,lat,cvs[ni],\
                    vmin=vmin,vmax=vmax,levels=levels,\
                        cmap=cmap,extend='both',norm=norm)
cbar = plt.colorbar(m,pad=0.12,ax=ax,ticks=ticks,orientation='horizontal',shrink=0.8)
cbar.ax.set_xlabel('km/yr',fontsize=12)#.set_rotation(0)
cbar.ax.tick_params(labelsize=11)
cbar.ax.set_xticklabels(ticks)#,rotation = 45)
ax.coastlines()
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=1, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
    
    
for ni in range(len(depths)):
    """ PLOT CV + gradS + gradT """
    fig, ax = plt.subplots(1,3,figsize=(8,4),subplot_kw={'projection':ccrs.PlateCarree()})
    plt.subplots_adjust(top=0.95,left=0.1,right=0.95)

    # plt.suptitle(str(depths[ni])+'m')
    # CV
    vmin = -10; vmax = 10
    levels = np.linspace(vmin,vmax,11)# np.arange(vmin,vmax+1,2.5)
    ticks = np.linspace(vmin,vmax,6) #np.arange(vmin,vmax+1,5)
    # norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
    cmap = 'RdYlBu_r'
    # vmin= -25;vmax=25
    # levels = [-20,-15,-10,-5,-2,-1,0,1,2,5,10,15,20]
    # ticks = [-20,-15,-10,-5,-2,-1,0,1,2,5,10,15,20]#[-20,-10,-5,-1,0,1,5,10,20]
    # bounds = [-20,-15,-10,-5,-2,-1,0,1,2,5,10,15,20]
    norm =  DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
    m = ax[0].contourf(lon,lat,cvs[ni],\
                        vmin=vmin,vmax=vmax,levels=levels,\
                            cmap=cmap,extend='both',norm=norm)
    cbar = plt.colorbar(m,pad=0.12,ax=ax[0],ticks=ticks,orientation='horizontal')
    cbar.ax.set_xlabel('km/an',fontsize=12)#.set_rotation(0)
    cbar.ax.tick_params(labelsize=11)
    cbar.ax.set_xticklabels(ticks)#,rotation = 45)
    ax[0].coastlines()
    gl = ax[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=1, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    
    # gradS
    vmin = 0; vmax = 1
    levels = np.linspace(vmin,vmax,11)
    norm = DivergingNorm(vmin=vmin,vcenter=0.001,vmax=vmax)
    ticks1 = np.round(np.linspace(vmin,vmax,6),1)
    cmap = 'coolwarm'

    m1 = ax[1].contourf(lon,lat, gS[ni],\
                        vmin=vmin,vmax=vmax,levels=levels,\
                            cmap=cmap,extend='max',norm=norm)
    ax[1].quiver(ulons[ni][::10],ulats[ni][::10],uns[ni][::10],vns[ni][::10],transform = ccrs.PlateCarree(),angles ="uv",\
                  width=0.002)
    cbar = plt.colorbar(m1,pad=0.12,ax=ax[1],orientation='horizontal',ticks=ticks1)
    cbar.ax.set_xlabel('mmol/m'+'\u00B3 / km',fontsize=12)#.set_rotation(0)
    cbar.ax.set_xticklabels(ticks1)#,rotation = 45)
    cbar.ax.tick_params(labelsize=11)
    ax[1].coastlines()
    gl = ax[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=1, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylabels_left = False

    # gradT
    vmin = -0.5; vmax=0.5
    levels = np.linspace(vmin,vmax,11)  
    ticks2 = np.round(np.linspace(vmin,vmax,6) ,1)
    cmap='coolwarm'    
    m1 = ax[2].contourf(lon,lat, gT[ni],\
                        vmin=vmin,vmax=vmax,levels=levels,\
                            cmap=cmap,extend='both')
    cbar = plt.colorbar(m1,pad=0.12,ax=ax[2],orientation='horizontal',ticks=ticks2)
    cbar.ax.set_xlabel('mmol/m'+'\u00B3 / an',fontsize=12)#.set_rotation(0)
    cbar.ax.set_xticklabels(ticks2)#,rotation = 45)

    cbar.ax.tick_params(labelsize=11)
    ax[2].coastlines()
    gl = ax[2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=1, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylabels_left = False

        
    path_fig = '/home/p/parouffe/Documents/figures/finales/'
    plt.savefig(path_fig+'CV_O2_'+str(depths[ni])+'m.png')  

""" tout CV + grad """
# cvs = np.ma.array(cvs)
# gS = np.ma.array(gS)
# gT = np.ma.array(gT)

# ll = [0,1,2,3,4,5]#0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5]
# cc = [0,1,2]#,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2]

# fig1, ax1 = plt.subplots(5,3,figsize=(7,9.5),subplot_kw={'projection':ccrs.PlateCarree()})
# plt.gcf().subplots_adjust(bottom=0.1)
# for ni in range(5):
#     ax1[ll[ni],0].set_ylabel(str(depths[ni])+'m')

#     for nj in range (3):
#         if nj == 0:
#             vmin = -10; vmax = 10
#             levels = np.linspace(vmin,vmax,11)
#             ticks = np.linspace(vmin,vmax,11)
#             norm =  DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
            
#             m = ax1[ll[ni],0].contourf(lon,lat,cvs[ni],cmap='RdYlBu_r',extend='both',
#                                             vmin=vmin,vmax=vmax,levels=levels,norm=norm)
#             ax1[ll[ni],0].coastlines()
#             titre = [str(depths[ni])+'m']
#             print(titre)
#             ax1[ni,0].set_ylabel(str(depths[ni])+'m')
#             if ni ==4 :
#                 cbar_ax = plt.gcf().add_axes([0.16, 0.072, 0.21, 0.015])
#                 cbar = plt.colorbar(m,pad=0.12,ax=ax1[ni,0],orientation='horizontal',cax=cbar_ax)
#                 cbar.ax.set_xlabel('km / an',fontsize=12)#.set_rotation(0)
#             gl = ax1[ll[ni],0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                           linewidth=1, color='gray', alpha=1, linestyle='--')
#         if nj == 1:
#             # gradS
#             ax1[ni,1].set_title(str(depths[ni])+'m',fontsize=10,pad=-1)
#             vmin = 0; vmax = 1
#             levels = np.linspace(vmin,vmax,11)
#             ticks1 = np.round(np.linspace(vmin,vmax,6),1)
#             cmap = 'coolwarm'
#             norm =  DivergingNorm(vmin=vmin,vcenter=0.001,vmax=vmax)
        
#             m1 = ax1[ni,1].contourf(lon,lat, gS[ni],\
#                                 vmin=vmin,vmax=vmax,levels=levels,\
#                                     cmap=cmap,extend='max',norm=norm)
#             # ax1[ni,1].quiver(ulons[ni][::10],ulats[ni][::10],uns[ni][::10],vns[ni][::10],transform = ccrs.PlateCarree(),angles ="uv",\
#             #               width=0.002)
#             if ni == 4:
#                 cbar_ax = plt.gcf().add_axes([0.42, 0.072, 0.21, 0.015])
#                 cbar = plt.colorbar(m1,pad=0.12,ax=ax1[ni,1],orientation='horizontal',ticks=ticks1,cax=cbar_ax)
#                 cbar.ax.set_xlabel('mmol/m'+'\u00B3 / km',fontsize=11)#.set_rotation(0)
#                 cbar.ax.set_xticklabels(ticks1,rotation = 45)
#             # cbar.ax.tick_params(labelsize=10)
#             ax1[ni,1].coastlines()
#             gl = ax1[ni,1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                               linewidth=1, color='gray', alpha=1, linestyle='--')
#         if nj ==2:
#             # gradT
#             vmin = -0.4; vmax=0.4
#             levels = np.linspace(vmin,vmax,9)  
#             ticks2 = np.round(np.linspace(vmin,vmax,5),1)
#             cmap='coolwarm'    
#             m1 = ax1[ni,2].contourf(lon,lat, gT[ni],\
#                                 vmin=vmin,vmax=vmax,levels=levels,\
#                                     cmap=cmap,extend='both')
#             if ni == 4:
#                 cbar_ax = plt.gcf().add_axes([0.67, 0.072, 0.21, 0.015])
#                 cbar = plt.colorbar(m1,pad=0.12,ax=ax1[ni,2],orientation='horizontal',ticks=ticks2,cax=cbar_ax)
#                 cbar.ax.set_xlabel('mmol/m'+'\u00B3 / an',fontsize=11)#.set_rotation(0)
#                 cbar.ax.set_xticklabels(ticks2,rotation = 45)
#                 # cbar.ax.tick_params(labelsize=12)
#             ax1[ni,2].coastlines()
#             gl = ax1[ni,2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                           linewidth=1, color='gray', alpha=1, linestyle='--')

#         if nj == 1 or nj == 2:
#             gl.ylabels_right = False
#             gl.ylabels_left = False
#             gl.xlabels_top = False
#             if ni != 4:
#                 gl.xlabels_bottom = False
#         if nj == 0:
#             gl.ylabels_right = False
#             gl.xlabels_top = False
#             if ni != 4:
#                 gl.xlabels_bottom = False

# fig1.subplots_adjust(bottom=0.12, top=0.98, left=0.15, right=0.90,
#                     wspace=0.095, hspace=0.095)          
 
# path = '/home/p/parouffe/Documents/figures/3.1./finales/'
# plt.savefig(path+'CV_grad_O2.png')
            
      # ax1[ni,nj].set_aspect('auto')

# cb_ax = fig1.add_axes([0.85, 0.25, 0.02, 0.5])
# cbar = plt.colorbar(m,cax=cb_ax,shrink=0.8,ticks=ticks,spacing='uniform')
# cbar.ax.set_ylabel('vitesses climatiques (km/an)',labelpad = 10,fontsize=14)


# fig1.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.85,
                     
# """ PLOT CV UNIQUEMENT 3x2"""

# path_fig = '/home/p/parouffe/Documents/figures/3.1./finales/'
# cvs = np.ma.array(cvs)
# ll = [0,0,1,1,2,2]
# cc = [0,1,0,1,0,1]


# # my_cmap = []
# # cmap = mpl.cm.get_cmap('', 7)    # PiYG
# # for i in range(cmap.N):
# #     my_cmap.append(cmap(i))
# # cmap = mpl.colors.ListedColormap(my_cmap)
# # from matplotlib.colors import BoundaryNorm, ListedColormap
# fig1, ax1 = plt.subplots(3,2,figsize=(7,9),subplot_kw={'projection':ccrs.PlateCarree()})
# # cmap = mpl.colors.ListedColormap(['grey','red','yellow','orange','forestgreen','cyan','grey','brown'])
# # norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)  
# for ni in range(6):
#     vmin = -10; vmax = 10

#     levels = np.linspace(-10,10,11)
#     ticks = np.linspace(-10,10,11)
#     bounds = np.linspace(-10,10,11)
#     norm =  DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
    
#     ax1[ll[ni],cc[ni]].set_title(str(depths[ni])+'m')
#     m = ax1[ll[ni],cc[ni]].contourf(lon,lat,cvs[ni],cmap='RdYlBu_r',extend='both',
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
# cbar.ax.set_ylabel('vitesses climatiques (km/an)',labelpad = 10,fontsize=14)


# fig1.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.85,
#                     wspace=0.06, hspace=0.15)  
# plt.savefig(path_fig+'CV_O2.png')

# for ni in range(6):
#     print(np.nanmax(stds[ni]))





# fig, ax = plt.subplots(6,3,figsize=(7,10),subplot_kw={'projection':ccrs.PlateCarree()},constrained_layout=True)
# for si in range(6):    
#     # CV
#     ax[si,1].set_xlabel('bouh')#str(depths[si])+'m'
#     vmin = -10; vmax = 10
#     levels = np.linspace(vmin,vmax,11)# np.arange(vmin,vmax+1,2.5)
#     ticks = np.linspace(vmin,vmax,6) #np.arange(vmin,vmax+1,5)
#     norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
#     cmap = 'RdYlBu_r'
    
#     m = ax[si,0].contourf(lon,lat,cvs[si],\
#                         vmin=vmin,vmax=vmax,norm=norm,levels=levels,\
#                             cmap=cmap,extend='both')
#     ax[si,0].coastlines()
#     gl = ax[si,0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                       linewidth=1, color='gray', alpha=1, linestyle='--') 
#     if si==5:
#         cbar = plt.colorbar(m,ax=ax[si,0],ticks=ticks,orientation='horizontal')
#         cbar.ax.set_xlabel('km/an',fontsize=10)#.set_rotation(0)
#         cbar.ax.tick_params(labelsize=10)
#         cbar.ax.set_xticklabels(ticks,rotation = 45,fontsize = 10)

#     gl.xlabels_top = False
#     gl.ylabels_right = False
#     gl.xlabels_bottom = False
#     gl.ylabels_left = False
    
    
#     # gradS
#     vmin = 0; vmax = 1
#     levels = np.linspace(vmin,vmax,11)
#     # norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
#     ticks = np.round(np.linspace(vmin,vmax,6),2)
#     cmap = 'Blues'
#     m1 = ax[si,1].contourf(lon,lat, gS[si],\
#                         vmin=vmin,vmax=vmax,levels=levels,\
#                             cmap=cmap,extend='both')
#     ax[si,1].coastlines()
#     gl = ax[si,1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                       linewidth=1, color='gray', alpha=1, linestyle='--')
#     if si==5:
#         # cbar = plt.colorbar(m,cax=cb_ax, ticks=ticks)
#         cbar = plt.colorbar(m1,ax=ax[si,1],ticks=ticks,location='bottom')
#         cbar.ax.set_xlabel('km/an',fontsize=10)#.set_rotation(0)
#         cbar.ax.tick_params(labelsize=10)
#         cbar.ax.set_xticklabels(ticks,rotation = 45,fontsize = 10)

#     gl.xlabels_top = False
#     gl.ylabels_right = False
#     gl.xlabels_bottom = False
#     gl.ylabels_left = False
#     #gradT
#     vmin=-0.5; vmax= 0.5
#     levels = np.linspace(vmin,vmax,11)
#     tickss = np.round(np.linspace(vmin,vmax,6),2)
#     cmap = 'RdBu_r'

#     m2 = ax[si,2].contourf(lon,lat, gT[si],\
#                         vmin=vmin,vmax=vmax,levels=levels,\
#                             cmap=cmap,extend='both')
#     ax[si,2].coastlines()
#     gl = ax[si,2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                       linewidth=1, color='gray', alpha=1, linestyle='--')
#     if si==5:
#         cbar = plt.colorbar(m2,ax=ax[si,2],ticks=tickss,orientation='horizontal')
#         cbar.ax.set_xlabel('km/an',fontsize=10)#.set_rotation(0)
#         cbar.ax.tick_params(labelsize=10)
#         cbar.ax.set_xticklabels(tickss,rotation = 45,fontsize = 10)

#     gl.xlabels_top = False
#     gl.ylabels_right = False
#     gl.xlabels_bottom = False
#     gl.ylabels_left = False
        
#     # ax[si,0].set_aspect('auto')
#     # ax[si,1].set_aspect('auto')
#     # ax[si,2].set_aspect('auto')
# ax[0,0].set_title('CV',fontsize=10)
# ax[0,1].set_title('gradient spatial',fontsize=10)
# ax[0,2].set_title('gradient temporel',fontsize=10)
# plt.subplots_adjust(left=0.1, right=0.9, bottom=0.14, top=0.97,hspace=0.04,wspace=0.04)

      
    # plt.savefig(path_fig+'CV_'+str(depth)+'m.png')
    
    
# """ lpot std """
# path_fig = '/home/p/parouffe/Documents/figures/3.1./'

# stds = np.ma.array(stds)
# ll = [0,0,1,1,2,2]
# cc = [0,1,0,1,0,1]


# fig1, ax1 = plt.subplots(3,2,figsize=(7,9),subplot_kw={'projection':ccrs.PlateCarree()})

# for ni in range(6):
#     vmin = 0; vmax = 2.5
#     levels = np.linspace(vmin,vmax,11)
#     ticks = np.linspace(vmin,vmax,6)
#     ax1[ll[ni],cc[ni]].set_title(str(depths[ni])+'m')
#     m = ax1[ll[ni],cc[ni]].contourf(lon,lat,stds[ni],cmap='GnBu',extend='max',\
#                                     vmin=vmin,vmax=vmax,levels = levels)
#     # cbar = plt.colorbar(m,ax=ax1[ll[ni],cc[ni]],ticks=ticks)
#     ax1[ll[ni],cc[ni]].coastlines()
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
# cbar = plt.colorbar(m,cax=cb_ax, ticks=ticks)
# cbar.ax.set_ylabel('Ecart-type',fontsize=12,labelpad = 10)

# fig1.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
#                     wspace=0.13, hspace=0.13)    # # bornes = np.arange(-15,345,30)
    
# # plt.savefig(path_fig+'stds.png')

# for ni in range(6):
#     print(np.nanmax(stds[ni]))
    # """ configuration of rosewind """
    
    # cat_E = []    # 345 - 15  
    # cat_EEN = []  # 15-45
    # cat_ENN = []  # 45-75
    # cat_N = []    # 75-105
    # cat_NNW = []  # 105-135
    # cat_NWW = []  # 135-165
    # cat_W = []    # 165-195
    # cat_WWS = []  # 195-225
    # cat_WSS = []  # 225-255
    # cat_S = []    # 255-285
    # cat_SSE = []  # 285-315
    # cat_SEE = []  # 315-345
    
    # p = []
    # borne = np.linspace(-np.pi,np.pi,13)
    # borne+=np.pi/12
    # deg = np.rad2deg(borne)
    # liste = [0]
    # for ang in liste_angles :
    #     if (ang>borne[5]) and (ang<=borne[6]):
    #         cat_E.append(ang)
    #     elif (ang>borne[6]) and (ang<=borne[7]):
    #         cat_EEN.append(ang)
    #     elif (ang>borne[7]) and (ang<=borne[8]):
    #         cat_ENN.append(ang)
    #     elif (ang>borne[8]) and (ang<=borne[9]):
    #         cat_N.append(ang)
    #     elif (ang>borne[9]) and (ang<=borne[10]):
    #         cat_NNW.append(ang)
    #     elif (ang>borne[10]) and (ang<=borne[11]):
    #         cat_NWW.append(ang)
    #     elif (ang<borne[0]) or (ang>=borne[11]):
    #         cat_W.append(ang)
    #     elif (ang>borne[0]) and (ang<=borne[1]):
    #         cat_WWS.append(ang)
    #     elif (ang>borne[1]) and (ang<=borne[2]):
    #         cat_WSS.append(ang)
    #     elif (ang>borne[2]) and (ang<=borne[3]):
    #         cat_S.append(ang)    
    #     elif (ang>borne[3]) and (ang<=borne[4]):
    #         cat_SSE.append(ang)                     
    #     elif (ang>borne[4]) and (ang<=borne[5]):
    #         cat_SEE.append(ang)                      
    #     else :
    #         p.append(ang)
    
    # freq = [len(cat_E),len(cat_EEN),len(cat_ENN),len(cat_N),len(cat_NNW),len(cat_NWW),\
    #         len(cat_W),len(cat_WWS),len(cat_WSS),len(cat_S),len(cat_SSE),len(cat_SEE)]
    # freq = np.array(freq)/len(liste_angles)*100                    
    
    # n = len(cat_E)+len(cat_EEN)+len(cat_ENN)+len(cat_N)+len(cat_NNW)+len(cat_NWW)+\
    #         len(cat_W)+len(cat_WWS)+len(cat_WSS)+len(cat_S)+len(cat_SSE)+len(cat_SEE)    
    
    # cat = np.linspace(0,2*np.pi,13)
    
    
    # # if zone == 'full':
    # #         fig = plt.figure(figsize=(6,6))
    # #         ax = fig.add_subplot(1,1,1, projection='polar')
        
    # #         colors = "seagreen"
    # #         bars = ax.bar(cat[:-1],freq,color=colors,width=0.5)
    # #         ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'])
    # #         ax.xaxis.grid('major')
    # #         plt.subplots_adjust(left=0.09, bottom=0.08, right=0.98, top=0.7 , hspace=0.1 )
    # #         plt.title("direction of climate velocities - Temperature \n "+str(depth)+'m \n \n occurences in % \n')         
    
    # # else :
    
    # if zone == 'omz' or zone =='large':
    #     fig = plt.figure(figsize=(10,6.5))
    #     if zone == 'omz':
    #         shrink=0.6
    #     else :
    #         shrink = 0.8
    # if zone == 'coast':
    #     fig = plt.figure(figsize=(9,7.5))
    #     shrink = 0.9
        
    # plt.suptitle('O2 Climate Velocity - RCP8.5 \n '+str(depth)+'m', fontsize=15)
    # ax = fig.add_subplot(121,projection=ccrs.PlateCarree())
    # vmin = -20; vmax = 10
    
    # levels = np.linspace(vmin,vmax,16)
    # norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
    # bar_ticks = np.linspace(vmin,vmax,7)
    # m = plt.contourf(lon,lat,cv_mean,vmin=vmin,vmax=vmax,\
    #                   levels=levels,extend='both', cmap='RdBu_r',\
    #                       norm=norm)
    # crs = ccrs.PlateCarree()
    # ax.quiver(ulon[::10],ulat[::10],u_norm[::10],v_norm[::10],transform = crs,angles ="uv",\
    #               width=0.002)
    # cbar = plt.colorbar(m,ticks=bar_ticks,pad=0.10,shrink=shrink)
    # cbar.ax.set_ylabel("km/yr",fontsize=15,labelpad=5)
    # cbar.ax.tick_params(labelsize=12)
    # ax.coastlines()
    # gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    #                   linewidth=1, color='gray', alpha=1, linestyle='--')
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    # # plt.subplots_adjust(left = 0.090,right=0.7,bottom=0.06, top = 0.9,hspace=0.9 )
    # ax.set_title('ensemble mean',pad=15)
    
    # vmin=0;vmax=10;
    # levels=np.linspace(0,vmax,11)
    # bar_ticks=np.linspace(0,vmax,6)
    # ax = fig.add_subplot(122,projection=ccrs.PlateCarree())
    # m = plt.contourf(lon,lat,cv_std,cmap='Reds',\
    #                      vmin=vmin,vmax=vmax,levels=levels,extend='max')
    # cbar = plt.colorbar(m,pad=0.10,shrink=shrink)
    # cbar.ax.set_ylabel("km/yr",fontsize=15,labelpad=10)
    # cbar.ax.tick_params(labelsize=12)
    # ax.coastlines()
    # gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    #                   linewidth=1, color='gray', alpha=1, linestyle='--')
    # # plt.subplots_adjust(left = 0.070,right=0.99,bottom=0.06, top = 0.9,hspace=0.9 )
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    # gl.ylabels_left = False
    # ax.set_title('standard deviation',pad=15)
    
    # if save=='yes':
    #     plt.savefig(path_fig+var_name+'_'+zone+'_'+str(depth)+'m.png')
    
    # # """fig rosewind """
    # # fig1 = plt.figure(figsize=(5,5))
    # # ax = fig1.add_subplot(111, projection='polar')
    # # if zone == 'coast':
    # #     color = "royalblue"

    # # bars = ax.bar(cat[:-1],freq,color=color,width=0.5)
    # # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'])
    # # ax.xaxis.grid('major')
    # # plt.subplots_adjust(left=0.09, bottom=0.08, right=0.98, top=0.7 , hspace=0.1 )
    # # ax.set_title(zone+'  '+str(depth)+'m \noccurences in % \n',fontsize=14,pad=2)                                      
    # # plt.subplots_adjust(left = 0.05,right = 0.9  )
    
    # # plt.savefig(path_fig+var_name+'_ROSEWIND_'+zone+'_'+str(depth)+'m.png')
                         
    
    
    # """  fig gradients """
    
    # if zone == 'coast':
    #     fig2 = plt.figure(figsize=(9,7.5))
    #     plt.suptitle(str(depth)+'m',fontsize=15)
    #     shrink = 0.90
    
    # ### gradS
    # ax = fig2.add_subplot(121,projection=ccrs.PlateCarree())
    # vmin = 0; vmax = 1
    # levels = np.linspace(vmin,vmax,11)
    # # norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
    # ticks = np.linspace(vmin,vmax,6)
    # m = plt.contourf(lon,lat,gradS_mean,vmin=vmin,vmax=vmax,levels=levels,\
    #                   extend='max', cmap='YlGnBu')
    #     #
    # crs = ccrs.PlateCarree()
    # ax.quiver(ulon[::5],ulat[::5],u_norm[::5],v_norm[::5],transform = crs,angles ="uv",\
    #               width=0.002)
    # cbar = plt.colorbar(m,pad=0.10,shrink=shrink,ticks=ticks)
    # # ticks=bar_ticks,
    # cbar.ax.set_ylabel('mmol/m'+'\u00B3'+' /km',fontsize=15,labelpad=5)
    # cbar.ax.tick_params(labelsize=12)
    # ax.coastlines()
    # gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    #                   linewidth=1, color='gray', alpha=1, linestyle='--')
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    # # plt.subplots_adjust(left = 0.090,right=0.7,bottom=0.06, top = 0.9,hspace=0.9 )
    # ax.set_title('Spatial gradient',fontsize=14,pad=15)
    
    # ### graT
    # ax = fig2.add_subplot(122,projection=ccrs.PlateCarree())
    # vmin = -0.5; vmax = 1
    # levels = np.linspace(vmin,vmax,11)
    # norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
    # ticks = np.linspace(vmin,vmax,6)
    # m = plt.contourf(lon,lat,gradT_mean,vmin=vmin,vmax=vmax,levels=levels,\
    #                   extend='both', cmap='RdBu_r',norm=norm)
    # crs = ccrs.PlateCarree()
    # ax.quiver(ulon[::5],ulat[::5],u_norm[::5],v_norm[::5],transform = crs,angles ="uv",\
    #               width=0.002)
    # cbar = plt.colorbar(m,pad=0.10,shrink=shrink,ticks=ticks)
    # cbar.ax.set_ylabel('mmol/m'+'\u00B3'+' /yr',fontsize=15,labelpad=5)
    # cbar.ax.tick_params(labelsize=12)
    # ax.coastlines()
    # gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    #                   linewidth=1, color='gray', alpha=1, linestyle='--')
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    # gl.ylabels_left = False
    
    # plt.subplots_adjust(left = 0.090,right=0.92,bottom=0.06, top = 0.9,hspace=0.9 )
    # ax.set_title('Temporal gradient',fontsize=14,pad=15)
    
    # if save=='yes':
    #     print('in')
    #     plt.savefig(path_fig+var_name+'_GRAD_'+zone+'_'+str(depth)+'m.png')