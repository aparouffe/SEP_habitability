#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 15:50:20 2021

@author: parouffe

this script calculates the ensemble mean of :
    PHI
    CV PHI + 1std
    gradS
    gradT
    directions of gradS
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

# depths = [0,100,300,500,800]    # 1 to 60
depths = [200]

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

"""  """

nc = Dataset(fO2[0])
zt = nc.variables["z_t"][:]
zt = zt/100

# time = nc.variables['time'][:]  # extract days since...
lat = nc.variables['TLAT'][:]   # lat of T points
lon = nc.variables['TLONG'][:]  # lon of T points  

lon2 = np.abs(lon[0,:]-(360-110)).argmin()

lat = lat[:,lon2:]
lon = lon[:,lon2:]

dxt = nc.variables["DXT"][:,lon2:]
dyt = nc.variables["DYT"][:,lon2:]    
dxt = dxt/10**5    # to km
dyt = dyt/10**5    # to km



""" initialisation of variables """

cvs = []; gS = []; gT = []
angss = []
ulats =[]; ulons = []
uns = []; vns = []
conc_o2 = []
stds = []


for depth in depths:
    
    print(depth)
    depth_lvl = np.abs(zt-depth).argmin()      

    cv = []
    gradS = []; gradT = []
    gradNS = []; gradWE = []
    angles = []
    conc = []
    
    for fi in range(0,n) :

        """   spatial variables + depth level"""
        print(fi)
        

        nc1 = Dataset(fO2[fi])
        nc2 = Dataset(ftemp[fi])
        nc3 = Dataset(fsalt[fi])
        # print(fO2[fi])
        # print(ftemp[fi])
        # print(fsalt[fi])
        
        o2 = nc1.variables['O2'][:,depth_lvl,:,lon2:]
        Temp = nc2.variables['TEMP'][:,depth_lvl,:,lon2:]
        salt = nc3.variables['SALT'][:,depth_lvl,:,lon2:]                                               # dim [time,z,lat=116,lon=97]

        
        conc.append(np.ma.mean(o2,axis=0))
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
    # print(depth)
    cv = np.array(cv)
    cv_std = np.nanstd(cv,axis=0)
    gradS = np.ma.array(gradS); gradS_mean = np.ma.mean(gradS,axis=0)
    gradNS = np.ma.array(gradNS); gradNS_mean = np.ma.mean(gradNS,axis=0)
    gradWE = np.array(gradWE); gradWE_mean = np.ma.mean(gradWE,axis=0)
    gradT = np.array(gradT); gradT_mean = np.mean(gradT,axis=0)
    angles = np.array(angles); angles = np.mean(angles,axis=0)
    cv_mean = np.ma.mean(cv,axis=0)
    conc = np.ma.mean(conc,axis=0)

    
    """ configuration of geographical area"""

    cv_mean = cv_mean[1:-1,1:-3]
    gradS_mean = gradS_mean[1:-1,1:-3]
    gradT_mean = gradT_mean[1:-1,1:-3]
    gradWE_mean = gradWE_mean[1:-1,1:-3]
    gradNS_mean = gradNS_mean[1:-1,1:-3]
    u = gradWE_mean.reshape(-1)
    v = gradNS_mean.reshape(-1)
    u_norm = np.divide(u, np.sqrt(u**2. + v**2.))  # normalize for same arrow size
    v_norm = np.divide(v, np.sqrt(u**2. + v**2.))
    ulon = lon[1:-1,1:-3].reshape(-1)
    ulat = lat[1:-1,1:-3].reshape(-1) 
    conc = conc[1:-1,1:-3]

    liste_angles = angles[1:-1,1:-3].reshape(-1)


    cvs.append(cv_mean)
    gS.append(gradS_mean)
    gT.append(gradT_mean)   
    stds.append(cv_std)
    conc_o2.append(conc)
    angss.append(liste_angles)
    ulats.append(ulat)
    ulons.append(ulon)
    uns.append(u_norm)
    vns.append(v_norm)
    
lon = lon[1:-1,1:-3]
lat = lat[1:-1,1:-3]

std = np.nanstd(np.ma.array(cvs[0]),axis=0)


""" PLOT CV + gradS + gradT """
# import matplotlib as mpl

# my_cmap = []
# cmap = mpl.cm.get_cmap('RdYlBu_r', 14)    # PiYG
# for i in range(cmap.N):
#     my_cmap.append(cmap(i))
# cmap_cv = mpl.colors.ListedColormap(my_cmap)
# from matplotlib.colors import BoundaryNorm, ListedColormap


# my_cmap = []
# cmap = mpl.cm.get_cmap('YlOrRd', 7)    # PiYG
# for i in range(cmap.N):
#     my_cmap.append(cmap(i))
# cmap = mpl.colors.ListedColormap(my_cmap)
# from matplotlib.colors import BoundaryNorm, ListedColormap
# vmin = 0; vmax = 100

# levels = [0,0.5,1,2,5,10,50,100]#np.linspace(vmin,vmax,9)
# ticks = [0,0.5,1,2,5,10,50,100]#np.linspace(vmin,vmax,5)
# bounds = [0,0.5,1,2,5,10,50,100]
# norm =  BoundaryNorm(bounds, cmap.N)

# vmin = 0; vmax=100
# levels = np.linspace(vmin,vmax,11)
# ticks = np.linspace(vmin,vmax,6)
# cmap='hot'
# fig,ax = plt.subplots(1,1,subplot_kw={'projection':ccrs.PlateCarree()})
# plt.title('std cv phi')
# m = ax.contourf(lon,lat,cv_std,vmin=vmin,vmax=vmax,levels=levels,extend='max',cmap=cmap,norm=norm)
# cbar = plt.colorbar(m,ax=ax,ticks=ticks,spacing='uniform')
# ax.coastlines()
# gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                   linewidth=1, color='gray', alpha=1, linestyle='--')
# gl.xlabels_top = False
# gl.ylabels_right = False
# cbar.ax.set_ylabel('km/yr',fontsize=14)        
# plt.show()
for ni in range(len(depths)):
    fig, ax = plt.subplots(1,1,figsize=(3.5,3.5),subplot_kw={'projection':ccrs.PlateCarree()})
    plt.subplots_adjust(top=0.98,left=0.05,right=1,bottom=0.05)
    vmin = -10; vmax = 10
    levels = np.linspace(vmin,vmax,11)# np.arange(vmin,vmax+1,2.5)
    ticks = np.linspace(vmin,vmax,6) #np.arange(vmin,vmax+1,5)
    norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
    cmap = 'RdYlBu_r'  
    m = ax.contourf(lon,lat,cvs[ni],\
                        vmin=vmin,vmax=vmax,levels=levels,\
                            cmap=cmap,extend='both',norm=norm)
    m1 = plt.contour(lon,lat,conc_o2[ni],levels=[45],colors='k')
    # ax1.clabel(m3,fmt='',colors='r',fontsize=13)

    cbar = plt.colorbar(m,ax=ax,ticks=ticks,shrink=0.7,orientation='horizontal',pad=0.1)
    cbar.ax.set_xlabel('km/yr',fontsize=10)#.set_rotation(0)
    cbar.ax.tick_params(labelsize=10)
    # cbar.ax.set_yticklabels(ticks)#,rotation = 45)
    ax.coastlines()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=1, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    plt.show()

crs=ccrs.PlateCarree()
    
for ni in range(len(depths)):
    """ PLOT CV + gradS + gradT """
    fig, ax = plt.subplots(1,3,figsize=(8,2.7),subplot_kw={'projection':ccrs.PlateCarree()})
    plt.subplots_adjust(top=0.95,left=0.05,right=0.9,wspace=0.28)
    # CV
    vmin = -10; vmax = 10
    levels = np.linspace(vmin,vmax,11)# np.arange(vmin,vmax+1,2.5)
    ticks = np.linspace(vmin,vmax,6) #np.arange(vmin,vmax+1,5)
    norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
    cmap = 'RdYlBu_r'

    m = ax[0].contourf(lon,lat,cvs[ni],\
                        vmin=vmin,vmax=vmax,levels=levels,\
                            cmap=cmap,extend='both',norm=norm)
    m1 = plt.contour(lon,lat,conc_o2[ni],levels=[45],colors='k')
    cbar = plt.colorbar(m,pad=0.05,ax=ax[0],ticks=ticks,shrink=0.7)#,orientation='horizontal')
    cbar.ax.set_ylabel('km/yr',fontsize=10)#.set_rotation(0)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_yticklabels(ticks)#,rotation = 45)
    ax[0].coastlines()
    gl = ax[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=1, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    ax[0].scatter(281,-34,marker='*',color='red',transform=crs)
    ax[0].scatter(281,-25,marker='*',color='red',transform=crs)
    ax[0].scatter(251,-27,marker='*',color='red',transform=crs)      
    # gradS
    vmin = 0; vmax = 0.02
    levels = np.linspace(vmin,vmax,11)
    norm = DivergingNorm(vmin=vmin,vcenter=0.001,vmax=vmax)
    ticks1 = np.linspace(vmin,vmax,6)
    cmap = 'coolwarm'

    m1 = ax[1].contourf(lon,lat, gS[ni],\
                        vmin=vmin,vmax=vmax,levels=levels,\
                            cmap=cmap,extend='max',norm=norm)
    ax[1].quiver(ulons[ni][::10],ulats[ni][::10],uns[ni][::10],vns[ni][::10],transform = ccrs.PlateCarree(),angles ="uv",\
                  width=0.002)
    cbar = plt.colorbar(m1,pad=0.05,ax=ax[1],ticks=ticks1,shrink=0.7)#,orientation='horizontal')
    cbar.ax.set_ylabel('unit / km',fontsize=10)#.set_rotation(0)
    # cbar.ax.set_yticklabels(ticks1)#,rotation = 45)
    cbar.ax.tick_params(labelsize=10)
    ax[1].coastlines()
    gl = ax[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=1, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylabels_left = False

    # gradT
    vmin = -0.01; vmax=0.01
    levels = np.linspace(vmin,vmax,11)  
    ticks2 = np.linspace(vmin,vmax,6) 
    cmap='coolwarm'    
    m1 = ax[2].contourf(lon,lat, gT[ni],\
                        vmin=vmin,vmax=vmax,levels=levels,\
                            cmap=cmap,extend='both')
    cbar = plt.colorbar(m1,pad=0.05,ax=ax[2],ticks=ticks2,shrink=0.7)#,orientation='horizontal')
    cbar.ax.set_ylabel('unit / yr',fontsize=10)#.set_rotation(0)
    # cbar.ax.set_xticklabels(ticks2)#,rotation = 45)

    cbar.ax.tick_params(labelsize=10)
    ax[2].coastlines()
    gl = ax[2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=1, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.ylabels_left = False
    
    ax[0].set_title('climate velocity')
    ax[1].set_title('spatial gradient')
    ax[2].set_title('temporal gradient')
    
    plt.show()
    # path_fig = '/home/p/parouffe/Documents/figures/finales/'
    # plt.savefig(path_fig+'fig5c_CV_PHI_'+str(depths[ni])+'m.png')  
""" tout """
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
#             levels = np.linspace(-10,10,11)
#             ticks = np.linspace(-10,10,11)
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
#             vmin = 0; vmax = 0.02
#             levels = np.linspace(vmin,vmax,11)
#             ticks1 = np.linspace(vmin,vmax,6)
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
#                 cbar.ax.set_xlabel('point / km',fontsize=11)#.set_rotation(0)
#                 cbar.ax.set_xticklabels(ticks1,rotation = 45)
#             # cbar.ax.tick_params(labelsize=10)
#             ax1[ni,1].coastlines()
#             gl = ax1[ni,1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                               linewidth=1, color='gray', alpha=1, linestyle='--')
#         if nj ==2:
#             # gradT
#             vmin = -0.01; vmax=0.01
#             levels = np.linspace(vmin,vmax,11)  
#             ticks2 = np.linspace(vmin,vmax,6)
#             cmap='coolwarm'    
#             m1 = ax1[ni,2].contourf(lon,lat, gT[ni],\
#                                 vmin=vmin,vmax=vmax,levels=levels,\
#                                     cmap=cmap,extend='both')
#             if ni == 4:
#                 cbar_ax = plt.gcf().add_axes([0.67, 0.072, 0.21, 0.015])
#                 cbar = plt.colorbar(m1,pad=0.12,ax=ax1[ni,2],orientation='horizontal',ticks=ticks2,cax=cbar_ax)
#                 cbar.ax.set_xlabel('point / an',fontsize=11)#.set_rotation(0)
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
# plt.savefig(path+'CV_grad_PHI.png')
      

    
# freq, cat = rosewind(liste_angles)
# for ni in range (0,6):
#     freq, cat = rosewind(angss[ni])

#     fig1 = plt.figure(figsize=(5,5))
#     ax = fig1.add_subplot(111, projection='polar')
#     color = 'darkorange'
#     bars = ax.bar(cat[:-1],freq,color=color,width=0.5)
#     ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'])
#     ax.xaxis.grid('major')
#     plt.subplots_adjust(left=0.09, bottom=0.08, right=0.98, top=0.7 , hspace=0.1 )
#     ax.set_title(str(depths[ni])+'m \noccurences in % \n',fontsize=14,pad=2)                                      
#     plt.subplots_adjust(left = 0.05,right = 0.9  )
    
#     fig2 = plt.figure()
#     ax = fig2.add_subplot(121,projection=ccrs.PlateCarree())
#     vmin = -10; vmax = 10
#     levels = np.linspace(vmin,vmax,11)# np.arange(vmin,vmax+1,2.5)
#     ticks = np.linspace(vmin,vmax,6) #np.arange(vmin,vmax+1,5)
#     cmap = 'RdYlBu_r'

#     m = ax.contourf(lon,lat,cvs[ni],\
#                         vmin=vmin,vmax=vmax,levels=levels,\
#                             cmap=cmap,extend='both')
#     ax.quiver(ulons[ni][::10],ulats[ni][::10],uns[ni][::10],vns[ni][::10],transform = ccrs.PlateCarree(),angles ="uv",\
#                   width=0.002)
#     ax.coastlines()
#     gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                       linewidth=1, color='gray', alpha=1, linestyle='--')
#     gl.xlabels_top = False
#     gl.ylabels_right = False
#     # plt.subplots_adjust(left = 0.090,right=0.7,bottom=0.06, top = 0.9,hspace=0.9 )
#     ax.set_title(str(depths[ni])+'m    direction')
    
    
    # plt.savefig(path_fig+str(depth)+'m_ROSEWIND_'+zone+'_'+str(depth)+'.png')


    # elif zone == "large":
    #     u = np.array(gradS_WE[latmin:latmax,lon1:lon2]).reshape(-1)
    #     v = np.array(gradS_NS[latmin:latmax,lon1:lon2]).reshape(-1)
    #     ulon = np.array(lon[latmin:latmax,lon1:lon2]).reshape(-1)
    #     ulat = np.array(lat[latmin:latmax,lon1:lon2]).reshape(-1) 
    #     liste_angles = angle[latmin:latmax,lon1:lon2].reshape(-1)
    #     var = gradS_df[latmin:latmax,lon1:lon2]
    #     lat = lat[latmin:latmax,lon1:lon2]
    #     lon = lon[latmin:latmax,lon1:lon2]   
    # elif zone == 'full':
    #     u = np.array(gradS_WE).reshape(-1)
    #     v = np.array(gradS_NS).reshape(-1)   
    #     ulon = np.array(lon).reshape(-1)
    #     ulat = np.array(lat).reshape(-1) 
    #     liste_angles = angle.reshape(-1)
    # if zone =='full':
    #     var = np.mean(var,axis=0)
                 
    # u_norm = np.divide(u, np.sqrt(u**2. + v**2.))   # normalize for same arrow size
    # v_norm = np.divide(v, np.sqrt(u**2. + v**2.))
    
    
    
    # """ configuration of rosewind """
    
    # freq, cat = rosewind(liste_angles)

    # # if zone == 'full':
    # #         fig = plt.figure(figsize=(6,6))
    # #         ax = fig.add_subplot(1,1,1, projection='polar')
        
    # #         colors = "seagreen"
    # #         bars = ax.bar(cat[:-1],freq,color=colors,width=0.5)
    # #         ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'])
    # #         ax.xaxis.grid('major')
    # #         plt.subplots_adjust(left=0.09, bottom=0.08, right=0.98, top=0.7 , hspace=0.1 )
    # #         plt.title("direction of climate velocities - Temperature \n "+str(depth)+'m \n \n occurences in % \n')         
    
    # else :
    
    # if zone == 'omz' or zone =='large':
    #     fig = plt.figure(figsize=(10,6.5))
    #     if zone == 'omz':
    #         shrink=0.6
    #     else :
            # shrink = 0.8
    # if zone == 'coast':
#     fig = plt.figure(figsize=(9,7.5))
#     shrink = 0.9
        
#     plt.suptitle('PHI Climate Velocity - RCP8.5 \n '+str(depth)+'m', fontsize=15)
#     ax = fig.add_subplot(121,projection=ccrs.PlateCarree())
#     vmin = -20; vmax = 10

#     levels = np.arange(vmin,vmax+2.5,2.5)
#     norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
#     bar_ticks = np.linspace(vmin,vmax,7)
#     m = plt.contourf(lon,lat,cv_mean,vmin=vmin,vmax=vmax,\
#                       levels=levels,extend='both', cmap='RdYlBu_r',\
#                           norm=norm)
#     crs = ccrs.PlateCarree()
#     ax.quiver(ulon[::10],ulat[::10],u_norm[::10],v_norm[::10],transform = crs,angles ="uv",\
#                   width=0.002)
#     cbar = plt.colorbar(m,ticks=bar_ticks,pad=0.10,shrink=shrink)
#     cbar.ax.set_ylabel("km/yr",fontsize=15,labelpad=5)
#     cbar.ax.tick_params(labelsize=12)
#     ax.coastlines()
#     gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                       linewidth=1, color='gray', alpha=1, linestyle='--')
#     gl.xlabels_top = False
#     gl.ylabels_right = False
#     # plt.subplots_adjust(left = 0.090,right=0.7,bottom=0.06, top = 0.9,hspace=0.9 )
#     ax.set_title('ensemble mean',pad=15)
    
#     vmin=0;vmax=10;
#     levels=np.linspace(0,vmax,11)
#     bar_ticks=np.linspace(0,vmax,6)
#     ax = fig.add_subplot(122,projection=ccrs.PlateCarree())
#     m = plt.contourf(lon,lat,cv_std,cmap='Reds',\
#                           vmin=vmin,vmax=vmax,levels=levels,extend='max')
#     cbar = plt.colorbar(m,pad=0.10,shrink=shrink)
#     cbar.ax.set_ylabel("km/yr",fontsize=15,labelpad=10)
#     cbar.ax.tick_params(labelsize=12)
#     ax.coastlines()
#     gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                       linewidth=1, color='gray', alpha=1, linestyle='--')
#     # plt.subplots_adjust(left = 0.070,right=0.99,bottom=0.06, top = 0.9,hspace=0.9 )
#     gl.xlabels_top = False
#     gl.ylabels_right = False
#     gl.ylabels_left = False
#     ax.set_title('standard deviation',pad=15)
    
#     if save=='yes':
#         plt.savefig(path_fig+str(depth)+'m_CV_'+zone+'_'+str(depth)+'.png')
    
#     """fig rosewind """
    
#     # fig1 = plt.figure(figsize=(5,5))
#     # ax = fig1.add_subplot(111, projection='polar')
#     # if zone == 'coast':
#     #     color = 'darkorange'
#     # bars = ax.bar(cat[:-1],freq,color=color,width=0.5)
#     # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'])
#     # ax.xaxis.grid('major')
#     # plt.subplots_adjust(left=0.09, bottom=0.08, right=0.98, top=0.7 , hspace=0.1 )
#     # ax.set_title(zone+'  '+str(depth)+'m \noccurences in % \n',fontsize=14,pad=2)                                      
#     # plt.subplots_adjust(left = 0.05,right = 0.9  )
#     # plt.savefig(path_fig+str(depth)+'m_ROSEWIND_'+zone+'_'+str(depth)+'.png')
                         
    
    
#     """  fig gradients """
    
# # if zone == 'coast':

#     fig2 = plt.figure(figsize=(9.5,7.5))
#     plt.suptitle(str(depth)+'m',fontsize=15)
#     shrink = 0.90

#     ### gradS
#     ax = fig2.add_subplot(121,projection=ccrs.PlateCarree())
#     vmin = 0;
#     if depth == 100:
#         vmax=10
#     else:
#         vmax=10
#     levels = np.linspace(vmin,vmax,11)
#     # norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
#     ticks = np.linspace(vmin,vmax,6)
#     m = plt.contourf(lon,lat,gradS_mean*1e3,vmin=vmin,vmax=vmax,levels=levels,\
#                       extend='max', cmap='YlGnBu')
#     crs = ccrs.PlateCarree()
#     ax.quiver(ulon[::5],ulat[::5],u_norm[::5],v_norm[::5],transform = crs,angles ="uv",\
#                   width=0.002)
#     cbar = plt.colorbar(m,pad=0.10,shrink=shrink,ticks=ticks)
#     # ticks=bar_ticks,
#     cbar.ax.set_ylabel('$ 10^{-3} $ unit' +'/km',fontsize=15,labelpad=5)
#     cbar.ax.tick_params(labelsize=12)
#     ax.coastlines()
#     gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                       linewidth=1, color='gray', alpha=1, linestyle='--')
#     gl.xlabels_top = False
#     gl.ylabels_right = False
#     # plt.subplots_adjust(left = 0.090,right=0.7,bottom=0.06, top = 0.9,hspace=0.9 )
#     ax.set_title('Spatial gradient',fontsize=14,pad=15)
    
#     ### graT
#     ax = fig2.add_subplot(122,projection=ccrs.PlateCarree())
#     vmin = -10; vmax = 10
#     levels = np.linspace(vmin,vmax,11)
#     ticks = np.linspace(vmin,vmax,11)
#     norm = DivergingNorm(vmin=vmin,vcenter=0,vmax=vmax)
#     m = plt.contourf(lon,lat,gradT_mean*1e3,vmin=vmin,vmax=vmax,levels=levels,\
#                       extend='both', cmap='RdBu_r',norm=norm)
#     crs = ccrs.PlateCarree()
#     ax.quiver(ulon[::5],ulat[::5],u_norm[::5],v_norm[::5],transform = crs,angles ="uv",\
#                   width=0.002)
#     cbar = plt.colorbar(m,pad=0.10,shrink=shrink,ticks=ticks)
#     cbar.ax.set_ylabel('$ 10^{-3} $ unit /yr',fontsize=15,labelpad=5)
#     cbar.ax.tick_params(labelsize=12)
#     ax.coastlines()
#     gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                       linewidth=1, color='gray', alpha=1, linestyle='--')
#     gl.xlabels_top = False
#     gl.ylabels_right = False
#     gl.ylabels_left = False
    
#     plt.subplots_adjust(left = 0.090,right=0.92,bottom=0.06, top = 0.9,hspace=0.9 )
#     ax.set_title('Temporal gradient',fontsize=14,pad=15)
    
#     if save=='yes':
#         print('in')
#         plt.savefig(path_fig+str(depth)+'m_GRAD_'+zone+'_'+str(depth)+'.png')