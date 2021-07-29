#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 16:46:28 2021

@author: parouffe
"""

import numpy as np

def metabolic_index_ST(o2,Temp,salt,zt,depth_lvl,a0,e0):
                       
    
    """
    z_t in meters
    a0 in atm-1
    e0 in eV
    """
    print('Phi')
    # Parameters of Metabolic Index
    # a0 = 1/a0
    A0 = a0/1.01325e3     # 1/atm to 1/mbar with A0=25 atm-1 #*1.01325e3
    # print('A0 (mbar)-1 ',A0)

    E0 = e0  *1.60218e-19                # J
    # print('E0 ',E0)

    kB = 1.38064852*1e-23 #8.617333262e-5       # J/K # eV/K
    Tref = 15 + 273.15        # Celsius to Kelvin
    Tabs = Temp + 273.15      # Kelvin
    co2 = o2                  # mmol/m3 == micromol/L
    # S = 34                    # sans unité
    Spreset = 0               # sans unité
    
    # salinity coorection
    Ts = np.log( (298.15-Temp)/(273.15+Temp))  # sans unité , np.log(x) = ln(x)
                                               # with Temp in celsius
    b0 = -6.24523e-3
    b1 = -7.37614e-3
    b2 = -1.03410e-2
    b3 = -8.17083e-3
    c0 = -4.88682e-7
    
    Scorr = np.exp( (salt-Spreset)*(b0 + b1*Ts + b2*Ts**2. + b3*Ts**3.) + c0*(salt**2. - Spreset**2.) )
    # Scorr sans dimension
    
    # pH20
    d0 = 24.4553
    d1 = -67.4509
    d2 = -4.8489
    d3 = -5.44e-4
    
    # ph20 en mbar
    pH2O = 1013.25*np.exp(d0 + d1*(100/Tabs) + d2*np.log(Tabs/100) + d3*salt) # mbar
    pH2O_preset = 1013.25*np.exp(d0 + d1*(100/Tabs) + d2*np.log(Tabs/100) + d3*Spreset) # mbar
    co2_corr = Scorr*((1013.25-pH2O_preset)/(1013.25-pH2O)) * co2
    
    # Temperature correction
    xo2 = 0.20946  # sans dimension
    a0 = 2.00907
    a1 = 3.22014
    a2 = 4.05010
    a3 = 4.94457
    a4 = -2.56848e-1
    a5 = 3.88767
    
    Tcorr = 44.6596 * np.exp(a0 + a1*Ts + a2*(Ts**2) + a3*(Ts**3.) + a4*(Ts**4.) + a5*(Ts**5.))
    # sans dim ?
    
    # calcul pO2:
    grav = 980.616   # cm/s2
    rho_sw = 1.026   # g/cm3
    
    # pression hydrostatique
    P = zt * grav/100 *rho_sw*1000  # (vector) : hydrostatic pressure at midpoint of layer
        #    m    * cm/s2   * g/cm3 = g/(s2*cm) = 0.01 Pa 
        # to  m    *  m/s2   * kg/m3 = kg/(m/s2) = Pa
    P = P*(10**-4)     # hydrostatic pressure in Pa to dbar
    Vm = 0.317         # molar volume of oxygen (m3 mol-1 Pa dbar-1)
    R =  8.314462618   # universal gaz constant (J /mol /K) 
    
    a = xo2*(1013.25-pH2O)
    b = (Tcorr)*Scorr
    
    c = (Vm*P[depth_lvl])/(R*Tabs)
    d = np.exp(c)
    
    pO2 = co2_corr * a/b * d   # mbar # co2 salinity corrected 
    
    # Metabolic Index
    supply = A0*pO2
    exponant = np.array(-E0/kB*(1/Tabs-1/Tref))
    # exponant2 = np.array(exponant,dtype=np.float64)
    demand = np.exp(exponant)
    phi = supply / demand
    
    return phi