#Compare painted data with observed data - at different radii, for different distributions of GS gas
# Works but could do with some tidying up 
import numpy as np
import h5py
import pandas as pd
import math
from astropy.io import fits
from astropy.table import Table, join
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import matplotlib.colors as colors
import sys
sys.path.append('./scripts/')
from chemevo import *

fl = chem_evo_data('./nomerger.hdf5')
gl = chem_evo_data('./centremerger.hdf5')
hl = chem_evo_data('./solargaussian.hdf5')
il = chem_evo_data('./doublemerger.hdf5')



data_file_1 = '/data/ktfm2/apogee_data/apogee_astroNN_DR16.fits'


hdu_list_1 = fits.open(data_file_1, memmap=True)
apogee_data = Table(hdu_list_1[1].data)
hdu_list_1.close()

def betw(x,l,u):
    return (x>l)&(x<u)

def outs(x,l,u):
    return (x<l)|(x>u)

#===================================================================================================================

#Radial Migration filter

fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['age_lowess_correct']))&(apogee_data['age_lowess_correct']>0.0)&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(~pd.isna(apogee_data['FE_H_ERR']))&(apogee_data['LOGG']<3.5)&(outs(apogee_data['GALZ'],-1.0,1.0))&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],7.6,8.6))

lower_fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['age_lowess_correct']))&(apogee_data['age_lowess_correct']>0.0)&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(~pd.isna(apogee_data['FE_H_ERR']))&(apogee_data['LOGG']<3.5)&(outs(apogee_data['GALZ'],-1.0,1.0))&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],4.6,5.6))

upper_fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['age_lowess_correct']))&(apogee_data['age_lowess_correct']>0.0)&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(~pd.isna(apogee_data['FE_H_ERR']))&(apogee_data['LOGG']<3.5)&(outs(apogee_data['GALZ'],-1.0,1.0))&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],10.6,11.6))


#=====================================================================================================================

def stretchy(x):
    y = [13.7*(np.tanh((a-4.0)/3.0)) if a > 5.198 else a for a in x]
    return np.array(y)



#====================================================================================================================

#Radial migration, randomness from gaussian
rl_random = np.random.normal(0.0, 4.0*((apogee_data['age_lowess_correct'][fltr]/13.7)**0.5))
rl_scatter = apogee_data['rl'][fltr]+rl_random

lower_rl_random = np.random.normal(0.0, 4.0*((apogee_data['age_lowess_correct'][lower_fltr]/13.7)**0.5))
lower_rl_scatter = apogee_data['rl'][lower_fltr]+lower_rl_random

upper_rl_random = np.random.normal(0.0, 4.0*((apogee_data['age_lowess_correct'][upper_fltr]/13.7)**0.5))
upper_rl_scatter = apogee_data['rl'][upper_fltr]+upper_rl_random



age_stretch_rl_random = np.random.normal(0.0, 4.0*((stretchy(apogee_data['age_lowess_correct'][fltr])/13.7)**0.5))
age_stretch_rl_scatter = apogee_data['rl'][fltr]+age_stretch_rl_random

lower_age_stretch_rl_random = np.random.normal(0.0, 4.0*((stretchy(apogee_data['age_lowess_correct'][lower_fltr])/13.7)**0.5))
lower_age_stretch_rl_scatter = apogee_data['rl'][lower_fltr]+lower_age_stretch_rl_random

upper_age_stretch_rl_random = np.random.normal(0.0, 4.0*((stretchy(apogee_data['age_lowess_correct'][upper_fltr])/13.7)**0.5))
upper_age_stretch_rl_scatter = apogee_data['rl'][upper_fltr]+upper_age_stretch_rl_random


#====================================================================================================================


#Radial migration painting

rl_scatter_model_fe = fl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Fe'])
rl_scatter_model_mg = fl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Mg']) 

lower_rl_scatter_model_fe = fl.paint(lower_rl_scatter,(13.7-apogee_data['age_lowess_correct'][lower_fltr]),['Fe'])
lower_rl_scatter_model_mg = fl.paint(lower_rl_scatter,(13.7-apogee_data['age_lowess_correct'][lower_fltr]),['Mg']) 

upper_rl_scatter_model_fe = fl.paint(upper_rl_scatter,(13.7-apogee_data['age_lowess_correct'][upper_fltr]),['Fe'])
upper_rl_scatter_model_mg = fl.paint(upper_rl_scatter,(13.7-apogee_data['age_lowess_correct'][upper_fltr]),['Mg']) 



age_stretch_model_fe = fl.paint(age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][fltr])),['Fe'])
age_stretch_model_mg = fl.paint(age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][fltr])),['Mg']) 

lower_age_stretch_model_fe = fl.paint(lower_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][lower_fltr])),['Fe'])
lower_age_stretch_model_mg = fl.paint(lower_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][lower_fltr])),['Mg']) 

upper_age_stretch_model_fe = fl.paint(upper_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][upper_fltr])),['Fe'])
upper_age_stretch_model_mg = fl.paint(upper_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][upper_fltr])),['Mg']) 


#===========================================================================================================================

rl_scatter_model_fe1 = gl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Fe'])
rl_scatter_model_mg1 = gl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Mg']) 

lower_rl_scatter_model_fe1 = gl.paint(lower_rl_scatter,(13.7-apogee_data['age_lowess_correct'][lower_fltr]),['Fe'])
lower_rl_scatter_model_mg1 = gl.paint(lower_rl_scatter,(13.7-apogee_data['age_lowess_correct'][lower_fltr]),['Mg']) 

upper_rl_scatter_model_fe1 = gl.paint(upper_rl_scatter,(13.7-apogee_data['age_lowess_correct'][upper_fltr]),['Fe'])
upper_rl_scatter_model_mg1 = gl.paint(upper_rl_scatter,(13.7-apogee_data['age_lowess_correct'][upper_fltr]),['Mg']) 



age_stretch_model_fe1 = gl.paint(age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][fltr])),['Fe'])
age_stretch_model_mg1 = gl.paint(age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][fltr])),['Mg']) 

lower_age_stretch_model_fe1 = gl.paint(lower_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][lower_fltr])),['Fe'])
lower_age_stretch_model_mg1 = gl.paint(lower_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][lower_fltr])),['Mg']) 

upper_age_stretch_model_fe1 = gl.paint(upper_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][upper_fltr])),['Fe'])
upper_age_stretch_model_mg1 = gl.paint(upper_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][upper_fltr])),['Mg']) 

#=====================================================================================================================


rl_scatter_model_fe2 = hl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Fe'])
rl_scatter_model_mg2 = hl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Mg']) 

lower_rl_scatter_model_fe2 = hl.paint(lower_rl_scatter,(13.7-apogee_data['age_lowess_correct'][lower_fltr]),['Fe'])
lower_rl_scatter_model_mg2 = hl.paint(lower_rl_scatter,(13.7-apogee_data['age_lowess_correct'][lower_fltr]),['Mg']) 

upper_rl_scatter_model_fe2 = hl.paint(upper_rl_scatter,(13.7-apogee_data['age_lowess_correct'][upper_fltr]),['Fe'])
upper_rl_scatter_model_mg2 = hl.paint(upper_rl_scatter,(13.7-apogee_data['age_lowess_correct'][upper_fltr]),['Mg']) 



age_stretch_model_fe2 = hl.paint(age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][fltr])),['Fe'])
age_stretch_model_mg2 = hl.paint(age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][fltr])),['Mg']) 

lower_age_stretch_model_fe2 = hl.paint(lower_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][lower_fltr])),['Fe'])
lower_age_stretch_model_mg2 = hl.paint(lower_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][lower_fltr])),['Mg']) 

upper_age_stretch_model_fe2 = hl.paint(upper_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][upper_fltr])),['Fe'])
upper_age_stretch_model_mg2 = hl.paint(upper_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][upper_fltr])),['Mg']) 

#===================================================================================================================

rl_scatter_model_fe3 = il.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Fe']) 
rl_scatter_model_mg3 = il.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Mg']) 

lower_rl_scatter_model_fe3 = il.paint(lower_rl_scatter,(13.7-apogee_data['age_lowess_correct'][lower_fltr]),['Fe'])
lower_rl_scatter_model_mg3 = il.paint(lower_rl_scatter,(13.7-apogee_data['age_lowess_correct'][lower_fltr]),['Mg']) 

upper_rl_scatter_model_fe3 = il.paint(upper_rl_scatter,(13.7-apogee_data['age_lowess_correct'][upper_fltr]),['Fe'])
upper_rl_scatter_model_mg3 = il.paint(upper_rl_scatter,(13.7-apogee_data['age_lowess_correct'][upper_fltr]),['Mg']) 



age_stretch_model_fe3 = il.paint(age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][fltr])),['Fe'])
age_stretch_model_mg3 = il.paint(age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][fltr])),['Mg']) 

lower_age_stretch_model_fe3 = il.paint(lower_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][lower_fltr])),['Fe'])
lower_age_stretch_model_mg3 = il.paint(lower_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][lower_fltr])),['Mg']) 

upper_age_stretch_model_fe3 = il.paint(upper_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][upper_fltr])),['Fe'])
upper_age_stretch_model_mg3 = il.paint(upper_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][upper_fltr])),['Mg']) 


#==================================================================================================================



#model_fe_uncert = np.random.normal(0.0, apogee_data['FE_H_ERR'][fltr])
#model_mg_uncert = np.random.normal(0.0, apogee_data['MG_H_ERR'][fltr])
for_fe_scatter = np.ones(len(apogee_data['FE_H_ERR'][fltr]))
for_mg_scatter = np.ones(len(apogee_data['MG_H_ERR'][fltr]))

lower_for_fe_scatter = np.ones(len(apogee_data['FE_H_ERR'][lower_fltr]))
lower_for_mg_scatter = np.ones(len(apogee_data['MG_H_ERR'][lower_fltr]))

upper_for_fe_scatter = np.ones(len(apogee_data['FE_H_ERR'][upper_fltr]))
upper_for_mg_scatter = np.ones(len(apogee_data['MG_H_ERR'][upper_fltr]))


model_fe_uncert = np.random.normal(0.0, (0.05*for_fe_scatter))
model_mg_uncert = np.random.normal(0.0, (0.005*for_mg_scatter))

lower_model_fe_uncert = np.random.normal(0.0, (0.05*lower_for_fe_scatter))
lower_model_mg_uncert = np.random.normal(0.0, (0.005*lower_for_mg_scatter))

upper_model_fe_uncert = np.random.normal(0.0, (0.05*upper_for_fe_scatter))
upper_model_mg_uncert = np.random.normal(0.0, (0.005*upper_for_mg_scatter))


#Changed for radial effects data set
model_fe_random = rl_scatter_model_fe['Fe_H'] + model_fe_uncert
model_mg_random = rl_scatter_model_mg['Mg_H'] + model_mg_uncert

lower_model_fe_random = lower_rl_scatter_model_fe['Fe_H'] + lower_model_fe_uncert
lower_model_mg_random = lower_rl_scatter_model_mg['Mg_H'] + lower_model_mg_uncert

upper_model_fe_random = upper_rl_scatter_model_fe['Fe_H'] + upper_model_fe_uncert
upper_model_mg_random = upper_rl_scatter_model_mg['Mg_H'] + upper_model_mg_uncert


age_stretch_model_fe_random = age_stretch_model_fe['Fe_H'] + model_fe_uncert
age_stretch_model_mg_random = age_stretch_model_mg['Mg_H'] + model_mg_uncert

lower_age_stretch_model_fe_random = lower_age_stretch_model_fe['Fe_H'] + lower_model_fe_uncert
lower_age_stretch_model_mg_random = lower_age_stretch_model_mg['Mg_H'] + lower_model_mg_uncert

upper_age_stretch_model_fe_random = upper_age_stretch_model_fe['Fe_H'] + upper_model_fe_uncert
upper_age_stretch_model_mg_random = upper_age_stretch_model_mg['Mg_H'] + upper_model_mg_uncert


#=======================================================================================================
model_fe_random1 = rl_scatter_model_fe1['Fe_H'] + model_fe_uncert
model_mg_random1 = rl_scatter_model_mg1['Mg_H'] + model_mg_uncert

lower_model_fe_random1 = lower_rl_scatter_model_fe1['Fe_H'] + lower_model_fe_uncert
lower_model_mg_random1 = lower_rl_scatter_model_mg1['Mg_H'] + lower_model_mg_uncert

upper_model_fe_random1 = upper_rl_scatter_model_fe1['Fe_H'] + upper_model_fe_uncert
upper_model_mg_random1 = upper_rl_scatter_model_mg1['Mg_H'] + upper_model_mg_uncert


age_stretch_model_fe_random1 = age_stretch_model_fe1['Fe_H'] + model_fe_uncert
age_stretch_model_mg_random1 = age_stretch_model_mg1['Mg_H'] + model_mg_uncert

lower_age_stretch_model_fe_random1 = lower_age_stretch_model_fe1['Fe_H'] + lower_model_fe_uncert
lower_age_stretch_model_mg_random1 = lower_age_stretch_model_mg1['Mg_H'] + lower_model_mg_uncert

upper_age_stretch_model_fe_random1 = upper_age_stretch_model_fe1['Fe_H'] + upper_model_fe_uncert
upper_age_stretch_model_mg_random1 = upper_age_stretch_model_mg1['Mg_H'] + upper_model_mg_uncert

#======================================================================================================
model_fe_random2 = rl_scatter_model_fe2['Fe_H'] + model_fe_uncert
model_mg_random2 = rl_scatter_model_mg2['Mg_H'] + model_mg_uncert

lower_model_fe_random2 = lower_rl_scatter_model_fe2['Fe_H'] + lower_model_fe_uncert
lower_model_mg_random2 = lower_rl_scatter_model_mg2['Mg_H'] + lower_model_mg_uncert

upper_model_fe_random2 = upper_rl_scatter_model_fe2['Fe_H'] + upper_model_fe_uncert
upper_model_mg_random2 = upper_rl_scatter_model_mg2['Mg_H'] + upper_model_mg_uncert


age_stretch_model_fe_random2 = age_stretch_model_fe2['Fe_H'] + model_fe_uncert
age_stretch_model_mg_random2 = age_stretch_model_mg2['Mg_H'] + model_mg_uncert

lower_age_stretch_model_fe_random2 = lower_age_stretch_model_fe2['Fe_H'] + lower_model_fe_uncert
lower_age_stretch_model_mg_random2 = lower_age_stretch_model_mg2['Mg_H'] + lower_model_mg_uncert

upper_age_stretch_model_fe_random2 = upper_age_stretch_model_fe2['Fe_H'] + upper_model_fe_uncert
upper_age_stretch_model_mg_random2 = upper_age_stretch_model_mg2['Mg_H'] + upper_model_mg_uncert
#====================================================================================================================
model_fe_random3 = rl_scatter_model_fe3['Fe_H'] + model_fe_uncert
model_mg_random3 = rl_scatter_model_mg3['Mg_H'] + model_mg_uncert

lower_model_fe_random3 = lower_rl_scatter_model_fe3['Fe_H'] + lower_model_fe_uncert
lower_model_mg_random3 = lower_rl_scatter_model_mg3['Mg_H'] + lower_model_mg_uncert

upper_model_fe_random3 = upper_rl_scatter_model_fe3['Fe_H'] + upper_model_fe_uncert
upper_model_mg_random3 = upper_rl_scatter_model_mg3['Mg_H'] + upper_model_mg_uncert


age_stretch_model_fe_random3 = age_stretch_model_fe3['Fe_H'] + model_fe_uncert
age_stretch_model_mg_random3 = age_stretch_model_mg3['Mg_H'] + model_mg_uncert

lower_age_stretch_model_fe_random3 = lower_age_stretch_model_fe3['Fe_H'] + lower_model_fe_uncert
lower_age_stretch_model_mg_random3 = lower_age_stretch_model_mg3['Mg_H'] + lower_model_mg_uncert

upper_age_stretch_model_fe_random3 = upper_age_stretch_model_fe3['Fe_H'] + upper_model_fe_uncert
upper_age_stretch_model_mg_random3 = upper_age_stretch_model_mg3['Mg_H'] + upper_model_mg_uncert

#=======================================================================================================================
#Density plots with Radial Migration and scatter in both
f,a=plt.subplots(1,5,figsize=[15.,4.],sharex=True,sharey=True)
plt.sca(a[0])
plt.hist2d(apogee_data['FE_H'][fltr],apogee_data['MG_H'][fltr]-apogee_data['FE_H'][fltr],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.ylabel('[Mg/Fe]',fontsize=12)
plt.title('Density plot of data',fontsize=12)
plt.sca(a[1])
plt.hist2d(age_stretch_model_fe_random,age_stretch_model_mg_random-age_stretch_model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('No merger',fontsize=12)
plt.sca(a[2])
plt.hist2d(age_stretch_model_fe_random1,age_stretch_model_mg_random1-age_stretch_model_fe_random1,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Central Exp',fontsize=12)
plt.sca(a[3])
plt.hist2d(age_stretch_model_fe_random2,age_stretch_model_mg_random2-age_stretch_model_fe_random2,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Solar Gaussian',fontsize=12)
plt.sca(a[4])
plt.hist2d(age_stretch_model_fe_random3,age_stretch_model_mg_random3-age_stretch_model_fe_random3,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Double Gaussian',fontsize=12)
#f.suptitle('Comparison between observed and model predictions, w/ RM, w/Both Scatter')
#plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/RMObsModBothDensity.pdf', bbox_inches='tight')
plt.show()

#==============================================================================================================
#lower radii
#Density plots with Radial Migration and scatter in both
f,a=plt.subplots(1,5,figsize=[15.,4.],sharex=True,sharey=True)
plt.sca(a[0])
plt.hist2d(apogee_data['FE_H'][lower_fltr],apogee_data['MG_H'][lower_fltr]-apogee_data['FE_H'][lower_fltr],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.ylabel('[Mg/Fe]',fontsize=12)
plt.title('Density plot of data L',fontsize=12)
plt.sca(a[1])
plt.hist2d(lower_age_stretch_model_fe_random,lower_age_stretch_model_mg_random-lower_age_stretch_model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('No Merger L',fontsize=12)
plt.sca(a[2])
plt.hist2d(lower_age_stretch_model_fe_random1,lower_age_stretch_model_mg_random1-lower_age_stretch_model_fe_random1,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Central Exp L',fontsize=12)
plt.sca(a[3])
plt.hist2d(lower_age_stretch_model_fe_random2,lower_age_stretch_model_mg_random2-lower_age_stretch_model_fe_random2,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Solar Gaussian L',fontsize=12)
plt.sca(a[4])
plt.hist2d(lower_age_stretch_model_fe_random3,lower_age_stretch_model_mg_random3-lower_age_stretch_model_fe_random3,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Double Gaussian L',fontsize=12)
#f.suptitle('Comparison between observed and model predictions, w/ RM, w/Both Scatter')
#plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/RMObsModBothDensity.pdf', bbox_inches='tight')
plt.show()

#==============================================================================================================
#upper radii
#Density plots with Radial Migration and scatter in both
f,a=plt.subplots(1,5,figsize=[15.,4.],sharex=True,sharey=True)
plt.sca(a[0])
plt.hist2d(apogee_data['FE_H'][upper_fltr],apogee_data['MG_H'][upper_fltr]-apogee_data['FE_H'][upper_fltr],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.ylabel('[Mg/Fe]',fontsize=12)
plt.title('Density plot of data',fontsize=12)
plt.sca(a[1])
plt.hist2d(upper_age_stretch_model_fe_random,upper_age_stretch_model_mg_random-upper_age_stretch_model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('No Merger U',fontsize=12)
plt.sca(a[2])
plt.hist2d(upper_age_stretch_model_fe_random1,upper_age_stretch_model_mg_random1-upper_age_stretch_model_fe_random1,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Central Exp U',fontsize=12)
plt.sca(a[3])
plt.hist2d(upper_age_stretch_model_fe_random2,upper_age_stretch_model_mg_random2-upper_age_stretch_model_fe_random2,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Solar Gaussian U',fontsize=12)
plt.sca(a[4])
plt.hist2d(upper_age_stretch_model_fe_random3,upper_age_stretch_model_mg_random3-upper_age_stretch_model_fe_random3,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Double Gaussian U',fontsize=12)
#f.suptitle('Comparison between observed and model predictions, w/ RM, w/Both Scatter')
#plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/RMObsModBothDensity.pdf', bbox_inches='tight')
plt.show()

#==============================================================================================================



