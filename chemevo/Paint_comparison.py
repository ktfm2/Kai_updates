#Compare painted data with observed data - for three different sets of ages
#Works but could do with some tidying up of the code
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

#fl = chem_evo_data('./comparison.hdf5')
#fl = chem_evo_data('./output.hdf5')
fl = chem_evo_data('./KSsfr.hdf5')


hdf5_file = '/data/ktfm2/apogee_data/gaia_spectro.hdf5'
data_file_1 = '/data/ktfm2/apogee_data/apogee_astroNN_DR16.fits'
data_file_2 = '/data/jls/apokasc_astroNN.fits'

hdf = h5py.File(hdf5_file, "r")
dataset = hdf['data']

log_age_data = dataset["log10_age"]
ID_data = dataset["APOGEE_ID"]

SD_table = Table([ID_data, log_age_data], names=('apogee_id','log_age_data'))


hdu_list_1 = fits.open(data_file_1, memmap=True)
apogee_data = Table(hdu_list_1[1].data)
hdu_list_1.close()

hdu_list_2 = fits.open(data_file_2, memmap=True)
apokasc_data = Table(hdu_list_2[1].data)
hdu_list_2.close()

#print(apokasc_data.colnames)

def betw(x,l,u):
    return (x>l)&(x<u)

def outs(x,l,u):
    return (x<l)|(x>u)

#Join the APOGEE data table and table from Sanders & Das - will have less rows
#Comment out if not using S&D ages
#full_table = join(apogee_data, SD_table)
#apogee_data = full_table

#Use APOKASC data, comment out if not using
#apogee_data = apokasc_data


#===================================================================================================================

#Radial Migration filter
#fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['LogAge']))&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(~pd.isna(apogee_data['FE_H_ERR']))&(~pd.isna(apogee_data['MG_H_ERR']))&(apogee_data['LOGG']<3.5)&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],7.6,8.6))
#fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['log_age_data']))&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(~pd.isna(apogee_data['FE_H_ERR']))&(apogee_data['LOGG']<3.5)&(outs(apogee_data['GALZ'],-1.0,1.0))&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],7.6,8.6))


fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['age_lowess_correct']))&(apogee_data['age_lowess_correct']>0.0)&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(~pd.isna(apogee_data['FE_H_ERR']))&(apogee_data['LOGG']<3.5)&(outs(apogee_data['GALZ'],-1.0,1.0))&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],7.6,8.6))

lower_fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['age_lowess_correct']))&(apogee_data['age_lowess_correct']>0.0)&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(~pd.isna(apogee_data['FE_H_ERR']))&(apogee_data['LOGG']<3.5)&(outs(apogee_data['GALZ'],-1.0,1.0))&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],4.6,5.6))

upper_fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['age_lowess_correct']))&(apogee_data['age_lowess_correct']>0.0)&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(~pd.isna(apogee_data['FE_H_ERR']))&(apogee_data['LOGG']<3.5)&(outs(apogee_data['GALZ'],-1.0,1.0))&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],10.6,11.6))


#fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['age_linear_correct']))&(apogee_data['age_linear_correct']>0.0)&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(~pd.isna(apogee_data['FE_H_ERR']))&(apogee_data['LOGG']<3.5)&(outs(apogee_data['GALZ'],-1.0,1.0))&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],7.6,8.6))

#print(len(apogee_data['LogAge']))
#print(len(apogee_data['LogAge'][fltr]))
print(len(apogee_data['FE_H'][fltr]))


#=====================================================================================================================

def stretchy(x):
#    y = [1.4*a if a > 5.0 else a for a in x]
#    y = [13.7*(np.tanh(a/8.0)) if a < 12.5 else a for a in x]
    y = [13.7*(np.tanh((a-4.0)/3.0)) if a > 5.198 else a for a in x]
    return np.array(y)



#====================================================================================================================

#Radial migration, randomness from gaussian
#rl_random = np.random.normal(0.0, 4.0*(((10.0**(apogee_data['LogAge'][fltr]-3.0))/13.7)**0.5))
#rl_random = np.random.normal(0.0, 4.0*(((10.0**apogee_data['log_age_data'][fltr])/13.7)**0.5))
#rl_random = np.random.normal(0.0, 4.0*((apogee_data['age_linear_correct'][fltr]/13.7)**0.5))
rl_random = np.random.normal(0.0, 4.0*((apogee_data['age_lowess_correct'][fltr]/13.7)**0.5))
rl_scatter = apogee_data['rl'][fltr]+rl_random

lower_rl_random = np.random.normal(0.0, 4.0*((apogee_data['age_lowess_correct'][lower_fltr]/13.7)**0.5))
lower_rl_scatter = apogee_data['rl'][lower_fltr]+lower_rl_random

upper_rl_random = np.random.normal(0.0, 4.0*((apogee_data['age_lowess_correct'][upper_fltr]/13.7)**0.5))
upper_rl_scatter = apogee_data['rl'][upper_fltr]+upper_rl_random



#age_stretch_rl_random = np.random.normal(0.0, 4.0*(((1.3*apogee_data['age_lowess_correct'][fltr])/13.7)**0.5))
#age_stretch_rl_random = np.random.normal(0.0, 4.0*((stretchy(apogee_data['age_lowess_correct'][fltr])/13.7)**0.5))
age_stretch_rl_random = np.random.normal(0.0, 4.0*((stretchy(apogee_data['age_lowess_correct'][fltr])/13.7)**0.5))
age_stretch_rl_scatter = apogee_data['rl'][fltr]+age_stretch_rl_random

lower_age_stretch_rl_random = np.random.normal(0.0, 4.0*((stretchy(apogee_data['age_lowess_correct'][lower_fltr])/13.7)**0.5))
lower_age_stretch_rl_scatter = apogee_data['rl'][lower_fltr]+lower_age_stretch_rl_random

upper_age_stretch_rl_random = np.random.normal(0.0, 4.0*((stretchy(apogee_data['age_lowess_correct'][upper_fltr])/13.7)**0.5))
upper_age_stretch_rl_scatter = apogee_data['rl'][upper_fltr]+upper_age_stretch_rl_random


#====================================================================================================================

#Radial migration painting
#rl_scatter_model_fe = fl.paint(rl_scatter,(13.7-(10.0**(apogee_data['LogAge'][fltr]-3.0))),['Fe'])
#rl_scatter_model_mg = fl.paint(rl_scatter,(13.7-(10.0**(apogee_data['LogAge'][fltr]-3.0))),['Mg']) 

#rl_scatter_model_fe = fl.paint(rl_scatter,(13.7-(10.0**apogee_data['log_age_data'][fltr])),['Fe'])
#rl_scatter_model_mg = fl.paint(rl_scatter,(13.7-(10.0**apogee_data['log_age_data'][fltr])),['Mg']) 

rl_scatter_model_fe = fl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Fe'])
rl_scatter_model_mg = fl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Mg']) 

lower_rl_scatter_model_fe = fl.paint(lower_rl_scatter,(13.7-apogee_data['age_lowess_correct'][lower_fltr]),['Fe'])
lower_rl_scatter_model_mg = fl.paint(lower_rl_scatter,(13.7-apogee_data['age_lowess_correct'][lower_fltr]),['Mg']) 

upper_rl_scatter_model_fe = fl.paint(upper_rl_scatter,(13.7-apogee_data['age_lowess_correct'][upper_fltr]),['Fe'])
upper_rl_scatter_model_mg = fl.paint(upper_rl_scatter,(13.7-apogee_data['age_lowess_correct'][upper_fltr]),['Mg']) 


#rl_scatter_model_fe = fl.paint(rl_scatter,(13.7-apogee_data['age_linear_correct'][fltr]),['Fe'])
#rl_scatter_model_mg = fl.paint(rl_scatter,(13.7-apogee_data['age_linear_correct'][fltr]),['Mg']) 

#age_stretch_model_fe = fl.paint(age_stretch_rl_scatter,(13.7-(1.3*apogee_data['age_lowess_correct'][fltr])),['Fe'])
#age_stretch_model_mg = fl.paint(age_stretch_rl_scatter,(13.7-(1.3*apogee_data['age_lowess_correct'][fltr])),['Mg']) 

age_stretch_model_fe = fl.paint(age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][fltr])),['Fe'])
age_stretch_model_mg = fl.paint(age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][fltr])),['Mg']) 

lower_age_stretch_model_fe = fl.paint(lower_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][lower_fltr])),['Fe'])
lower_age_stretch_model_mg = fl.paint(lower_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][lower_fltr])),['Mg']) 

upper_age_stretch_model_fe = fl.paint(upper_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][upper_fltr])),['Fe'])
upper_age_stretch_model_mg = fl.paint(upper_age_stretch_rl_scatter,(13.7-stretchy(apogee_data['age_lowess_correct'][upper_fltr])),['Mg']) 


#=======================================================================================================================

#Plot radial migration - no scattering abundances
plt.scatter(apogee_data['FE_H'][fltr],apogee_data['MG_H'][fltr]-apogee_data['FE_H'][fltr],s=2.0,color='b')
plt.scatter(rl_scatter_model_fe['Fe_H'],rl_scatter_model_mg['Mg_H']-rl_scatter_model_fe['Fe_H'],s=2.0,color='r')
plt.title('Radial Migration Effects',fontsize=12)
plt.xlabel('[Fe/H]',fontsize=12)
plt.ylabel('[Mg/Fe]',fontsize=12)
plt.xlim(-1.5,0.5)
plt.ylim(-0.2,0.5)
plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/RMObsModNo.pdf', bbox_inches='tight')
plt.show()

#=======================================================================================================================

#Changed for with radial effects set
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


#==============================================================================================================

#Scatter both 
plt.scatter(apogee_data['FE_H'][fltr],apogee_data['MG_H'][fltr]-apogee_data['FE_H'][fltr],s=2.0, color='b')
plt.scatter(model_fe_random,model_mg_random-model_fe_random,s=2.0,color='r')
plt.title('Scatter both Fe and Mg',fontsize=12)
plt.xlabel('[Fe/H]',fontsize=12)
plt.ylabel('[Mg/Fe]',fontsize=12)
plt.xlim(-1.5,0.5)
plt.ylim(-0.2,0.5)
plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/RMObsModBoth.pdf', bbox_inches='tight')
plt.show()

#=============================================================================================================

#Scatter both and age stretch 
plt.scatter(apogee_data['FE_H'][fltr],apogee_data['MG_H'][fltr]-apogee_data['FE_H'][fltr],s=2.0, color='b')
plt.scatter(age_stretch_model_fe_random,age_stretch_model_mg_random-age_stretch_model_fe_random,s=2.0,color='r')
plt.title('Age stretch and scatter both Fe and Mg',fontsize=12)
plt.xlabel('[Fe/H]',fontsize=12)
plt.ylabel('[Mg/Fe]',fontsize=12)
plt.xlim(-1.5,0.5)
plt.ylim(-0.2,0.5)
#plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/RMObsModBoth.pdf', bbox_inches='tight')
plt.show()

#=============================================================================================================


#Scatter Only Fe
#plt.scatter(apogee_data['FE_H'][fltr],apogee_data['MG_H'][fltr]-apogee_data['FE_H'][fltr],s=2.0, color='b')
#plt.scatter(model_fe_random,rl_scatter_model_mg['Mg_H']-model_fe_random,s=2.0,color='r')
#plt.title('Scatter only Fe, but in both axes')
#plt.xlabel('[Fe/H]')
#plt.ylabel('[Mg/Fe]')
#plt.show()

#==============================================================================================================

#Scattering Mg only
#plt.scatter(apogee_data['FE_H'][fltr],apogee_data['MG_H'][fltr]-apogee_data['FE_H'][fltr],s=2.0, color='b')
#plt.scatter(rl_scatter_model_fe['Fe_H'],model_mg_random-rl_scatter_model_fe['Fe_H'],s=2.0,color='r')
#plt.title('Scatter only Mg, only affects y-axis')
#plt.xlabel('[Fe/H]')
#plt.ylabel('[Mg/Fe]')
#plt.show()

#==============================================================================================================

#Density plots with Radial Migration but no scatter
#f,a=plt.subplots(1,2,figsize=[15.,3.],sharex=True,sharey=True)
#plt.sca(a[0])
#plt.hist2d(apogee_data['FE_H'][fltr],apogee_data['MG_H'][fltr]-apogee_data['FE_H'][fltr],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.xlabel('[Fe/H]')
#plt.ylabel('[Mg/H]')
#plt.title('Density plot of data')
#plt.sca(a[1])
#plt.hist2d(rl_scatter_model_fe['Fe_H'],rl_scatter_model_mg['Mg_H']-rl_scatter_model_fe['Fe_H'],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.xlabel('[Fe/H]')
#plt.title('Density plot of painted data')
#f.suptitle('Comparison between observed and model predictions, w/ RM, wo/Scatter')
#plt.show()

#===============================================================================================================

#Density plots with Radial Migration and scatter in Fe
#f,a=plt.subplots(1,2,figsize=[15.,3.],sharex=True,sharey=True)
#plt.sca(a[0])
#plt.hist2d(apogee_data['FE_H'][fltr],apogee_data['MG_H'][fltr]-apogee_data['FE_H'][fltr],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.xlabel('[Fe/H]',fontsize=12)
#plt.ylabel('[Mg/H]',fontsize=12)
#plt.title('Density plot of data')
#plt.sca(a[1])
#plt.hist2d(model_fe_random,rl_scatter_model_mg['Mg_H']-model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.xlabel('[Fe/H]',fontsize=12)
#plt.title('Density plot of painted data')
#plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/RMObsModBothDensity.pdf', bbox_inches='tight')
#f.suptitle('Comparison between observed and model predictions, w/ RM, w/Fe Scatter')
#plt.show()

#===============================================================================================================

#Density plots with Radial Migration and scatter in both
f,a=plt.subplots(1,3,figsize=[15.,4.],sharex=True,sharey=True)
plt.sca(a[0])
plt.hist2d(apogee_data['FE_H'][fltr],apogee_data['MG_H'][fltr]-apogee_data['FE_H'][fltr],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.ylabel('[Mg/Fe]',fontsize=12)
plt.title('Density plot of data',fontsize=12)
plt.sca(a[1])
plt.hist2d(model_fe_random,model_mg_random-model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.hist2d(rl_scatter_model_fe['Fe_H'],rl_scatter_model_mg['Mg_H']-rl_scatter_model_fe['Fe_H'],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Density plot of painted data',fontsize=12)
plt.sca(a[2])
plt.hist2d(age_stretch_model_fe_random,age_stretch_model_mg_random-age_stretch_model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.hist2d(age_stretch_model_fe['Fe_H'],age_stretch_model_mg['Mg_H']-age_stretch_model_fe['Fe_H'],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Density plot of painted data with stretched ages',fontsize=12)
#f.suptitle('Comparison between observed and model predictions, w/ RM, w/Both Scatter')
plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/RMObsModBothDensity.pdf', bbox_inches='tight')
plt.show()

#==============================================================================================================
#lower radii
#Density plots with Radial Migration and scatter in both
f,a=plt.subplots(1,3,figsize=[15.,4.],sharex=True,sharey=True)
plt.sca(a[0])
plt.hist2d(apogee_data['FE_H'][lower_fltr],apogee_data['MG_H'][lower_fltr]-apogee_data['FE_H'][lower_fltr],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.ylabel('[Mg/Fe]',fontsize=12)
plt.title('Density plot of data',fontsize=12)
plt.sca(a[1])
plt.hist2d(lower_model_fe_random,lower_model_mg_random-lower_model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.hist2d(rl_scatter_model_fe['Fe_H'],rl_scatter_model_mg['Mg_H']-rl_scatter_model_fe['Fe_H'],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Density plot of painted data',fontsize=12)
plt.sca(a[2])
plt.hist2d(lower_age_stretch_model_fe_random,lower_age_stretch_model_mg_random-lower_age_stretch_model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.hist2d(age_stretch_model_fe['Fe_H'],age_stretch_model_mg['Mg_H']-age_stretch_model_fe['Fe_H'],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Density plot of painted data with stretched ages',fontsize=12)
#f.suptitle('Comparison between observed and model predictions, w/ RM, w/Both Scatter')
#plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/RMObsModBothDensity.pdf', bbox_inches='tight')
plt.show()

#==============================================================================================================
#upper radii
#Density plots with Radial Migration and scatter in both
f,a=plt.subplots(1,3,figsize=[15.,4.],sharex=True,sharey=True)
plt.sca(a[0])
plt.hist2d(apogee_data['FE_H'][upper_fltr],apogee_data['MG_H'][upper_fltr]-apogee_data['FE_H'][upper_fltr],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.ylabel('[Mg/Fe]',fontsize=12)
plt.title('Density plot of data',fontsize=12)
plt.sca(a[1])
plt.hist2d(upper_model_fe_random,upper_model_mg_random-upper_model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.hist2d(rl_scatter_model_fe['Fe_H'],rl_scatter_model_mg['Mg_H']-rl_scatter_model_fe['Fe_H'],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Density plot of painted data',fontsize=12)
plt.sca(a[2])
plt.hist2d(upper_age_stretch_model_fe_random,upper_age_stretch_model_mg_random-upper_age_stretch_model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.hist2d(age_stretch_model_fe['Fe_H'],age_stretch_model_mg['Mg_H']-age_stretch_model_fe['Fe_H'],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]',fontsize=12)
plt.title('Density plot of painted data with stretched ages',fontsize=12)
#f.suptitle('Comparison between observed and model predictions, w/ RM, w/Both Scatter')
#plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/RMObsModBothDensity.pdf', bbox_inches='tight')
plt.show()

#==============================================================================================================


#Density plots with Radial Migration and scatter in Mg
#f,a=plt.subplots(1,2,figsize=[15.,3.],sharex=True,sharey=True)
#plt.sca(a[0])
#plt.hist2d(apogee_data['FE_H'][fltr],apogee_data['MG_H'][fltr]-apogee_data['FE_H'][fltr],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.xlabel('[Fe/H]')
#plt.ylabel('[Mg/H]')
#plt.title('Density plot of data')
#plt.sca(a[1])
#plt.hist2d(rl_scatter_model_fe['Fe_H'],model_mg_random-rl_scatter_model_fe['Fe_H'],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
#plt.xlabel('[Fe/H]')
#plt.title('Density plot of painted data')
#f.suptitle('Comparison between observed and model predictions, w/ RM, w/Mg Scatter')
#plt.show()


