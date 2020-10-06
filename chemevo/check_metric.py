#Create various plots of chem evo models against data
import numpy as np
import pandas as pd
import math
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import matplotlib.colors as colors
import sys
sys.path.append('./scripts/')
from chemevo import *

#fl = chem_evo_data('./search.hdf5')
#fl = chem_evo_data('./multioutput.hdf5')
fl = chem_evo_data('./output.hdf5')


#astroNN VAC file from APOGEE DR16
data_file_1 = '/data/ktfm2/apogee_data/apogee_astroNN_DR16.fits'

hdu_list_1 = fits.open(data_file_1, memmap=True)
apogee_data = Table(hdu_list_1[1].data)

def betw(x,l,u):
    return (x>l)&(x<u)

def outs(x,l,u):
    return (x<l)|(x>u)


#Radial Migration filter
Remove_nans = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['age_lowess_correct']))&(apogee_data['age_lowess_correct']>0.0)&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(~pd.isna(apogee_data['FE_H_ERR']))&(apogee_data['LOGG']<3.5)&(outs(apogee_data['GALZ'],-1.0,1.0))&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],7.6,8.6))



#=====================================================================================================================

#Radial migration, randomness from gaussian
rl_random = np.random.normal(0.0, 4.0*((apogee_data['age_lowess_correct'][Remove_nans]/13.7)**0.5))
rl_scatter = apogee_data['rl'][Remove_nans]+rl_random

#====================================================================================================================

#Radial migration painting
rl_scatter_model_fe = fl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][Remove_nans]),['Fe'])
rl_scatter_model_mg = fl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][Remove_nans]),['Mg']) 


#=======================================================================================================================

#Plot radial migration - no scattering abundances
plt.scatter(apogee_data['FE_H'][Remove_nans],apogee_data['MG_H'][Remove_nans]-apogee_data['FE_H'][Remove_nans],s=2.0,color='b')
plt.scatter(rl_scatter_model_fe['Fe_H'],rl_scatter_model_mg['Mg_H']-rl_scatter_model_fe['Fe_H'],s=2.0,color='r')
plt.title('Radial Migration Effects')
plt.xlabel('[Fe/H]')
plt.ylabel('[Mg/Fe]')
plt.show()

#=======================================================================================================================

#Changed for with radial effects set
model_fe_uncert = np.random.normal(0.0, apogee_data['FE_H_ERR'][Remove_nans])
model_mg_uncert = np.random.normal(0.0, apogee_data['MG_H_ERR'][Remove_nans])

#Changed for radial effects data set
model_fe_random = rl_scatter_model_fe['Fe_H'] + model_fe_uncert
model_mg_random = rl_scatter_model_mg['Mg_H'] + model_mg_uncert


#==============================================================================================================

#Scatter both 
plt.scatter(apogee_data['FE_H'][Remove_nans],apogee_data['MG_H'][Remove_nans]-apogee_data['FE_H'][Remove_nans],s=2.0, color='b')
plt.scatter(model_fe_random,model_mg_random-model_fe_random,s=2.0,color='r')
plt.title('Scatter both Fe and Mg')
plt.xlabel('[Fe/H]')
plt.ylabel('[Mg/Fe]')
plt.show()

#=============================================================================================================

#Scatter Only Fe
plt.scatter(apogee_data['FE_H'][Remove_nans],apogee_data['MG_H'][Remove_nans]-apogee_data['FE_H'][Remove_nans],s=2.0, color='b')
plt.scatter(model_fe_random,rl_scatter_model_mg['Mg_H']-model_fe_random,s=2.0,color='r')
plt.title('Scatter only Fe, but in both axes')
plt.xlabel('[Fe/H]')
plt.ylabel('[Mg/Fe]')
plt.show()

#==============================================================================================================

#Scattering Mg only
plt.scatter(apogee_data['FE_H'][Remove_nans],apogee_data['MG_H'][Remove_nans]-apogee_data['FE_H'][Remove_nans],s=2.0, color='b')
plt.scatter(rl_scatter_model_fe['Fe_H'],model_mg_random-rl_scatter_model_fe['Fe_H'],s=2.0,color='r')
plt.title('Scatter only Mg, only affects y-axis')
plt.xlabel('[Fe/H]')
plt.ylabel('[Mg/Fe]')
plt.show()

#==============================================================================================================

#Density plots with Radial Migration but no scatter
f,a=plt.subplots(1,2,figsize=[15.,3.],sharex=True,sharey=True)
plt.sca(a[0])
plt.hist2d(apogee_data['FE_H'][Remove_nans],apogee_data['MG_H'][Remove_nans]-apogee_data['FE_H'][Remove_nans],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]')
plt.ylabel('[Mg/H]')
plt.title('Density plot of data')
plt.sca(a[1])
plt.hist2d(rl_scatter_model_fe['Fe_H'],rl_scatter_model_mg['Mg_H']-rl_scatter_model_fe['Fe_H'],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]')
plt.title('Density plot of painted data')
f.suptitle('Comparison between observed and model predictions, w/ RM, wo/Scatter')
plt.show()

#===============================================================================================================

#Density plots with Radial Migration and scatter in Fe
f,a=plt.subplots(1,2,figsize=[15.,3.],sharex=True,sharey=True)
plt.sca(a[0])
plt.hist2d(apogee_data['FE_H'][Remove_nans],apogee_data['MG_H'][Remove_nans]-apogee_data['FE_H'][Remove_nans],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]')
plt.ylabel('[Mg/H]')
plt.title('Density plot of data')
plt.sca(a[1])
plt.hist2d(model_fe_random,rl_scatter_model_mg['Mg_H']-model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]')
plt.title('Density plot of painted data')
f.suptitle('Comparison between observed and model predictions, w/ RM, w/Fe Scatter')
plt.show()

#===============================================================================================================

#Density plots with Radial Migration and scatter in both
f,a=plt.subplots(1,2,figsize=[15.,3.],sharex=True,sharey=True)
plt.sca(a[0])
plt.hist2d(apogee_data['FE_H'][Remove_nans],apogee_data['MG_H'][Remove_nans]-apogee_data['FE_H'][Remove_nans],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]')
plt.ylabel('[Mg/H]')
plt.title('Density plot of data')
plt.sca(a[1])
plt.hist2d(model_fe_random,model_mg_random-model_fe_random,bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]')
plt.title('Density plot of painted data')
f.suptitle('Comparison between observed and model predictions, w/ RM, w/Both Scatter')
plt.show()

#==============================================================================================================

#Density plots with Radial Migration and scatter in Mg
f,a=plt.subplots(1,2,figsize=[15.,3.],sharex=True,sharey=True)
plt.sca(a[0])
plt.hist2d(apogee_data['FE_H'][Remove_nans],apogee_data['MG_H'][Remove_nans]-apogee_data['FE_H'][Remove_nans],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]')
plt.ylabel('[Mg/H]')
plt.title('Density plot of data')
plt.sca(a[1])
plt.hist2d(rl_scatter_model_fe['Fe_H'],model_mg_random-rl_scatter_model_fe['Fe_H'],bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.xlabel('[Fe/H]')
plt.title('Density plot of painted data')
f.suptitle('Comparison between observed and model predictions, w/ RM, w/Mg Scatter')
plt.show()


