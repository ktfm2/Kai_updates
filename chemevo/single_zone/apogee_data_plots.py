#Use to plot models on top of data
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
import math
from matplotlib.colors import PowerNorm
import matplotlib.colors as colors
import pandas as pd
import sys
from scipy.interpolate import RectBivariateSpline, CubicSpline
sys.path.append('../scripts/')
from chemevo import *
hl = chem_evo_data('../output.hdf5')
fl = chem_evo_data('../comparison.hdf5')
#fl = chem_evo_data('../multioutput.hdf5')
gl = chem_evo_data('../comparison.hdf5')

#plt.plot(fl.t,fl.SFR)
#plt.show()

data_file_1 = '/data/ktfm2/apogee_data/apogee_astroNN_DR16.fits'  #The vac file
data_file_2 = '/data/ktfm2/apogee_data/allStar_r12_l33.fits'  #The all star file from apogeee

hdu_list_1 = fits.open(data_file_1, memmap=True)  #Open the fits file
apogee_data = Table(hdu_list_1[1].data) #Creates table from the fits file
#print(apogee_data.colnames) #Prints column names
#print(apogee_data['apogee_id','GALR','MG_H','FE_H','e','GALZ']) #Prints columns
#hdu_list_1.info() #Overview of file, 473307 rows 
#print(hdu_list_1[1].columns) #Prints column details, including name and format 

#hdu_list_2 = fits.open(data_file_2, memmap=True)
#hdu_list_2.info() #473307 rows - match!
#print(hdu_list_2[1].columns)

apogee_data.sort(['e'])

def betw(x,l,u):
	return (x>l)&(x<u)

def outs(x,l,u):
	return (x<l)|(x>u)

#Create individual plot at solar radius
#fltr = (~pd.isna(apogee_data['GALR']))&(~pd.isna(apogee_data['GALZ']))&(~pd.isna(apogee_data['e']))&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['FE_H_ERR']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(apogee_data['LOGG']<3.5)&outs(apogee_data['GALZ'],-1.0,1.0)&(apogee_data['FE_H_ERR']<0.2)&betw(apogee_data['GALR'],8.0,8.2)


#Solar radius with guiding radius
#fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['GALZ']))&(~pd.isna(apogee_data['e']))&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['FE_H_ERR']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(apogee_data['LOGG']<3.5)&outs(apogee_data['GALZ'],-1.0,1.0)&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&betw(apogee_data['rl'],7.6,8.6)
fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['GALZ']))&(~pd.isna(apogee_data['e']))&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['FE_H_ERR']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(apogee_data['LOGG']<3.5)&(betw(apogee_data['GALZ'],-5.0,5.0))&(outs(apogee_data['GALZ'],-1.0,1.0))&(apogee_data['FE_H_ERR']<0.2)&betw(apogee_data['rl'],7.6,8.6)


#==========================================================================================================================

c=plt.scatter([0.,1.],[0.,1.],c=[0.,1.])
plt.clf()


#Set up to plot [alpha/Fe] against [Fe/H] from model at above radius
radius=8.1;
dat = fl.abund['Fe']
t = np.linspace(fl.t[0],fl.t[-1],100)
rbs = RectBivariateSpline(fl.R,fl.t,dat)
a = rbs(radius,t)
dat=fl.abund['Mg']-fl.abund['Fe']
rbs = RectBivariateSpline(fl.R,fl.t,dat)
b = rbs(radius,t)

#Paint on time stamps
#timestamp =[1.0];
#radiusstamp =[8.1];
timestamp = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0];
radiusstamp= [8.1,8.1,8.1,8.1,8.1,8.1,8.1,8.1,8.1,8.1,8.1,8.1,8.1];
time_fe = fl.paint(radiusstamp,timestamp,['Fe'])
time_mg = fl.paint(radiusstamp,timestamp,['Mg'])
	
#plt.scatter(apogee_data['FE_H'][fltr],(apogee_data['MG_H'][fltr] - apogee_data['FE_H'][fltr]),c=plt.cm.viridis(apogee_data['e'][fltr]),s=2.0,zorder=-1);
#plt.plot(a.T,b.T,color='black',linewidth=3.0,zorder=0, label='2.7Gyr, 6.4x10^9')
plt.hist2d(apogee_data['FE_H'][fltr],
           (apogee_data['MG_H']-apogee_data['FE_H'])[fltr],
          bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
plt.plot(a.T,b.T,color='black',linewidth=3.0, label='2.7Gyr, 6.4x10^9')
#plt.scatter(time_fe['Fe_H'],time_mg['Mg_H']-time_fe['Fe_H'],s=8.0,color='r',zorder=1)
plt.title('Model against data, radius = 8.1 kpc')
#plt.colorbar(c,label='Eccentricity')
plt.legend()
plt.xlim(-1.5,0.6) #Usually -2.5 but -1.5 to cut off initial rise
plt.ylim(-0.1,0.4)
plt.xlabel('[Fe/H]')
plt.ylabel('[Mg/Fe]')
plt.savefig('../../../../Project_Images/ForProject/GridSearchResult.pdf', bbox_inches='tight')
plt.show()


#plt.plot(gl.t,gl.SFR)
#plt.xlabel(r'$t/\,\mathrm{Gyr}$')
#plt.ylabel(r'Rate /$\,\mathrm{M}_\odot\,\mathrm{pc}^{-2}\,\mathrm{Gyr}^{-1}$')
#plt.ylim(0.,13.0)
#plt.title('SFR for gas dump scenario')
#plt.show()

#=========================================================================================================================

#Density Plots
ssM=[np.argmin(np.abs(fl.R-radii)) for radii in [6.1,8.1,10.1]]
ssG=[np.argmin(np.abs(gl.R-radii)) for radii in [6.1,8.1,10.1]]
ssC=[np.argmin(np.abs(hl.R-radii)) for radii in [6.1,8.1,10.1]]
rranges = [[5.6,6.6],[7.6,8.6],[9.6,10.6]]
#rranges = [[5,7],[7,9],[9,11]]

#Selection for GALR
NewRadii = [6.1,8.1,10.1]
#selectionR = [
#        (apogee_data['FE_H_ERR']<0.2)&(apogee_data['LOGG']<3.5)&betw(apogee_data['GALR'],rd,ru)&(np.abs(apogee_data['GALZ'])>1) for rd,ru in rranges
#]

#Selection with guiding radius
#selectionR = [
#        (apogee_data['FE_H_ERR']<0.2)&(apogee_data['LOGG']<3.5)&betw(apogee_data['rl'],rd,ru)&(np.abs(apogee_data['GALZ'])>1.0)&(np.abs(apogee_data['GALZ'])<5.0) for rd,ru in rranges
#]
selectionR = [
        (apogee_data['FE_H_ERR']<0.2)&(apogee_data['LOGG']<3.5)&betw(apogee_data['rl'],rd,ru)&(np.abs(apogee_data['GALZ']>1.0))&(np.abs(apogee_data['GALZ'])<5.0) for rd,ru in rranges
]


f,a=plt.subplots(1,3,figsize=[15.,3.],sharex=True,sharey=True)
e='Mg'
plt.subplots_adjust(wspace=0.05)
for rr in range(3):
    plt.sca(a[rr])
    plt.plot(hl.abund['Fe'][ssC[rr]]-hl.abund['H'][ssC[rr]],
                 (hl.abund[e]-hl.abund['Fe'])[ssC[rr]],c='r',lw=3, label='2.7Gyr, 6.4x10^9')
#    plt.plot(gl.abund['Fe'][ssG[rr]]-gl.abund['H'][ssG[rr]],
#                 (gl.abund[e]-gl.abund['Fe'])[ssG[rr]],c='k',lw=3, label='Fiducial model')
#    plt.plot(fl.abund['Fe'][ssM[rr]]-fl.abund['H'][ssM[rr]],
#                 (fl.abund[e]-fl.abund['Fe'])[ssM[rr]],c='b',lw=3, label='10^10')
    plt.hist2d(apogee_data['FE_H'][selectionR[rr]],
               (apogee_data['%s_H'%e.upper()]-apogee_data['FE_H'])[selectionR[rr]],
              bins=50,range=[[-1.5,0.6],[-0.1,0.4]],norm=colors.LogNorm(),cmap=plt.cm.Spectral_r);
    plt.xlim(-1.5,0.5)
    plt.ylim(-0.1,0.4)
    plt.xlabel(r'$\mathrm{[Fe/H]}$', fontsize=12)
    if rr==0:
        plt.ylabel(r'$\mathrm{[%s/Fe]}$'%e, fontsize=12)
    plt.title('$%.1f<R/\mathrm{kpc}<%.1f$'%(rranges[rr][0],rranges[rr][1]),fontsize=12)
#    plt.title('Galactocentric Radius $R/\mathrm{kpc}=%.1f$'%(NewRadii[rr]),fontsize=12)
#    plt.legend(loc='upper right')
#    f.suptitle('Different models against data', fontsize=12)
f.figsize = [8.0,6.0]
plt.legend()
#f.savefig('../../../../Project_Images/ForProject/VaryMassDensity.pdf', bbox_inches='tight')
#f.savefig('../../../../Project_Images/ForProject/GasNoGas.pdf', bbox_inches='tight')
#f.savefig('../../../../Project_Images/ForProject/GridSearchResult.pdf', bbox_inches='tight')
plt.show()




#========================================================================================================================

radii = [6.1,8.1,10.1]; #Can switch to whatever
#lower_radii = [6.0,8.0,10.0];
#upper_radii = [6.2,8.2,10.2];

lower_radii = [5.6,7.6,9.6];
upper_radii = [6.6,8.6,10.6];


f,ax = plt.subplots(1,3,figsize=[15.,3.], sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.05)

for x in range(3):
	radius =radii[x]; #Radius for the models
	lower_radius = lower_radii[x]; #Lower radius for data
	upper_radius = upper_radii[x]; #Upper radius for data
#Filter with GALR
#	fltr = (~pd.isna(apogee_data['GALR']))&(~pd.isna(apogee_data['GALZ']))&(~pd.isna(apogee_data['e']))&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['FE_H_ERR']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(apogee_data['LOGG']<3.5)&outs(apogee_data['GALZ'],-1.0,1.0)&(apogee_data['FE_H_ERR']<0.2)&betw(apogee_data['GALR'],lower_radius,upper_radius)

#Filter with rl	
	#fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['GALZ']))&(~pd.isna(apogee_data['e']))&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['FE_H_ERR']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(apogee_data['LOGG']<3.5)&outs(apogee_data['GALZ'],-1.0,1.0)&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&betw(apogee_data['rl'],lower_radius,upper_radius)
	fltr = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['GALZ']))&(~pd.isna(apogee_data['e']))&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['FE_H_ERR']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['LOGG']))&(apogee_data['LOGG']<3.5)&(betw(apogee_data['GALZ'],-5.0,5.0))&(apogee_data['FE_H_ERR']<0.2)&betw(apogee_data['rl'],lower_radius,upper_radius)

	counter = len(apogee_data['FE_H_ERR'][fltr]);
	
	plt.sca(ax[x])
	#plt.subplot(1,3,(x+1))
	c=plt.scatter([0.,1.],[0.,1.],c=[0.,1.])
	#plt.clf()


	#Set up to plot [alpha/Fe] against [Fe/H] from model at above radius
	dat = fl.abund['Fe']
	t = np.linspace(fl.t[0],fl.t[-1],100)
	rbs = RectBivariateSpline(fl.R,fl.t,dat)
	a = rbs(radius,t)
	dat=fl.abund['Mg']-fl.abund['Fe']
	rbs = RectBivariateSpline(fl.R,fl.t,dat)
	b = rbs(radius,t)
	
	#alpha_obs = ((apogee_data['MG_H'][x])*pow(apogee_data['MG_H_ERR'][x],-2.0) + (apogee_data['CA_H'][x])*pow(apogee_data['CA_H_ERR'][x],-2.0) + (apogee_data['O_H'][x])*pow(apogee_data['O_H_ERR'][x],-2.0) + (apogee_data['SI_H'][x])*pow(apogee_data['SI_H_ERR'][x],-2.0) + (apogee_data['S_H'][x])*pow(apogee_data['S_H_ERR'][x],-2.0))/(pow(apogee_data['MG_H_ERR'][x],-2.0) + pow(apogee_data['CA_H_ERR'][x],-2.0) + pow(apogee_data['O_H_ERR'][x],-2.0) + pow(apogee_data['SI_H_ERR'][x],-2.0) + pow(apogee_data['S_H_ERR'][x],-2.0));   
	plt.scatter(apogee_data['FE_H'][fltr],(apogee_data['MG_H'][fltr] - apogee_data['FE_H'][fltr]),c=plt.cm.viridis(apogee_data['e'][fltr]),s=2.0);
	plt.plot(a.T,b.T,color='black',linewidth=3.0)

	plt.title('%i APOGEE giants, r = %.1f kpc' %(counter, round(radius,1)))
	if x == 2:
		plt.colorbar(c,label='Eccentricity')
	plt.xlim(-1.6,0.5) #Usually -2.5 but -1.6 to cut off initial rise
	plt.ylim(-0.2,0.5)
	plt.xlabel('[Fe/H]')
	if x == 0:
		plt.ylabel('[Mg/Fe]')
	#plt.show()
#f.constrained_layout()	
#f.suptitle('Models against data, annulus width of 0.2 kpc, 1<|Z|<5, logg<3.5, EFeH<0.2')
f.figsize = [8.0,6.0]
#f.savefig('testfig1.pdf', bbox_inches='tight')
plt.show()

hdu_list_1.close()
#hdu_list_2.close()
