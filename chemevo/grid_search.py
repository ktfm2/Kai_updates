#Script to do a grid search of gas dump mass and gas dump time
#Compares against 4 different sets of ages - linear correct form astroNN; lowess correct from astroNN; Sanders & Das; APOKASC
import numpy as np
import matplotlib.pyplot as plt
import math
import h5py
import json
from astropy.io import fits
from astropy.table import Table, join
import pandas as pd
import subprocess
import os
import sys
sys.path.append('./scripts/')
from chemevo import *


data_file_1 = '/data/ktfm2/apogee_data/apogee_astroNN_DR16.fits' #The astroNN VAC for APOGEE DR16
hdf5_file = '/data/ktfm2/apogee_data/gaia_spectro.hdf5' #The hdf5 file for Sanders and Das
data_file_2 = '/data/jls/apokasc_astroNN.fits' #The APOKASC data file joined with AstroNN

hdf = h5py.File(hdf5_file, "r")
dataset = hdf['data']

log_age_data = dataset["log10_age"]
ID_data = dataset["APOGEE_ID"]

SD_table = Table([ID_data, log_age_data], names=('apogee_id', 'log_age_data'))

hdu_list_1 = fits.open(data_file_1, memmap=True) #open fits file
apogee_data = Table(hdu_list_1[1].data) #Creates table from fits file
hdu_list_1.close() #Close the fits file

hdu_list_2 = fits.open(data_file_2, memmap=True) #open fits file
apokasc_data = Table(hdu_list_2[1].data) #Creates table from fits file
hdu_list_2.close() #Close the fits file


#Join tables together
full_table = join(apogee_data, SD_table)

#Define functions for the filter
def betw(x,l,u):
    return (x>l)&(x<u)

def outs(x,l,u):
    return (x<l)|(x>u)

#Define filter for apogee data, use guiding centre radius RL, galactic height GALZ, surface gravity LOGG
#Have 4 different filters and so on for linear age, lowess age, S&D age, APOKASC age - this extends to have disc stars

NaN_bit1 = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['age_lowess_correct']))&(~pd.isna(apogee_data['GALZ']))&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['FE_H_ERR']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['MG_H_ERR']))&(~pd.isna(apogee_data['LOGG'])) 

fltr1 = NaN_bit1&(apogee_data['age_lowess_correct']>0.0)&(apogee_data['LOGG']<3.5)&(betw(apogee_data['GALZ'],-5.0,5.0))&(outs(apogee_data['GALZ'],-1.0,1.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],7.6,8.6))


NaN_bit2 = (~pd.isna(apogee_data['rl']))&(~pd.isna(apogee_data['age_linear_correct']))&(~pd.isna(apogee_data['GALZ']))&(~pd.isna(apogee_data['FE_H']))&(~pd.isna(apogee_data['FE_H_ERR']))&(~pd.isna(apogee_data['MG_H']))&(~pd.isna(apogee_data['MG_H_ERR']))&(~pd.isna(apogee_data['LOGG'])) 

fltr2 = NaN_bit2&(apogee_data['age_linear_correct']>0.0)&(apogee_data['LOGG']<3.5)&(betw(apogee_data['GALZ'],-5.0,5.0))&(outs(apogee_data['GALZ'],-1.0,1.0))&(apogee_data['FE_H_ERR']<0.2)&(betw(apogee_data['rl'],7.6,8.6))


NaN_bit3 = (~pd.isna(full_table['rl']))&(~pd.isna(full_table['log_age_data']))&(~pd.isna(full_table['GALZ']))&(~pd.isna(full_table['FE_H']))&(~pd.isna(full_table['FE_H_ERR']))&(~pd.isna(full_table['MG_H']))&(~pd.isna(full_table['MG_H_ERR']))&(~pd.isna(full_table['LOGG'])) 

fltr3 = NaN_bit3&(full_table['LOGG']<3.5)&(betw(full_table['GALZ'],-5.0,5.0))&(outs(full_table['GALZ'],-1.0,1.0))&(full_table['FE_H_ERR']<0.2)&(betw(full_table['rl'],7.6,8.6))

NaN_bit4 = (~pd.isna(apokasc_data['rl']))&(~pd.isna(apokasc_data['LogAge']))&(~pd.isna(apokasc_data['GALZ']))&(~pd.isna(apokasc_data['FE_H']))&(~pd.isna(apokasc_data['FE_H_ERR']))&(~pd.isna(apokasc_data['MG_H']))&(~pd.isna(apokasc_data['MG_H_ERR']))&(~pd.isna(apokasc_data['LOGG'])) 

fltr4 = NaN_bit4&(apokasc_data['LOGG']<3.5)&(betw(apokasc_data['GALZ'],-5.0,5.0))&(apokasc_data['FE_H_ERR']<0.2)&(betw(apokasc_data['rl'],7.6,8.6))



#Randomly scatter rl
rl_random1 = np.random.normal(0.0, 4.0*((apogee_data['age_lowess_correct'][fltr1]/13.7)**0.5));
rl_scatter1 = apogee_data['rl'][fltr1]+rl_random1; 

rl_random2 = np.random.normal(0.0, 4.0*((apogee_data['age_linear_correct'][fltr2]/13.7)**0.5));
rl_scatter2 = apogee_data['rl'][fltr2]+rl_random2; 

rl_random3 = np.random.normal(0.0, 4.0*(((10.0**full_table['log_age_data'][fltr3])/13.7)**0.5));
rl_scatter3 = full_table['rl'][fltr3]+rl_random3; 

rl_random4 = np.random.normal(0.0, 4.0*(((10.0**(apokasc_data['LogAge'][fltr4]-3.0))/13.7)**0.5));
rl_scatter4 = apokasc_data['rl'][fltr4]+rl_random4; 

#Set up mock data to test the program
#gl = chem_evo_data('./comparison.hdf5')
#mock_fe = gl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Fe'])
#mock_mg = gl.paint(rl_scatter,(13.7-apogee_data['age_lowess_correct'][fltr]),['Mg'])

#Define points on grid
ageranges = np.arange(2.7,5.8,0.3);
massranges = range(100,10601,700);

#Points on grid for mock test
#ageranges = [3.8, 4.0, 4.2, 4.4, 4.6];
#massranges = [4000.0, 4500.0, 5000.0, 5500.0, 6000.0];


#Set counter for storing results in table later
counter = 0; 

#Set up initial dataframe
df1 = pd.DataFrame()
df2 = pd.DataFrame()
df3 = pd.DataFrame()
df4 = pd.DataFrame()


#Double loop through gas dump time and gas dump mass
for age in ageranges:
    for mass in massranges:
        
        #Change JSON file to the age and density
        json_file = open("./params/searchparams.json", "r")
        json_object = json.load(json_file)
        json_file.close()

        json_object["flows"]["gasdump"]["mass"] = mass
        json_object["flows"]["gasdump"]["time"] = age

        json_file = open("./params/searchparams.json", "w")
        json.dump(json_object, json_file)
        json_file.close()


        #Create new hdf5 file with this model
        subprocess.run('./mains/run.exe search.hdf5 ./params/searchparams.json', shell=True)


        #Paint data using created model
        f = chem_evo_data('./search.hdf5')

        model_fe1 = f.paint(rl_scatter1,(13.7 - apogee_data['age_lowess_correct'][fltr1]),['Fe']);
        model_mg1 = f.paint(rl_scatter1,(13.7 - apogee_data['age_lowess_correct'][fltr1]),['Mg']);

        model_fe2 = f.paint(rl_scatter2,(13.7 - apogee_data['age_linear_correct'][fltr2]),['Fe']);
        model_mg2 = f.paint(rl_scatter2,(13.7 - apogee_data['age_linear_correct'][fltr2]),['Mg']);
        
        model_fe3 = f.paint(rl_scatter3,(13.7 - 10.0**full_table['log_age_data'][fltr3]),['Fe']);
        model_mg3 = f.paint(rl_scatter3,(13.7 - 10.0**full_table['log_age_data'][fltr3]),['Mg']);

        model_fe4 = f.paint(rl_scatter4,(13.7 - 10.0**(apokasc_data['LogAge'][fltr4]-3.0)),['Fe']);
        model_mg4 = f.paint(rl_scatter4,(13.7 - 10.0**(apokasc_data['LogAge'][fltr4]-3.0)),['Mg']);
 
        #Work out each metric for each combination
        metric1 = 0.0;
        metric1 += ((model_fe1['Fe_H'] - apogee_data['FE_H'][fltr1])**2.0)/((apogee_data['FE_H_ERR'][fltr1])**2.0); 
        metric1 += ((model_mg1['Mg_H'] - model_fe1['Fe_H'] - apogee_data['MG_H'][fltr1] + apogee_data['FE_H'][fltr1])**2.0)/((apogee_data['FE_H_ERR'][fltr1])**2.0 + (apogee_data['MG_H_ERR'][fltr1])**2.0); 

        metric2 = 0.0;
        metric2 += ((model_fe2['Fe_H'] - apogee_data['FE_H'][fltr2])**2.0)/((apogee_data['FE_H_ERR'][fltr2])**2.0); 
        metric2 += ((model_mg2['Mg_H'] - model_fe2['Fe_H'] - apogee_data['MG_H'][fltr2] + apogee_data['FE_H'][fltr2])**2.0)/((apogee_data['FE_H_ERR'][fltr2])**2.0 + (apogee_data['MG_H_ERR'][fltr2])**2.0); 

        metric3 = 0.0;
        metric3 += ((model_fe3['Fe_H'] - full_table['FE_H'][fltr3])**2.0)/((full_table['FE_H_ERR'][fltr3])**2.0); 
        metric3 += ((model_mg3['Mg_H'] - model_fe3['Fe_H'] - full_table['MG_H'][fltr3] + full_table['FE_H'][fltr3])**2.0)/((full_table['FE_H_ERR'][fltr3])**2.0 + (full_table['MG_H_ERR'][fltr3])**2.0); 

        metric4 = 0.0;
        metric4 += ((model_fe4['Fe_H'] - apokasc_data['FE_H'][fltr4])**2.0)/((apokasc_data['FE_H_ERR'][fltr4])**2.0); 
        metric4 += ((model_mg4['Mg_H'] - model_fe4['Fe_H'] - apokasc_data['MG_H'][fltr4] + apokasc_data['FE_H'][fltr4])**2.0)/((apokasc_data['FE_H_ERR'][fltr4])**2.0 + (apokasc_data['MG_H_ERR'][fltr4])**2.0); 

        #Work out metric for mock data test
        #metric += ((model_fe['Fe_H'] - mock_fe['Fe_H'])**2.0)/((apogee_data['FE_H_ERR'][fltr])**2.0); 
        #metric += ((model_mg['Mg_H'] - model_fe['Fe_H'] - mock_mg['Mg_H'] + mock_fe['Fe_H'])**2.0)/((apogee_data['FE_H_ERR'][fltr])**2.0 + (apogee_data['MG_H_ERR'][fltr])**2.0); 
        

        #Get the metric by taking sum of vector elements - metric is stored initially as vector
        single_metric1 = 0;
        single_metric1 = np.sum(metric1);
        single_metric2 = 0;
        single_metric2 = np.sum(metric2);
        single_metric3 = 0;
        single_metric3 = np.sum(metric3);
        single_metric4 = 0;
        single_metric4 = np.sum(metric4);

        #Store metric in file
        new_df1 = pd.DataFrame({"Mass": mass, "Metric": single_metric1, "Time": age}, index=[0])
        df1 = df1.append(new_df1, ignore_index = True)
        new_df2 = pd.DataFrame({"Mass": mass, "Metric": single_metric2, "Time": age}, index=[0])
        df2 = df2.append(new_df2, ignore_index = True)
        new_df3 = pd.DataFrame({"Mass": mass, "Metric": single_metric3, "Time": age}, index=[0])
        df3 = df3.append(new_df3, ignore_index = True)
        new_df4 = pd.DataFrame({"Mass": mass, "Metric": single_metric4, "Time": age}, index=[0])
        df4 = df4.append(new_df4, ignore_index = True)


        #Update counter for checking progress
        counter += 1;
        print(counter)


#Identify which scenario is best
df1 = df1.sort_values(by=['Metric'])
df1.to_pickle("./search1.pkl")
df2 = df2.sort_values(by=['Metric'])
df2.to_pickle("./search2.pkl")
df3 = df3.sort_values(by=['Metric'])
df3.to_pickle("./search3.pkl")
df4 = df4.sort_values(by=['Metric'])
df4.to_pickle("./search4.pkl")

#df.to_pickle("./mock_test.pkl")

