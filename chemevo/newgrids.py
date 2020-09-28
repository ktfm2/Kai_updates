#Script to create different hdf5 files for plotting different scenarios
import numpy as np
import matplotlib.pyplot as plt
import math
import h5py
import json
import pandas as pd
import subprocess
import os
import sys

        
#No Merger Scenario
json_file = open("./params/newgridparams.json", "r")
json_object = json.load(json_file)
json_file.close()

json_object["flows"]["gasdump"]["mass"] = 0.0
json_object["flows"]["gasdump"]["central_radius"] = 0.0
json_object["flows"]["gasdump"]["radial_width"] = 0.1 #Avoid singularity

json_object["flows"]["alternategasdump"]["mass"] = 0.0
json_object["flows"]["alternategasdump"]["central_radius"] = 0.0
json_object["flows"]["alternategasdump"]["radial_width"] = 0.1



json_file = open("./params/newgridparams.json", "w")
json.dump(json_object, json_file)
json_file.close()


#Create new hdf5 file with this model
subprocess.run('./mains/run.exe nomerger.hdf5 ./params/newgridparams.json', shell=True)


#==============================

#Exp profile from centre Scenario
json_file = open("./params/newgridparams.json", "r")
json_object = json.load(json_file)
json_file.close()

json_object["flows"]["gasdump"]["mass"] = 5000.0
json_object["flows"]["gasdump"]["central_radius"] = 0.0
json_object["flows"]["gasdump"]["radial_width"] = 8.0

json_object["flows"]["alternategasdump"]["mass"] = 0.0
json_object["flows"]["alternategasdump"]["central_radius"] = 0.0
json_object["flows"]["alternategasdump"]["radial_width"] = 0.1



json_file = open("./params/newgridparams.json", "w")
json.dump(json_object, json_file)
json_file.close()


#Create new hdf5 file with this model
subprocess.run('./mains/run.exe centremerger.hdf5 ./params/newgridparams.json', shell=True)



#==============================

#Gaussian around Solar Scenario
json_file = open("./params/newgridparams.json", "r")
json_object = json.load(json_file)
json_file.close()

json_object["flows"]["gasdump"]["mass"] = 5000.0
json_object["flows"]["gasdump"]["central_radius"] = 8.1
json_object["flows"]["gasdump"]["radial_width"] = 1.0

json_object["flows"]["alternategasdump"]["mass"] = 0.0
json_object["flows"]["alternategasdump"]["central_radius"] = 0.0
json_object["flows"]["alternategasdump"]["radial_width"] = 0.1



json_file = open("./params/newgridparams.json", "w")
json.dump(json_object, json_file)
json_file.close()


#Create new hdf5 file with this model
subprocess.run('./mains/run.exe solargaussian.hdf5 ./params/newgridparams.json', shell=True)


#==============================

#Double Gaussian Scenario
json_file = open("./params/newgridparams.json", "r")
json_object = json.load(json_file)
json_file.close()

json_object["flows"]["gasdump"]["mass"] = 2500.0
json_object["flows"]["gasdump"]["central_radius"] = 8.1
json_object["flows"]["gasdump"]["radial_width"] = 1.0

json_object["flows"]["alternategasdump"]["mass"] = 2500.0
json_object["flows"]["alternategasdump"]["central_radius"] = 4.0
json_object["flows"]["alternategasdump"]["radial_width"] = 1.0



json_file = open("./params/newgridparams.json", "w")
json.dump(json_object, json_file)
json_file.close()


#Create new hdf5 file with this model
subprocess.run('./mains/run.exe doublemerger.hdf5 ./params/newgridparams.json', shell=True)

















