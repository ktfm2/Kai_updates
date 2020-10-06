import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../scripts/')
from chemevo import *

#f = chem_evo_data('../output.hdf5')
#g = chem_evo_data('../comparison.hdf5')
#h = chem_evo_data('../multioutput.hdf5')
#sg = chem_evo_data('../Sausage.hdf5')

f = chem_evo_data('../KSsfr.hdf5')

#f.plot_elements_against_time('Al','Mg',el_u='Fe',el_u2='Mn',radius=8.0,color='k')
#f.plot_elements_against_time('Al','Mg',el_u='Fe',el_u2='Mn',radius=4.0,color='r')
#g.plot_elements_against_time('Al','Mg',el_u='Fe',el_u2='Mn',radius=8.0,color='y')
#g.plot_elements_against_time('Al','Mg',el_u='Fe',el_u2='Mn',radius=4.0,color='g')
#sg.plot_elements_against_time('Al','Mg',el_u='Fe',el_u2='Mn',radius=4.0,color='b')
#plt.show(block=True)

f.summary_plot()
plt.show()

f.plot_elements_against_time('Fe','Mg',el_u2='Fe')
plt.show(block=True)

f.plot_elements_against_time('Fe','Mg',el_u2='Fe', radius=3.0)
plt.show(block=True)
#f.summary_plot()
#plt.show()

f.plot_radial(el='Mg', el2='Fe')
plt.show()

f.plot_radial(el='Fe')
plt.show()

f.plot_time_range(el='Mg', el2='Fe')
plt.show()

f.plot_time_range(el='Fe')
plt.show()

plt.plot(f.t,f.Inflow,color='k')
plt.plot(f.t,f.abund['rSFR'][np.argmin(np.abs(f.R-8.3))],color='b')
plt.plot(f.t,f.SFR,color='r')
plt.show()


plt.plot(f.Mgas.T[np.argmin(np.abs(f.R-8.3))],f.SFR,color='r')
xx=np.linspace(1,100,100)
plt.plot(xx,0.067*(xx**1.4),color='k')
plt.xscale('log')
plt.yscale('log')
#plt.scatter(f.t,timespots,s=3.0)
#plt.ylim(0.0,15.0);
#plt.xlim(1.5,4.5);
#plt.plot(f.t,f.Mgas)
#plt.plot(f.t,f.SNII)
#plt.xlabel('Time (Gyr)')
#plt.ylabel('Rate (Solar mass pc^-2 Gyr^-1)')
#plt.title('SFR')
#plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/SFR.pdf')
plt.show()
