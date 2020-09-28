#Display the data from the grid search as a contour plot/print out values
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#df = pd.read_pickle("./search.pkl")
#df = pd.read_pickle("./mock_test.pkl")
#print(df.loc[:])



#df1 = pd.read_pickle("./search1.pkl")
#print(df1.loc[:])
#df2 = pd.read_pickle("./search2.pkl")
#print(df2.loc[:])
#df3 = pd.read_pickle("./search3.pkl")
#print(df3.loc[:])
#df4 = pd.read_pickle("./search4.pkl")
#print(df4.loc[:])


df1 = pd.read_pickle("./search3.pkl")


#print(df1['Metric'].type())
#print(df1.dtypes)

df1 = df1.sort_index()

x = df1['Time']
y = df1['Mass']
z = df1['Metric']

#print(x)

X = np.ones(11) #for the time axis
Y = np.ones(16) #for the mass axis

rows, cols = (16,11);

Z = [[1 for i in range(cols)] for j in range(rows)];


for i in range(11): #For the times
    for j in range(16): #For the masses
        X[i] = x[(16*i)];
        Y[j] = y[j];
        Z[j][i] = z[(16*i)+j];

plt.contourf(X,Y,Z,40,cmap='RdGy')
plt.xlabel('Time after birth of MW (Gyr)')
plt.ylabel('Total gas mass (x10^6 Solar Masses)')
plt.title('Contour plot of Chi^2 - SD18')
#plt.colorbar(label='Chi^2')
plt.savefig('/home/ktfm2/Documents/Project_Images/ForProject/GridSearchContourSD18.pdf', bbox_inches='tight')
plt.show()



