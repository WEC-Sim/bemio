"""
Created on Thu Jan 15 10:42:54 2015

@author: mlawson
"""
import aqwaio as aio
import hydroData as hd
import matplotlib.pyplot as plt
plt.close('all')
plt.interactive(True)

# Load AQWA output data file
aq = aio.AqwaOutput(outFile='./data/aqwa-example-data.lis')

# Plot diag components of added mass and damping
componentsToPlot = [[0,0],[1,1],[2,2]]
aq.data[0].plotAddedMassAndDamping(componentsToPlot)
aq.data[1].plotAddedMassAndDamping(componentsToPlot)
aq.data[2].plotAddedMassAndDamping(componentsToPlot)

# Write hydrodynamic data for WEC-Sim
# hd.writeWecSimHydroData(aq.data) # Not currently working in hydroData

# Write hydrodynamic data to HDF5 file format
hd.writeHdf5(aq.data,aq.files['hdf5'])
hd.writePickle(aq.data,aq.files['pickle'])