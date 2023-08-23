# Plot bottom leaflet displacements
# Gil, et al 2010. The Immersed Structural Potential Method for haemodynamic applications
# Example III - Etop = 5.6 Mpa, Ebottom = 11.2 Mpa

# importing the required modules
import matplotlib.pyplot as plt
import numpy as np
import math

plt.style.use('seaborn-whitegrid')
# seaborn-whitegrid
# seaborn-deep
# seaborn-dark-palette
# seaborn-colorblind
# fivethirtyeight
# ggplot
# seaborn
# seaborn-pastel
# seaborn-white
# tableau-colorblind10
# Solarize_Light2

# User set
tmin = 0.0
tmax = 2.0
xmin = 0.0
xmax = 0.007
ymin = -0.006
ymax = 0.0

# Discretization
dx = 1e-4

# Font sizes
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Read file data
data = np.genfromtxt("../reference_data/horizontal_flaps/BottomLeafletObserver_Position_0000000000.dat")
# Time
t = data[:,0]
# Displacements
x = data[:,1] - data[1,1]
y = data[:,2] - data[1,2]

# Plotting the X displacement
plt.plot(t, x, color='C0')
plt.xlim([tmin, tmax])
plt.ylim([xmin, xmax])
plt.xlabel('Time (s)')
plt.ylabel('X-displacement (m)')
fig = plt.gcf() # get current figure
fig.savefig('x-disp-bottom-2-flaps-2mm-1.png')
plt.show()

# Plotting the Y displacement
plt.plot(t, y, color='C1')
plt.xlim([tmin, tmax])
plt.ylim([ymin, ymax])
plt.xlabel('Time (s)')
plt.ylabel('Y-displacement (m)')
fig = plt.gcf() # get current figure
fig.savefig('y-disp-bottom-2-flaps-2mm-1.png')
plt.show()