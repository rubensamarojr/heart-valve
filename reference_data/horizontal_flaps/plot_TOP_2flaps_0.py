# Plot top leaflet displacements
# Gil, et al 2010. The Immersed Structural Potential Method for haemodynamic applications
# Example II - Etop = Ebottom = 5.6 Mpa

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
tmax = 3.0
xmin = 0.0
xmax = 0.007
ymin = 0.006
ymax = 0.0

# Font sizes
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14
SPH_LINE = 2

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Read files data
data0 = np.genfromtxt("TopLeafletObserver_Position_0000000000.dat")
data1x = np.genfromtxt("output_2_flaps_2mm_0_2015_x.tsv")
data1y = np.genfromtxt("output_2_flaps_2mm_0_2015_y.tsv")
data2x = np.genfromtxt("output_2_flaps_2mm_0_2016_x.tsv")
data2y = np.genfromtxt("output_2_flaps_2mm_0_2016_y.tsv")
# Time
t0 = data0[:,0]
t1x = data1x[:,0]
t1y = data1y[:,0]
t2x = data2x[:,0]
t2y = data2y[:,0]
# Displacements
x0 = data0[:,1] - data0[1,1]
y0 = data0[:,2] - data0[1,2]
x1 = data1x[:,1]/100
y1 = data1y[:,1]/100
x2 = data2x[:,1]/100
y2 = data2y[:,1]/100

# Plotting the X displacement
plt.plot(t1x, x1, color='C2', linestyle=':', marker='+', label='Immersogeometric Kamensky et al. 2015')
plt.plot(t2x, x2, color='C4', linestyle=':', marker='*', label='B-Spline mesh Kadapa et al. 2016')
plt.plot(t0, x0, color='C1', linewidth=SPH_LINE, label='SPH Linear Elastic Isotropic material')
plt.xlim([tmin, tmax])
plt.ylim([xmin, xmax])
plt.xlabel('Time (s)')
plt.ylabel('X-displacement (m)')
legend = plt.legend(loc='best', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('C0') # Put a nicer background color on the legend
fig = plt.gcf() # get current figure
fig.savefig('x-disp-top-2-flaps-2mm-0-Linear.png')
plt.show()

# Plotting the Y displacement
plt.plot(t1y, y1, color='C2', linestyle=':', marker='+', label='Immersogeometric Kamensky et al. 2015')
plt.plot(t2y, y2, color='C4', linestyle=':', marker='*', label='B-Spline mesh Kadapa et al. 2016')
plt.plot(t0, y0, color='C1', linewidth=SPH_LINE, label='SPH Linear Elastic Isotropic material')
plt.xlim([tmin, tmax])
plt.ylim([ymin, ymax])
plt.xlabel('Time (s)')
plt.ylabel('Y-displacement (m)')
legend = plt.legend(loc='best', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('C0') # Put a nicer background color on the legend
fig = plt.gcf() # get current figure
fig.savefig('y-disp-top-2-flaps-2mm-0-Linear.png')
plt.show()