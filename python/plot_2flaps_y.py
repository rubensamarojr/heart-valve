# Plot top leaflet displacements
# Thomas Wick, 2014. Flapping and contact FSI computations with the fluid–solid
# interface-tracking/interface-capturing technique and mesh adaptivity
# Example 4.2 - Thin flap

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
tmax = 0.8
tmax = 0.8
xmin = 0.037
xmax = 0.040
ymin = -0.005
ymax = 0.005
flap_thickness = 0.065e-2
Umax = 0.0135
speed_sound = 10*Umax
shift_time = -0.0185/speed_sound
t0x = 0.165 # Point to ref shift line
y0x = 0.039 # Point to ref shift line
t0y = 0.28 # Point to ref shift line
y0y = -0.004 # Point to ref shift line

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
data0b = np.genfromtxt("../reference_data/curved_flaps/BottomLeafletObserver_Position_0000000000.dat")
data0t = np.genfromtxt("../reference_data/curved_flaps/TopLeafletObserver_Position_0000000000.dat")
#data1x = np.genfromtxt("../reference_data/curved_flaps/output_2_curvedFlaps_thin_2014_x.tsv")
data1y = np.genfromtxt("../reference_data/curved_flaps/output_2_curvedFlaps_thin_2014_y.tsv")
data2x = np.genfromtxt("../reference_data/curved_flaps/output_2_curvedFlaps_thin_2016_x.tsv")
data2y = np.genfromtxt("../reference_data/curved_flaps/output_2_curvedFlaps_thin_2016_y.tsv")
# Time
t0b = data0b[:,0]
t0t = data0t[:,0]
#t1x = data1x[:,0]
t1y = data1y[:,0]
t2x = data2x[:,1]
t2y = data2y[:,0]
# Y coordinate
x0b = data0b[:,1] - 0.0*data0b[1,1] + 1.0e-4
y0b = data0b[:,2] - 0.0*data0b[1,2] + flap_thickness/2.0
x0t = data0t[:,1] - 0.0*data0t[1,1] + 1.0e-4
y0t = data0t[:,2] - 0.0*data0t[1,2] - flap_thickness/2.0
#x1 = data1x[:,1]/100
y1b = data1y[:,1]/100
y1t = data1y[:,2]/100
x2 = data2x[:,0]/100
y2b = data2y[:,1]/100
y2t = data2y[:,2]/100

# Plotting the X displacement
#plt.plot(t1x, x1, color='C2', linestyle=':', marker='+', label='FSITICT Wick 2014')
plt.plot(t2x, x2, color='C4', linestyle='none', marker='*', label='arbitrary Lagrangian–Eulerian (ALE) Liu 2016')
plt.plot(t0b, x0b, color='C1', linewidth=SPH_LINE, label='SPH Linear Elastic Isotropic material')
plt.plot(t0t, x0t, color='C1', linewidth=SPH_LINE)
plt.plot(t0b + shift_time, x0b, linestyle=':', marker='', color='C1', linewidth=SPH_LINE)
plt.plot(t0t + shift_time, x0t, linestyle=':', marker='', color='C1', linewidth=SPH_LINE)
# Plot Vertical line
plt.vlines(t0x, -0.1, 0.1, color='gray', linestyle=':')
plt.vlines(t0x + shift_time, -0.1, 0.1, color='gray', linestyle=':')
# Plot Horizontal arrow
plt.annotate('', xy=(t0x,y0x), xytext=(t0x + shift_time,y0x), color='gray', arrowprops=dict(arrowstyle='<->'))
# Annotate 
plt.annotate('$\delta t = \dfrac{X_{D\'E\'}}{c_s}$',
            xy = (t0x*1.02 + shift_time, y0x*1.002), 
            xytext =(t0x*1.02 + shift_time, y0x*1.002), 
            color='black',
            fontweight='bold',
            fontsize = 10)
plt.xlim([tmin, tmax])
plt.ylim([xmin, xmax])
plt.xlabel('Time (s)')
plt.ylabel('X coordinate (m)')
legend = plt.legend(loc='best', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('C0') # Put a nicer background color on the legend
fig = plt.gcf() # get current figure
fig.savefig('x-disp-top-2-flaps-2mm-0-Linear-ElasticWall.png')
plt.show()

# Plotting the Y displacement
plt.plot(t1y, y1b, color='C2', linestyle=':', marker='+', label='FSITICT Wick 2014')
plt.plot(t1y, y1t, color='C2', linestyle=':', marker='+')
plt.plot(t2y, y2b, color='C4', linestyle=':', marker='*', label='arbitrary Lagrangian–Eulerian (ALE) Liu 2016')
plt.plot(t2y, y2t, color='C4', linestyle=':', marker='*')
plt.plot(t0b, y0b, color='C1', linewidth=SPH_LINE, label='SPH Linear Elastic Isotropic material')
plt.plot(t0t, y0t, color='C1', linewidth=SPH_LINE)
plt.plot(t0b + shift_time, y0b, linestyle=':', marker='', color='C1', linewidth=SPH_LINE)
plt.plot(t0t + shift_time, y0t, linestyle=':', marker='', color='C1', linewidth=SPH_LINE)
# Plot Vertical line
plt.vlines(t0y, -0.1, 0.1, color='gray', linestyle=':')
plt.vlines(t0y + shift_time, -0.1, 0.1, color='gray', linestyle=':')
# Plot Horizontal arrow
plt.annotate('', xy=(t0y,y0y), xytext=(t0y + shift_time,y0y), color='gray', arrowprops=dict(arrowstyle='<->'))
# Annotate 
plt.annotate('$\delta t = \dfrac{X_{D\'E\'}}{c_s}$',
            xy = (t0y*1.1 + shift_time,y0y*0.94), 
            xytext =(t0y*1.05 + shift_time,y0y*0.94), 
            color='black',
            fontweight='bold',
            fontsize = 10)
plt.xlim([tmin, tmax])
plt.ylim([ymin, ymax])
plt.xlabel('Time (s)')
plt.ylabel('Y coordinate (m)')
legend = plt.legend(loc='best', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('C0') # Put a nicer background color on the legend
fig = plt.gcf() # get current figure
fig.savefig('y-disp-top-2-flaps-2mm-0-Linear-ElasticWall.png')
plt.show()