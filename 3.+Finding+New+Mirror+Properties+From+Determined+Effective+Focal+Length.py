
# coding: utf-8

# In[1]:

import numpy as np
import sympy as sp
import matplotlib.pyplot as pl
get_ipython().magic('matplotlib inline')
from matplotlib import pylab
import FTS_Module as fts
x = sp.Symbol('x') #define general variable for equation

from matplotlib import rc
rc('font', size =18)
rc('font', family = 'serif')
rc('font', style ='normal')
rc('font', serif ='DejaVu Serif')
rc('xtick', direction ='in')
rc('ytick', direction ='in')

#Here we can try out new effective focal lengths/distances from 800m and up for the properties of the new mirror.


# In[2]:

f_eff = float(input('What is the effective focal length (dist between FTS and Mirror)? (mm) '))
y_parab, f, ROC_parab, ROC_mirror = fts.newmirror(f_eff) #Evaluate mirror parameters with function in module for eq of mirror, f, and ROC

A_co = y_parab.coeff(x**2) #extract A coefficent from equation for the parabola
A_co_py = float(A_co) #convert to float so python can read

#Print these values:
print('')
print('Equation for Parabola:',y_parab,'mm')
print('Primary Focal Length:',f,'mm')
print('Radius of Curvature of Parent Parabola:',ROC_parab,'mm')
print('Radius of Curvature of Mirror:',ROC_mirror,'mm')


#So, given any f_eff (distance beam travels from input to new mirror),
#we can calculate the equation of the mirror, the primary focal length, and the ROC.
#Condition: angle of reflection must always 90 degrees.


#Find corresponding y-coordinate for the effective focal length chosen:
y_f_eff = A_co_py*(f_eff)**2
print('Effective Focal Point Location:','(',f_eff,y_f_eff,')')

print('')
#Find equation of the new mirror:
mirror_eq = fts.mirrorline(f_eff)
print('Equation for Mirror Line:',mirror_eq)

mirror_slope = mirror_eq.coeff(x) #extract slope from the new mirror line (for plotting)
mirror_int = mirror_eq.coeff(x,0)#extract y-intercept from the new mirror line (again for plotting)


# In[3]:

#Define the x and y coordinates of the parent parabola:
horiz_ax = np.linspace(-(f_eff + 300),f_eff + 300, 1000)
vert_ax = A_co_py * (horiz_ax**2)


#Plot the parabola:
fig, ax = pl.subplots(1, figsize=(10,8))

pl.plot(horiz_ax, vert_ax) #plot the parabola

#Plot relevant points:
pl.plot(0,f,'bo',linewidth = 100,label = 'Focal Point/FTS Input Location') #plot the focus

pl.plot(f_eff,y_f_eff,'ro',label = 'Point where Beam is Incident') #plot the point where beam is incident on mirror

pl.plot(f_eff,y_f_eff + 300,'go',label = 'Minimum TolTEC Clearance Distance') #plot the point where beam is free of TolTEC


#Plot the beams of light:
pl.plot(f_eff*np.ones(100),np.linspace(y_f_eff,y_f_eff + 300,100),'k--',label = 'Beam Path') #path of light from TolTEC window to mirror
pl.plot(np.linspace(0,f_eff,100),f*np.ones(100),'k--') #path of light from mirror to FTS


#Plot the mirror line:
pl.plot(horiz_ax[500:], (mirror_slope*horiz_ax[500:]) + mirror_int, 'm',linewidth = 3,label = 'New Mirror')


#Format Plot
pl.title('New FTS Mirror with Effective Focal Length = %s'%(f_eff))
pl.legend()       
pylab.legend(loc=9, bbox_to_anchor=(1.5, 0.5))
pl.xlabel('x (mm)')
pl.ylabel('y (mm)')
#pl.xlim(-f_eff - 200, 10)
pl.show()


# In[ ]:



