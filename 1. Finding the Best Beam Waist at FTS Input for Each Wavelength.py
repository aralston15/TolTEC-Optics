
# coding: utf-8

# In[2]:

import numpy as np
import matplotlib.pyplot as pl
get_ipython().magic('matplotlib inline')
from matplotlib import pylab
import FTS_Module as fts
import sympy as sp

from matplotlib import rc
rc('font', size =18)
rc('font', family = 'serif')
rc('font', style ='normal')
rc('font', serif ='DejaVu Serif')
rc('xtick', direction ='in')
rc('ytick', direction ='in')


#Note: Input source is at the focus of M1, and DH is at the focus of M2.

#Assumptions:
#The incoming beam waist to the FTS is located at the input.

#Goal:
#Determine best initial beam waist for coupling the FTS optics to TolTEC (for 1.1, 1.4 and 2.1).


# In[3]:

#Finding the best beam waist for each wavelength:

def best_w0_finder(wavelength):

    w0_guess = np.linspace(0.1,25,100000) #an array of random guesses for the beam waist at the input (in mm)
    #might not work if the # of guesses is not large enough, or if the wavelength doesn't fall in this range

    #Calculate Radius of Curvature for M1:
    y, f, ROC_M1  = fts.newmirror(fts.f_effM1) #the new mirror function works because we use its effective focal length


    ROC_beam_array = np.zeros(len(w0_guess)) #empty array to hold the new beam ROC values (to match with M1's ROC)
    for i in range(len(ROC_beam_array)): #evaluate ROC of the beam at M1 for all the different beam waist guesses

        ROC_beam_test = fts.ROC_beam(wavelength,w0_guess[i],fts.f_effM1) #gives ROC of the beam at that f_eff (propagates forward f_eff bc beam waist is at the focus)
        ROC_beam_array[i] = ROC_beam_test #fills empty ROC_beam_array
        
    #Now we can find where these arrays match in value:

    index_best_tuple = np.where(np.abs(ROC_beam_array - ROC_M1) <= 0.01) #use np.where to find the index of the best ROC's (and w0)
    index_best = int(index_best_tuple[0]) #convert to integer


    #Print final values:
    print('For', wavelength, 'mm wavelength:')
    print('')
    print('***Best Beam Waist:',w0_guess[index_best],'mm')
    print('Gives:')
    print('ROC of Mirror 1:',ROC_M1,'mm')
    print('ROC of Beam:',ROC_beam_array[index_best],'mm')

    print('')
    print('Difference in ROC:',np.abs(ROC_M1 - ROC_beam_array[index_best]),'mm')
    print('Percent Difference:',fts.percent_diff(ROC_M1,ROC_beam_array[index_best]),'%')
    print('')
    print('')


# In[4]:

#Finally, run the function for the three wavelengths we're interested in:
best_w0_finder(1.1)
best_w0_finder(1.4)
best_w0_finder(2.1)


# In[11]:

#Plots various parabolas at different angles (for the sake of better visualizing the new mirror.)

t = np.linspace(-5,5,50)

angle_array = np.linspace(0,360,4) #in degrees

x_array = np.zeros(shape=(len(angle_array),len(t)))
y_array = np.zeros(shape=(len(angle_array),len(t)))

for i in range(len(angle_array)):
    
    x = (t*fts.cos_deg(angle_array[i])) - ((t**2)*fts.sin_deg(angle_array[i]))
    x_array[i] = x
    
    y = ((t*fts.sin_deg(angle_array[i])) + ((t**2)*fts.cos_deg(angle_array[i])))
    y_array[i] = y

for i in range(len(angle_array)):
    pl.plot(x_array[i],y_array[i])

pl.show()


# In[12]:

pl.plot(x_array, y_array)
pl.show()


# In[ ]:



