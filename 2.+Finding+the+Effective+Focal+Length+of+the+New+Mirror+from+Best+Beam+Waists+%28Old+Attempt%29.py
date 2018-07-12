
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


# In[2]:

#Effective Focal Length = distance from FTS Input to New Mirror

#Beam waist at input for 1.1 mm (with equal path lengths):
w0_11 = fts.best_w0_finder(1.1) #mm (updated 6/28/18 for new M1 ROC)

#Beam waist at input for 1.4 mm (with equal path lengths):
w0_14 = fts.best_w0_finder(1.4) #mm (updated 6/28/18 for new M1 ROC)

#Beam waist at input for 2.1 mm (with equal path lengths):
w0_21 = fts.best_w0_finder(2.1) #mm (updated 6/28/18 for new M1 ROC)


# In[3]:

#Steps:
#Choose arbitrary f_eff. This determines the ROC of the mirror.
#Then, propagate 8.14mm beam from input f_eff that gives an ROC of the beam matching that of the mirror.


# In[4]:

#Want to plot the evolution of the beam's radius of curvature from the FTS, outward 
#Note: can do for any wavelength, I use 1.1

pl.ion()

#Define horizontal and vertical axes.
horiz_ROC_FTS = np.linspace(10,1000,1000)
vert_ROC_FTS = np.zeros(len(horiz_ROC_FTS)) #Empty because we fill it with the for loop below.


for i in range(len(horiz_ROC_FTS)): #Evaluates the ROC of the beam at various f_eff's.
    ROC_beam = fts.ROC_beam(1.1,w0_11,horiz_ROC_FTS[i])
    vert_ROC_FTS[i] = ROC_beam
    
    
#Plotting the beam's radius of curvature as a function of the effective focal length:
fig, ax = pl.subplots(1, figsize=(10,8))
pl.plot(horiz_ROC_FTS,vert_ROC_FTS,'k')


#Format plot:
pl.xlabel('Effective Focal Length/ Distance from FTS (mm)')
pl.ylabel('Radius of Curvature (mm)')
pl.title('ROC of the Beam from the FTS for 1.1 mm')

pl.vlines(300,0,3500,'g',label = 'Minimum Allowed Effective Focal Length for TolTEC Clearance') #300 mm is the min f_eff
pl.vlines(139.948423484,0,3500,'r',label = 'Calculated Best Effective Focal Length') #the best f_eff

pl.legend()
pylab.legend(loc=9, bbox_to_anchor=(1.6, 0.65))
pl.show()

#This behavior makes sense because at the beam waist (f_eff = 0), the ROC is infinite (plane wave).


# In[5]:

#Also want to plot the beam size as it evolves from the input of the FTS outward:

#Create beam size and corresponding effective focal length arrays:
horiz_size_FTS = np.linspace(10,1000,1000)
vert_size_FTS = np.zeros(len(horiz_size_FTS))

for i in range(len(horiz_size_FTS)):
    beamsize = fts.freeprop(w0_11, 1.1, horiz_size_FTS[i]) #Evaluates the size of the beam at various f_eff's.
    vert_size_FTS[i] = beamsize #fills the vert_size array with beam sizes
    
    
#Plotting the beam size as a function of the effective focal length:
fig, ax = pl.subplots(1, figsize=(10,8))
pl.plot(horiz_size_FTS,vert_size_FTS,'b')
pl.plot(horiz_size_FTS,-vert_size_FTS,'b') #Plot negative y-vals to get image of entire beam.


#Format plot:
pl.xlabel('Effective Focal Length (mm)')
pl.ylabel('Beam Size (mm)')
pl.title('Beam Size of the Beam from the FTS')


#Plotting Focal Length Limits:
pl.vlines(300,-50,50,'g',label = 'Minimum Allowed Effective Focal Length for TolTEC Clearance') #300 mm is the min f_eff
pl.vlines(139.948423484,-50,50,'r',label = 'Calculated Best Effective Focal Length') #the best f_eff

pl.legend()
pylab.legend(loc=9, bbox_to_anchor=(1.6, 0.65))
pl.show()


# In[6]:

print(vert_size_FTS)


# In[ ]:




# In[10]:


#Want to plot the beam size as it evolves from the Lyot Stop outward (to see if it matches with FTS propagation):

horiz_size_Lyot = np.linspace(10,np.abs(fts.distw0_M7) + fts.dist2M7 + 200,1000)
vert_size_Lyot = np.zeros(len(horiz_size_Lyot))

for i in range(len(horiz_size_Lyot)): #Evaluates the size of the beam at various f_eff's.
    beamsize = fts.freeprop(w0_11,1.1, horiz_size_Lyot[i])
    vert_size_Lyot[i] = beamsize
    
    
    
#Plotting the beam size as a function of the effective focal length:
fig, ax = pl.subplots(1, figsize=(10,8))
pl.plot(horiz_size_Lyot,vert_size_Lyot,'b')
pl.plot(horiz_size_Lyot,-vert_size_Lyot,'b') #Plot negative y-vals to get image of entire beam.


#Format plot:
pl.xlabel('Distance from the Beam Waist (mm)')
pl.ylabel('Beam Size (mm)')
pl.title('Beam Size of the Beam from the Lyot Stop')


#Plotting Focal Length Limits:
pl.vlines(0,-50,50,'m',label = 'Beam Waist due to M7')
pl.vlines(np.abs(fts.distw0_M7),-fts.diam_M7/2,fts.diam_M7/2,'m',label = 'M7')
pl.vlines(np.abs(fts.distw0_M7) + fts.dist2M7,-fts.LS_diam/2,fts.LS_diam/2,'m',label = 'Lyot Stop')

#pl.legend()
#pylab.legend(loc=9, bbox_to_anchor=(1.45, 0.65))
pl.show()


# # Finding the Best Effective Focal Length (based on matching ROC) at the New Mirror

# In[9]:

#This cell will return the Radii of Curvature for both the Mirror and the Beam based on any input f_eff.

#(Trial and Error Method)

f_eff = float(input('What is the effective focal length of the new mirror? (mm) ')) #Choose effective focal length.
y, f, ROC_mirror = fts.newmirror(f_eff)
print('Radius of Curvature of this Mirror:',ROC_mirror,'mm')

print('')
#Determine the ROC of the beam after propagating the chosen f_effective:
ROC_beam = fts.ROC_beam(1.1,w0,f_eff) #(wavelength, beam waist, distance to propagate)
print('Radius of Curvature of the Beam after Propagating %s'%(f_eff),'mm :',ROC_beam,'mm')

print('')
print('Difference in ROC:',np.abs(ROC_mirror - ROC_beam),'mm')

def percent_diff(first, second): #returns percent difference
    return 100*np.abs((first - second)/((first + second)/2))

print('Percent Difference:',percent_diff(ROC_mirror,ROC_beam),'%')

#The smaller the effective focal length, the more the ROC's match.


# # A More Exact Way to Find the Best f_eff:

# In[28]:

def best_f_eff_finder(wavelength, w0): #w0 is the beam waist that corresponds with the wavelength you choose
    
    #A More Exact Way to Find the Best f_eff:
    #For Loop Method, with Np.Where to Determine the Best f_eff for Matching Beam and Mirror ROC
    
    f_eff_guess = np.linspace(0.1,500,100000) #an array of random guesses for f_eff (in mm)
    
    ROC_mirror_array = np.zeros(len(f_eff_guess)) #empty array to fill with possible mirror ROCS
    ROC_beam_array = np.zeros(len(f_eff_guess)) #empty array to fill with possible beam ROCs
    
    
    for i in range(len(f_eff_guess)): #evaluates the ROC of the beam and mirror at the various f_eff values
    
        y_vals, f_vals, ROC_parab_vals, ROC_mirror_vals  = fts.newmirror(f_eff_guess[i]) #gives ROC of mirror at that f_eff
        ROC_mirror_array[i] = ROC_mirror_vals #fills empty ROC_mirror_array
    
        ROC_beam_vals = fts.ROC_beam(wavelength,w0,f_eff_guess[i]) #gives ROC of the beam at that f_eff
        ROC_beam_array[i] = ROC_beam_vals #fills empty ROC_beam_array
        
        
        
    #Now we can find where these arrays match in value:
    
    #use np.where to find the index of the best f_eff (and its corresponding mirror + beam ROCs)
    index_best_tuple = np.where(np.abs(ROC_mirror_array - ROC_beam_array) <= 0.01)
    index_best_array = index_best_tuple[0]#honestly wut but okay (selects the first array in the list of arrays)
    index_best = int(index_best_array[0]) #convert to integer


    #Print final values:
    print('Best Effective Focal Length:',f_eff_guess[index_best],'mm')
    print('Gives:')
    print('ROC of Mirror:',ROC_mirror_array[index_best],'mm')
    print('ROC of Beam:',ROC_beam_array[index_best],'mm')

    print('')
    print('Difference in ROC:',np.abs(ROC_mirror_array[index_best] - ROC_beam_array[index_best]),'mm')
    print('Percent Difference:',fts.percent_diff(ROC_mirror_array[index_best],ROC_beam_array[index_best]),'%')

    y_best, f_best, ROC_parab_best, ROC_mirror_best = fts.newmirror(f_eff_guess[index_best])
    print('')
    print('New Mirror Properties:')
    print('y = ',y_best)
    print('Primary Focal Length:',f_best)
    print('Radius of Curvature of Parent Parabola:', ROC_parab_best)
    print('Radius of Curvature of New Mirror:', ROC_mirror_best)
    
    return f_eff_guess[index_best] #returns the best effective focal length


# In[30]:

best_f_eff_finder(1.4, w0_14)


# In[31]:

best_f_eff_finder(2.1, w0_21)


# In[15]:

#A More Exact Way to Find the Best f_eff:
#For Loop Method, with Np.Where to Determine the Best f_eff for Matching Beam and Mirror ROC
    
f_eff_guess = np.linspace(0.1,500,100000) #an array of random guesses for f_eff (in mm)
    
ROC_mirror_array = np.zeros(len(f_eff_guess)) #empty array to fill with possible mirror ROCS
ROC_beam_array = np.zeros(len(f_eff_guess)) #empty array to fill with possible beam ROCs
    
    
for i in range(len(f_eff_guess)): 
    
    y_vals, f_vals, ROC_parab_vals, ROC_mirror_vals  = fts.newmirror(f_eff_guess[i]) #gives ROC of mirror at that f_eff
    ROC_mirror_array[i] = ROC_mirror_vals #fills empty ROC_mirror_array
    
    ROC_beam_vals = fts.ROC_beam(1.1,w0_11,f_eff_guess[i]) #gives ROC of the beam at that f_eff
    ROC_beam_array[i] = ROC_beam_vals #fills empty ROC_beam_array


# In[17]:

#Now we can find where these arrays match in value:
    
#use np.where to find the index of the best f_eff (and its corresponding mirror + beam ROCs)
index_best_tuple = np.where(np.abs(ROC_mirror_array - ROC_beam_array) <= 0.01) 
index_best = int(index_best_tuple[0]) #convert to integer


#Print final values:
print('Best Effective Focal Length:',f_eff_guess[index_best],'mm')
print('Gives:')
print('ROC of Mirror:',ROC_mirror_array[index_best],'mm')
print('ROC of Beam:',ROC_beam_array[index_best],'mm')

print('')
print('Difference in ROC:',np.abs(ROC_mirror_array[index_best] - ROC_beam_array[index_best]),'mm')
print('Percent Difference:',percent_diff(ROC_mirror_array[index_best],ROC_beam_array[index_best]),'%')

y_best, f_best, ROC_M_best = fts.newmirror(f_eff_guess[index_best])
print('')
print('New Mirror Properties:')
print('y = ',y_best)
print('Primary Focal Length:',f_best)
print('Radius of Curvature:',ROC_M_best)


# In[23]:

index_best_tuple


# In[ ]:




# In[13]:

#A More Exact Way to Find the Best f_eff:
#For Loop Method, with Np.Where to Determine the Best f_eff for Matching Beam and Mirror ROC

f_eff_guess = np.linspace(0.1,500,100000) #an array of random guesses for f_eff (in mm)

ROC_mirror_array = np.zeros(len(f_eff_guess))
ROC_beam_array = np.zeros(len(f_eff_guess))

for i in range(len(f_eff_guess)): 
    
    y_test, f_test, ROCROC_mirror_test  = fts.newmirror(f_eff_guess[i]) #gives ROC of mirror at that f_eff
    ROC_mirror_array[i] = ROC_mirror_test #fills empty ROC_mirror_array
    
    ROC_beam_test = fts.ROC_beam(2.1,w0,f_eff_guess[i]) #gives ROC of the beam at that f_eff
    ROC_beam_array[i] = ROC_beam_test #fills empty ROC_beam_array
    


# In[24]:

#Now we can find where these arrays match in value:

index_best_tuple = np.where(np.abs(ROC_mirror_array - ROC_beam_array) <= 0.01) #use np.where to find the index of the best ROC's (and f_eff)
index_best = int(index_best_tuple[0]) #convert to integer


#Print final values:
print('Best Effective Focal Length:',f_eff_guess[index_best],'mm')
print('Gives:')
print('ROC of Mirror:',ROC_mirror_array[index_best],'mm')
print('ROC of Beam:',ROC_beam_array[index_best],'mm')

print('')
print('Difference in ROC:',np.abs(ROC_mirror_array[index_best] - ROC_beam_array[index_best]),'mm')
print('Percent Difference:',percent_diff(ROC_mirror_array[index_best],ROC_beam_array[index_best]),'%')

y_best, f_best, ROC_M_best = fts.newmirror(f_eff_guess[index_best])
print('')
print('New Mirror Properties:')
print('y = ',y_best)
print('Primary Focal Length:',f_best)
print('Radius of Curvature:',ROC_M_best)


# In[ ]:

#WRONG (use new function in module)

#Results:
#--For 1.1 mm light--#: 
#Best Effective Focal Length: 139.948423484 mm
#ROC of Both the Beam and New Mirror: 395.83 mm
#Primary Focal Length of New Mirror: 69.9742117421 mm

#--For 1.4 mm light--#:
#Best Effective Focal Length: 109.959122591 mm
#ROC of Both the Beam and New Mirror: 311.011 mm
#Primary Focal Length: 54.9795612956 mm
 
#--For 2.1 mm light--#:
#Best Effective Focal Length: 73.3060880609 mm
#ROC of Both the Beam and New Mirror: 207.34 mm
#Primary Focal Length: 36.6530440304 mm

