
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as pl
get_ipython().magic('matplotlib inline')
from matplotlib import pylab
import FTS_Module as fts

from matplotlib import rc
rc('font', size =18)
rc('font', family = 'serif')
rc('font', style ='normal')
rc('font', serif ='DejaVu Serif')
rc('xtick', direction ='in')
rc('ytick', direction ='in')


#Note: Input source is at the focus of M1, and DH is at the focus of M2.

#Assumptions:
#We are evaluating the 1.1 mm wavelength light.
#The incoming beam waist to the FTS is located at the input.

#Goals:
#Determine best initial beam waist for coupling the FTS optics to TolTEC.
#Based on this beam waist, find the distance between our mirror and the input/the effective focal length of our mirror.
#Have three overall propagation plots for the different path lengths a beam could encounter:
#1. 7.5" (equal), 2. 3" (closest), 3. 12" (farthest)

Variable Key:

DH_bottom = bottom of dihedral
DH_top = top of dihedral

(front = element before dihedral, back = element after dihedral)

Outputs of Gaussian Propagation Function: (could be lens or free or combo as input prop matrix)
d(element)_front/back_out = distance between element and its beam waist
w(element)_front/back_out = size of beam waist after the propagation in question

Output of Free Propagation Function:
w(element)_front/back = size of the beam after propagating the given distance from the given beam waist
# In[2]:

#Current questions/ideas:
#1. Why does an asymmetry between the front and back occur for initial waists larger than 5mm?

#If we wanted to change the location of the initial beam waist (to not be at the input):
#we should change fts.distM1_in to be (140/cos_deg(4.8)) + distance from waist to input.


# In[3]:

#Determine the best beam waist for each wavelength:
w0_11 = fts.best_w0_finder(1.1)
w0_14 = fts.best_w0_finder(1.4)
w0_21 = fts.best_w0_finder(2.1)


# In[4]:

#We want these values to match for the best coupling, so we alter the input beam waist:

#Which wavelength we're interested in for this run (user input):
wavelength = eval(input("What is the wavelength of the beam of light going into the FTS? (mm)"))

#Initial Beam Waist [mm] (user input):
w0_var = eval(input("What is the initial beam waist at the input? (mm) ")) #Input the beam waist corresponding to the wavelength you want.


#Calculate Radius of Curvature for M1:
y, f, ROC_parab_M1, ROC_M1  = fts.newmirror(fts.f_effM1)
print('Radius of Curvature of M1:',ROC_M1,'mm')

#To determine ROC of the beam at M1, propagate initial beam waist to M1.
ROC_beam = fts.ROC_beam(wavelength,w0_var,fts.distM1_in)
print("Radius of Curvature of the Beam at M1:",ROC_beam,'mm')


# In[5]:

#First, we propagate from Input to Dihedral (one element at a time):

print('One by One Values (Front):')

#Propagate beam from input (assumed location of initial beam waist) to and through Mirror 1:
dM1_front_out, wM1_front_out = fts.Gauss_op(fts.lens(fts.f_effM1), wavelength,w0_var,fts.distM1_in)
print('1.1 Beam Waist after M1:',wM1_front_out,'mm')
print('Distance between M1 and its Output Beam Waist:',dM1_front_out,'mm')


print('')
#figure out beam size at M1
wM1_front = fts.freeprop(wM1_front_out, wavelength, dM1_front_out)
print(wavelength,'Beam Size at Mirror 1:',wM1_front,'mm')


print('')
#Propagate beam from M1's beam waist to and through Mirror 2:
dM2_front_out, wM2_front_out = fts.Gauss_op(fts.lens(fts.f_effM2), wavelength, wM1_front_out,fts.distM1_M2 - dM1_front_out)
print(wavelength, 'Beam Waist after M2:',wM2_front_out,'mm')
print('Distance between M2 and its Output Beam Waist:',dM2_front_out,'mm')


print('')
#figure out beam size at M2
wM2_front = fts.freeprop(wM2_front_out, wavelength, dM2_front_out)
print(wavelength, 'Beam Size at Mirror 2:',wM2_front,'mm')


print('')
#figure out beam size at bottom of Dihedral:
wDH_bottom = fts.freeprop(wM2_front_out, wavelength, fts.distM2_DH - dM2_front_out)
print(wavelength, 'Beam Size at Bottom of Dihedral:',wDH_bottom,'mm')


# In[6]:

#Create plot of beam size evolution from Input to bottom of Dihedral:

#horiz_ax_1 = [input, M1, M2, bottom of DH]
horiz_ax_1 = np.array([0,fts.distM1_in, fts.distM1_in + fts.distM1_M2, fts.distM1_in + fts.distM1_M2 + fts.distM2_DH]) #distance from input
vert_ax_1 = np.array([w0_var, wM1_front, wM2_front, wDH_bottom]) #beam size at each element

#Make the plot:
fig, ax = pl.subplots(1, figsize=(10,8))
pl.plot(horiz_ax_1, vert_ax_1,'b')
pl.plot(horiz_ax_1, -vert_ax_1,'b')


#Add element locations as vertical lines:
pl.vlines(horiz_ax_1[0],-5,5,'k',label = 'Input',linewidth = 3)
pl.vlines(horiz_ax_1[1],-fts.diam_M1/2,fts.diam_M1/2,'orchid',label = 'Mirror 1',linewidth = 3)
pl.vlines(horiz_ax_1[2],-fts.diam_M2/2,fts.diam_M2/2,'purple',label = 'Mirror 2',linewidth = 3)
pl.vlines(horiz_ax_1[3],-fts.length_DH/2,fts.length_DH/2,'maroon',label = 'Dihedral Bottom',linewidth = 3)

#Efective Beam Waist Locations also as vertical lines:
pl.vlines(fts.distM1_in + dM1_front_out,-wM1_front_out,wM1_front_out,'darkgreen',label = 'Mirror 1 BW',linewidth = 4)
pl.vlines(fts.distM1_in + fts.distM1_M2 + dM2_front_out,-wM2_front_out,wM2_front_out,'lightgreen',label = 'Mirror 2 BW',linewidth = 4)


#Format plot:
pl.xlabel('Distance from Input (mm)')
pl.ylabel('Beam Size (mm)')

title_string_1 = "FTS Beam Size Evolution From Input to Bottom of Dihedral"
subtitle_string_1 = 'Wavelength = {0} mm, Dihedral Offset from Center = 7.5"'

pl.suptitle(title_string_1.format(wavelength), y=0.96, fontsize=18)
pl.title(subtitle_string_1.format(wavelength), fontsize=17)

pl.legend()
pylab.legend(loc=9, bbox_to_anchor=(1.25, 0.65))

pl.show()


# In[7]:

#The last step is to propagate from bottom of Dihedral to Output (again one by one):

print('One by One Values Back to Output:')

#figure out beam size at top of Dihedral:
wDH_top = fts.freeprop(wM2_front_out,1.1, (fts.distM2_DH - dM2_front_out + fts.inch2mm(3)))
print('1.1 Beam Size at Top of Dihedral:',wDH_top,'mm')


print('')
#Propagate beam from (front) M2's beam waist to and through M2 (back):
dM2_back_out, wM2_back_out = fts.Gauss_op(fts.lens(fts.f_effM2),1.1,wM2_front_out,(fts.distM2_DH - dM2_front_out) + fts.inch2mm(3) + fts.distM1_M2)
print('1.1 Beam Waist after M2 (back):',wM2_back_out,'mm')
print('Distance between M2 (back) and its Beam Waist:',dM2_back_out,'mm')


print('')
#figure out beam size at M2 (back):
wM2_back = fts.freeprop(wM2_back_out,1.1,dM2_back_out)
print('1.1 Beam Size at Mirror 2 (back):', wM2_back,'mm')


print('')
#Propagate from (back) M2's beam waist to and through M1:
dM1_back_out, wM1_back_out = fts.Gauss_op(fts.lens(fts.f_effM1),1.1,wM2_back_out,dM2_back_out - fts.distM1_M2)
print('1.1 Beam Waist after M1 (back):',wM1_back_out,'mm')
print('Distance between M1 (back) and its beam waist:',dM1_back_out,'mm')


print('')
#figure out beam size on M1 (back):
wM1_back = fts.freeprop(wM1_back_out,1.1, dM1_back_out)
print('1.1 Beam Size at Mirror 1 (back):',wM1_back,'mm')


print('')
#figure out beam size at output
#propagate from M1's beam waist back to output (back propagation)
w_output = fts.freeprop(wM1_back_out,1.1, dM1_back_out - fts.distM1_out)
print('1.1 Beam Size at Output:',w_output,'mm')


# In[8]:

#Create plot of beam size evolution from bottom of dihedral to output:

#horiz_ax_2 = [top of DH, M2, M1, output]
horiz_ax_2 = np.array([fts.distM1_in + fts.distM1_M2 + fts.distM2_DH + fts.distDHb_DHt,fts.distM1_in + fts.distM1_M2 + fts.distM2_DH + fts.distDHb_DHt + fts.distM2_DH, fts.distM1_in + fts.distM1_M2 + fts.distM2_DH + fts.distDHb_DHt + fts.distM2_DH + fts.distM1_M2,fts.distM1_in + fts.distM1_M2 + fts.distM2_DH + fts.distDHb_DHt + fts.distM2_DH + fts.distM1_M2 + fts.distM1_out]) #distance from input
vert_ax_2 = np.array([wDH_top, wM2_back, wM1_back, w_output])


#Make the plot:
fig, ax = pl.subplots(1, figsize=(10,8))
pl.plot(horiz_ax_2, vert_ax_2,'b')
pl.plot(horiz_ax_2, -vert_ax_2,'b')


#Add element locations as vertical lines:
pl.vlines(horiz_ax_2[0],-fts.length_DH/2,fts.length_DH/2,'maroon',label = 'Dihedral Top',linewidth = 3)
pl.vlines(horiz_ax_2[1],-fts.diam_M2/2,fts.diam_M2/2,'purple',label = 'Mirror 2',linewidth = 3)
pl.vlines(horiz_ax_2[2],-fts.diam_M1/2,fts.diam_M1/2,'orchid',label = 'Mirror 1',linewidth = 3)
pl.vlines(horiz_ax_2[3],-5,5,'k',label = 'Output')


#Need to add effective beam waist locations.

#Format plot:
pl.xlabel('Distance from Input (mm)')
pl.ylabel('Beam Size (mm)')

title_string_3 = "FTS Beam Size Evolution From Top of Dihedral to Output"
subtitle_string_3 = 'Wavelength = {0} mm, Dihedral Offset from Center = 7.5"'

pl.suptitle(title_string_3.format(wavelength), y=0.96, fontsize=18)
pl.title(subtitle_string_3.format(wavelength), fontsize=17)

pl.legend()
pylab.legend(loc=9, bbox_to_anchor=(1.25, 0.65))

pl.show()


# In[9]:

#Create overall evolution plot by combining the two previous plots (using concatenate command):

horiz_tot = np.concatenate((horiz_ax_1,horiz_ax_2))
vert_tot = np.concatenate((vert_ax_1,vert_ax_2))


#Make the plot:
fig, ax = pl.subplots(1, figsize=(10,8))
pl.plot(horiz_tot, vert_tot,'b')
pl.plot(horiz_tot, -vert_tot,'b')


#Add element locations as vertical lines:
pl.vlines(horiz_tot[0],-5,5,'k',label = 'Input')
pl.vlines(horiz_tot[1],-fts.diam_M1/2,fts.diam_M1/2,'orchid',label = 'Mirror 1',linewidth = 3)
pl.vlines(horiz_tot[2],-fts.diam_M2/2,fts.diam_M2/2,'purple',label = 'Mirror 2',linewidth = 3)
pl.vlines(horiz_tot[3],-fts.length_DH/2,fts.length_DH/2,'maroon',label = 'Dihedral',linewidth = 3)
pl.vlines(horiz_tot[4],-fts.length_DH/2,fts.length_DH/2,'maroon',linewidth = 3)
pl.vlines(horiz_tot[5],-fts.diam_M2/2,fts.diam_M2/2,'purple',linewidth = 3)
pl.vlines(horiz_tot[6],-fts.diam_M1/2,fts.diam_M1/2,'orchid',linewidth = 3)
pl.vlines(horiz_tot[7],-5,5,'k',label = 'Output')

#Format plot:
pl.xlabel('Distance from Input (mm)')
pl.ylabel('Beam Size (mm)')

title_string_3 = "FTS Beam Size Evolution From Input to Output"
subtitle_string_3 = 'Wavelength = {0} mm, Dihedral Offset from Center = 7.5"'

pl.suptitle(title_string_3.format(wavelength), y=0.96, fontsize=18)
pl.title(subtitle_string_3.format(wavelength), fontsize=17)

pl.legend()
pylab.legend(loc=9, bbox_to_anchor=(1.25, 0.65))

pl.show()


# In[ ]:



