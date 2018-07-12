# coding: utf-8

import numpy as np
import sympy as sp
x = sp.Symbol('x') #define general variable for equation

#useful global variables and functions

pi = np.pi
def sqrt(a):
    return np.sqrt(a)
def exp(a):
    return np.exp(a)


#----Conversion Functions----#

def sin_deg(degrees): #evaluates sine given an input of degrees
    rad = degrees*(pi/180)
    return np.sin(rad)

def cos_deg(degrees): #evaluates cosine given an input of degrees
    rad = degrees*(pi/180)
    return np.cos(rad)

def tan_deg(degrees): #evaluates tangent given an input of degrees
    rad = degrees*(pi/180)
    return np.tan(rad)

def arctan_deg(opposite, adjacent): #evaluates inverse tangent given an input of degrees
    rad = np.arctan(opposite/adjacent)
    deg = radians*(180/np.pi)
    return rad, deg

def inch2mm(inches): #converts inches to mm
    return inches*25.4

def mm2inch(mm): #converts mm to inches
    return mm/25.4

def rad2deg(radians): #converts radians to degrees
    return radians*(180/pi)

def deg2rad(degrees): #converts degrees to radians
    return degrees*(pi/180)

    

#----Triangle Functions----#

#Side a = horizontal leg, side b = vertical leg, side c = hypotenuse
def pythag_a(b, c): #use when you want the length of a leg of a right triangle
    return np.sqrt((c**2) - (b**2))

def pythag_c(a, b): #use when you want the length of the hypotenuse of a right triangle
    return np.sqrt((a**2) + (b**2))


def lawofcos(): #a, b, and c are the sides (c must be the side opposite the angle you need)
    
    c = float(input('What is length of side C? ~opposite the angle you need~ (mm) ')) #asks length of side C (MUST BE OPPOSITE THE ANGLE)
    
    a = float(input('What is length of side A? (mm) ')) #asks length of side A (A and B are interchangable)
    
    b = float(input('What is length of side B? (mm) ')) #asks length of side B
    
    fraction = ((a**2) + (b**2) - (c**2))/(2*a*b) #computes the Law of Cosines
    rad = np.arccos(fraction) #gives angle in radians
    deg = rad*(180/np.pi) #gives angle in degrees
    return rad, deg #returns the angle opposite side C, in both radians and degrees


def lawofcos_alt(a,b,c): #Alternate version of Law of Cosine function where you give it a,b,c all at once
    #lawofcos() that takes no arguments is easier by itself for a human, whereas lawofcos_alt is better within functions
    
    fraction = ((a**2) + (b**2) - (c**2))/(2*a*b) #computes the Law of Cosines
    rad = np.arccos(fraction) #gives angle in radians
    deg = rad*(180/np.pi) #gives angle in degrees
    return rad, deg #returns the angle opposite side C, in radians and degrees


def lawofcos_SAS(angle, b, c): #function determines side A, then uses that to solve for angle B
    a_squared = b**2 + c**2 -2*b*c*fts.cos_deg(angle)
    a = np.sqrt(a_squared) #Law of Cosines to find side A
    sinB = b*fts.sin_deg(angle)/a #Law of Sines to find angle B
    degrees = np.arcsin(sinB) * (180/np.pi) #convert angle B to degrees
    return degrees #returns angle B in degrees



#----Optical Property Functions----#

#This function works for any off-axis parabolic mirror (with 90 degree reflection), in any unit:    
def newmirror(f_eff): #Given any f_eff, the function will return the ROC of the mirror and primary focal point.
    f = f_eff/2 #focal point of overall mirror (only true because 90 degree reflection)
   
    x = sp.Symbol('x') #define general variable for equation
    y = (1/(4*f))*(x**2) #general equation of the parabola
    
    ROC_mirror = np.sqrt(2**5)*f #radius of curvature **only works at the incident point**
    ROC_parabola = 2*f #radius of curvature **of the parent parabola**
    
    return y, f, ROC_parabola, ROC_mirror


def ROC_beam(wavelength, w0, z): #given a prop distance from the beam waist, returns the radius of curvature of the Gaussian beam
    z_r = (pi * (w0**2))/wavelength
    return z * (1 + (z_r/z)**2)


def focal(n, R1, R2): #provides the focal length of a lens (thin or thick) given the ROC
    return 1/((n -1)*(1/R1 + 1/R2)) #if thin use R1 (or R2) = np.infinity


def eff_focal(f1,f2,d): #returns the effective focal length of a combination lens system
    F_inv = (1/f1) + (1/f2) - (d/(f1*f2)) #d is the distance between the two lenses
    return 1/F_inv


#----General Functions----#

def free(d): #propagation matrix through free space
    return np.matrix([[1.,d],[0.,1.]])

def lens(f): #propagation matrix through a thin lens
    return np.matrix([[1.,0.],[-1./f,1.]])

def edgetaper(diameter,w): #returns edge taper (relative power density at a certain radius of the beam)
    r = diameter/2. #radius of the element
    Te = exp(-2.*(r/w)**2.) #edge taper
    Te_db = 8.686*(r/w)**2 #log value of edge taper
    
    return Te, Te_db

def FWHM(Te_db, wavelength, diameter): #returns full-width-half-maximum in arcseconds
    return (1.02 + 0.0135*Te_db) * wavelength/(diameter) * 57.3 * 3600.
    
#the lovely matrix multiplication equation created by the wonderful Natalie
def matMult(a): #input array of matricies in reverse order you want them multiplied (lin alg) [M3,M2,M1]
    a_f = a[-1] #final matrix in array
    for i in range(len(a)): #depending on length of array (# of matricies)
        if i !=0: #if i is not equal to 0, multiply the matrix immediately before by a_f
            a_f = np.dot(a[-i-1],a_f)
        else:
            continue
    return a_f


def deriv(y): #Evaluates derivatives symbolically (variable of differentiation must be 'x')
    x = sp.Symbol('x') #defines variable of differentiation
    return sp.diff(y) #y is the function you want to differentiate

#First attempt:
#def integrate(function,a_lim,b_lim): #Evaluates definite integrals numerically. (variable must again be x)
    #Need to give: function = lambda x: (function_in)
    
    #import scipy.integrate as integrate
    #result = integrate.quad(function,a_lim,b_lim)[0] #the integral of the function from a_lim to b_lim
    #abs_error = integrate.quad(function,a_lim,b_lim)[1] #estimate of absolute error in the result
    
    #return result #returns only the result

def integral(y): #Evaluates integrals symbolically (variable of integration must be 'x')
    x = sp.Symbol('x') #defines variable of integration
    return sp.integrate(y) #y is the function you want to integrate

def integral_def(y, x_lim, y_lim): #Evaluates definite integrals (variable again must be 'x')
    x = sp.Symbol('x') #defines variable of integration
    return sp.integrate(y, (x, x_lim, y_lim))        


def parab_length(parabola, a_lim, b_lim): #returns the length of a specific section of the parabola
    #function_in = equation of the parabola
    #a_lim and b_lim = limits of integration
    
    fprime = deriv(parabola) #takes derivative of parabola equation
    function_int = sp.sqrt(1 + (fprime)**2) #formula for the length of a curve
    return integral_def(function_int, a_lim, b_lim) #evaluates integral


def mirrorline(f_eff): #given an effective focal length, returns equation of mirror line
    f = f_eff/2 #only true for 90 degree reflection
    A = 1/(4*f) #coefficent for parabola equation
    x = sp.Symbol('x') #define x as a symbol in sympy
    
    #Define the point where the line should intersect the parabola (tangent line)
    x_point = f_eff
    y_point = A*(f_eff)**2
    
    #Find slope of the line:
    y_parabola = A*(x**2)
    yprime = deriv(y_parabola)
    
    slope_co = yprime.coeff(x)#extract coefficent
    slope_sym = slope_co*x_point #evaluate yprime at f_eff
    slope = float(slope_sym)
    
    b = y_point - (slope*x_point) #y-intercept/point of reflection
    y_line = (slope*x) + b #final equation of the line in symbol form
    return y_line


def percent_diff(first, second): #returns percent difference between two values
    return 100*np.abs((first - second)/((first + second)/2))



#----Gaussian Propagation Functions----#

def Gauss_op(ABCD, wavelength, w0, d_in):  #d_in is the distance from the input beam waist to the first propagation element
    #given any propagation matrix, can calculate *size and location of output beam waist*
    A = ABCD[0,0]
    B = ABCD[0,1]
    C = ABCD[1,0]
    D = ABCD[1,1]
    
    zc = (pi * (w0**2.))/wavelength #confocal distance/parameter (Goldsmith 2.41)
    
    #distance between the end of the propagation and the output waist (Goldsmith 3.19a)
    d_out = -( ((A*d_in) + B) * ((C*d_in) + D) + (A*C*(zc**2)) ) / ( ((C*d_in) + D)**2 + ((C**2)*(zc**2)) )
    
    #output beam waist size (Goldsmith 3.19b)
    w_out = w0 / ( sqrt( (((C*d_in) + D)**2)  + ((C**2)*(zc**2)) ) )
    
    return d_out, w_out

def freeprop(w0, wavelength, distance): #returns the beam size after free propagation for a chosen distance
    return w0 * sqrt( 1 + ((wavelength * distance) / ( pi * (w0**2))**2 ) )


#----Inputs----#

#Index of Refraction
n = 1.00 #Note: Chose this because index of refraction for a vaccuum: 1.00, for air: 1.0003. 

#Distances [mm]
distM1_in = 140/cos_deg(4.8) #distance the beam travels (since offset by 4.8 degrees for input-output difference) between input and M1
distM1_M2 = 400 #distance between M1 and M2
distM2_DH = 400 #distance between M2 and DH (FOR EQUAL PATH LENGTHS)
distDHb_DHt = inch2mm(3) #distance between the top and bottom of the DH
distM1_out = 140/cos_deg(4.8) #distance between output and M1 (also offset by 4.8 degrees)


#Focal Lengths [mm]
f_effM1 = 140 #effective focal length of M1
f_effM2 = 400 #effective focal length of M2

f_M1 = 76 #primary focal length of M1
f_M2 = 304 #primary focal length of M2

#we have these pairs of focal length values because both M1 and M2 are off-axis paraboloids

                     
#Element Sizes [mm] **all approximate from measurement**
diam_M1 = inch2mm(6.75) #(long edge) #diameter of M1
diam_M2 = inch2mm(7 + (14/16)) #diameter of M2
diam_IPOP = inch2mm(7.65) #diameter of vertical polarizer
diam_BS = inch2mm(8) #diameter of beam splitter
length_DH = inch2mm(7) #length of DH side



#--------Potential TolTEC Inputs Needed-----------#

#Beam Size at Lyot Stop for 1.1 mm:
w_Lyot = 208.6046822 #mm

#Beam Waist after M7 (optical element before Lyot Stop):
w0_Lyot = 13.3114267268 #mm

#Distance between beam waist and M7:
distw0_M7 = -7384.4334282

#Distance between M7 and Lyot Stop
dist2M7 = 530 #mm

LS_diam = 258.6
diam_M7 = 0.33e3

                     
                     
#BIG FUNCTIONS:

def best_w0_finder(wavelength): #returns the best beam waist for any wavelength [mm]

    w0_guess = np.linspace(0.1,25,100000) #an array of random guesses for the beam waist at the input (in mm)
    #might not work if the # of guesses is not large enough, or if the wavelength doesn't fall in this range

    #Calculate Radius of Curvature for M1:
    y, f, ROC_parab, ROC_M1  = newmirror(f_effM1) #the new mirror function works because we use its effective focal length
     
                     
    #Matching with M1's ROC:
    ROC_beam_array = np.zeros(len(w0_guess)) #empty array to hold the random beam ROC values (after running for loop)
                           
    for i in range(len(w0_guess)): #evaluates the ROCs of the beams at M1 for all the different beam waist guesses

        ROC_beam_guess = ROC_beam(wavelength, w0_guess[i], distM1_in) #gives ROC of the beam after propagating forward (we don't use f_eff because the beam waist is 4.8 degrees off the exact focus in order to obtain input/output separation)
        ROC_beam_array[i] = ROC_beam_guess #fills empty ROC_beam_array
        
                     
    #Now we can find where these arrays match in value:
    index_best_tuple = np.where(np.abs(ROC_beam_array - ROC_M1) <= 0.01) #we use np.where to find the index of the optimal ROC and w0 that a light beam could have
    index_best = int(index_best_tuple[0]) #convert to integer


   #Print final values that show how much the ROCs agree (if you choose):
    answer = input("Print detailed values for this wavelength's optimal beam waist? (yes or no) ")
    if answer == "yes" or answer == "Yes":

        print('For', wavelength, 'mm wavelength:')
        print('')
        print('***Best Beam Waist:',w0_guess[index_best],'mm')
        print('Gives:')
        print('ROC of Mirror 1:',ROC_M1,'mm')
        print('ROC of Beam:',ROC_beam_array[index_best],'mm')

        print('')
        print('Difference in ROC:',np.abs(ROC_M1 - ROC_beam_array[index_best]),'mm')
        print('Percent Difference:',percent_diff(ROC_M1,ROC_beam_array[index_best]),'%')
        
    else:
        print(' ')
    
    return w0_guess[index_best] #returns best beam waist
        
    
        

def best_f_eff_finder(wavelength, w0): #w0 is the beam waist that corresponds with the wavelength you choose
    
    #A More Exact Way to Find the Best f_eff:
    #For Loop Method, with Np.Where to Determine the Best f_eff for Matching Beam and Mirror ROC
    
    f_eff_guess = np.linspace(0.1,500,100000) #an array of random guesses for f_eff (in mm)
    
    ROC_mirror_array = np.zeros(len(f_eff_guess)) #empty array to fill with possible mirror ROCS
    ROC_beam_array = np.zeros(len(f_eff_guess)) #empty array to fill with possible beam ROCs
    
    
    for i in range(len(f_eff_guess)): #evaluates the ROC of the beam and mirror at the various f_eff values
    
        y_vals, f_vals, ROC_parab_vals, ROC_mirror_vals  = newmirror(f_eff_guess[i]) #gives ROC of mirror at that f_eff
        ROC_mirror_array[i] = ROC_mirror_vals #fills empty ROC_mirror_array
    
        ROC_beam_vals = ROC_beam(wavelength,w0,f_eff_guess[i]) #gives ROC of the beam at that f_eff
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
    print('Percent Difference:',percent_diff(ROC_mirror_array[index_best],ROC_beam_array[index_best]),'%')

    y_best, f_best, ROC_parab_best, ROC_mirror_best = newmirror(f_eff_guess[index_best])
    print('')
    print('New Mirror Properties:')
    print('y = ',y_best)
    print('Primary Focal Length:',f_best)
    print('Radius of Curvature of Parent Parabola:', ROC_parab_best)
    print('Radius of Curvature of New Mirror:', ROC_mirror_best)
    
    return f_eff_guess[index_best] #returns the best effective focal length
    

      
        
        
#This function returns the distance the beam travels from M2 to the DH when the path lengths aren't equal. 
#Relevant for letting dihedrals move (acting as an interferometer)!

def diffpathlength(distDH_center_new_inches): #Must give the new distance from DH to center of DH axis in inches.
    
    #Note: The separation of the dihedrals will always total 15".
    #So, when the path lengths are equal, both dihedrals are 7.5" away from the DH center.
    #This means that the argument of this function must be between 3" and 12", our limits, (excluding 7.5").
    
    distDH_center_new = inch2mm(distDH_center_new_inches) #converts input to mm
    
    
    #-------Given starting distances-------#:
    distM2_BS = inch2mm(7) #MEASURED: distance between M2 and BS (doesn't change)
    
    
    #When path lengths are equal:
    distM2_DH_eq = distM2_DH #FROM THESIS: distance between M2 and either DH ~when path lengths are equal~
    
    distDH_center_eq = inch2mm(7.5) #MEASURED: distance from a DH to center of DH axis ~when path lengths are equal~

    distBS_DH_eq = np.abs(distM2_DH_eq - distM2_BS) #CALCULATED: distance between BS and DH ~when path lengths are equal~
    
    
    #The path difference between old and new DH offset from center:
    pathdiff = np.abs(distDH_center_new - distDH_center_eq)
    
    
    #Use pythagorean theorem to find the distance between the center of the DH axis to the BS:
    distBS_center = pythag_a(distDH_center_eq,distBS_DH_eq) #(doesn't change)
   

    #-------FIRST LEG CALCULATION-------#:
    #Use pythagorean theorem to find the new distance between the BS and DH:
    distBS_DH_new = pythag_c(distBS_center, distDH_center_new)
    
    
    #In case we want the angle: (not relevant to M2/DH path distance if distM2_BS doesn't change, Grant says true)
    #Use the Law of Cosines to find the angle between the two path lengths:
    angle_diff_rad, angle_diff_deg  = lawofcos_alt(distBS_DH_eq, distBS_DH_new, pathdiff) #must input side lengths in order: 
    #(a,b,c) where: side C = side opposite angle between the two path lengths
    
    
    return distBS_DH_new + distM2_BS #returns the ~new~ distance the beam travels BETWEEN M2 AND DH (not just between BS and DH)

