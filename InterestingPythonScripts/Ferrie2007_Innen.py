import numpy as np
import matplotlib.pyplot as plt
import pylab
import scipy.integrate as integrate

x= -500
y= -500
z = 10

# Konstanten fuer CMZ
xc =-50			# Position Mitte in allg Koordinaten
yc = 50
TettaC = 70

#Konstanten fuer DISK
alpha = 13.5
beta = 20.
TettaD = 48.5

# Abmessungen in CMZ Koordinaten
XMAX=250		
XC = XMAX/2
LC = XMAX/(2*np.log(2)**0.25)
HC = 18.
HC2 = 54.

# Abmessungen in DISK Koordinaten
XD = 1200
LD = 438.
HD = 42.
HD2 = 120.

#Konstanten fuer HII -WIM-
y3 = -10
z3= -20
L3 = 145.
H3 = 26.
L2 = 3700.
H2 = 140.
L1 = 17000
H1=950.

#Konstanen fuer HII VHIM
alphaVH = 21
LVH=162
HVH = 90


def Bogenmass(x):			# Trafo ins Bogenmass fuer Winkel zur Berechnung
	return x*np.pi/180

def cos(x):					# Cos FKT fuer Gradmass
	x=Bogenmass(x)
	return np.cos(x)
def sin(x): 				# Sin FKT fuer Gradmass
	x=Bogenmass(x)
	return np.sin(x)
def sech2(x):
	return np.cosh(x)**2
def u(x):
	if x.all<0:
		return 0
	else:
		return 1	

def CMZ_X_Trafo(x,y):
	return (x-xc)*cos(TettaC) +(y-yc)*sin(TettaC)
	
def CMZ_Y_Trafo(x,y):
	return -(x-xc)*sin(TettaC) +(y-yc)*cos(TettaC)

def DISK_X_Trafo(x,y,z):
	return x*cos(beta)*cos(TettaD) - y*(sin(alpha)*sin(beta)*cos(TettaD) -cos(alpha)*sin(TettaD))-z*(cos(alpha)*sin(beta)*cos(TettaD) +sin(alpha)*sin(TettaD))

def DISK_Y_Trafo(x,y,z):
	xT= x*cos(beta)*sin(TettaD)
	yT = y*(sin(alpha)*sin(beta)*sin(TettaD) +cos(alpha)*cos(TettaD))
	zT = z*(cos(alpha)*sin(beta)*sin(TettaD) -sin(alpha)*sin(TettaD))
	return -xT+yT+zT

def DISK_Z_Trafo(x,y,z):
	xT = x*sin(beta)
	yT = y*sin(alpha)*cos(beta)
	zT = z*cos(alpha)*cos(beta)
	return xT+yT+zT
	
#Mollekularer Wasserstoff im CMZ,
def n_H2_CMZ(x0,y0,z0): 			# Eingabe in Urspruenglichen koordinaten
	x = CMZ_X_Trafo(x0,y0)
	y = CMZ_Y_Trafo(x0,y0)
	XY_Help = ((np.sqrt(x**2+(2.5*y)**2)-XC)/LC)**4
	return 150*np.exp(-XY_Help)*np.exp(-(z0/HC)**2)

#Atomarer Wasserstoff im CMZ 
def n_HI_CMZ(x0,y0,z0):			#Eingabe in Urspruenglichen Koordinaten
	x=CMZ_X_Trafo(x0,y0)
	y=CMZ_Y_Trafo(x0,y0)
	A=np.sqrt(x**2 +(2.5*y)**2)
	B= (A-XC)/LC
	XY_Help=B**4
	Z = (z0/HC2)**2
	return 8.8*np.exp(-XY_Help)*np.exp(-Z)
	
#Mollekularer Wasserstoff in der DISK
def n_H2_DISK(x0,y0,z0):
	x= DISK_X_Trafo(x0,y0,z0)
	y= DISK_Y_Trafo(x0,y0,z0)
	z=DISK_Z_Trafo(x0,y0,z0)
	return 4.8*np.exp(-((np.sqrt(x**2 + (3.1*y)**2) - XD)/LD)**4)*np.exp(-(z/HD)**2)

#Atomarer Wasserstoff in der DISK
def n_HI_DISK(x0,y0,z0):
	x= DISK_X_Trafo(x0,y0,z0)
	y= DISK_Y_Trafo(x0,y0,z0)
	z=DISK_Z_Trafo(x0,y0,z0)
	return 0.34*np.exp(-((np.sqrt(x**2 + (3.1*y)**2) - XD)/LD)**4)*np.exp(-(z/HD2)**2)

#Ioniesierter Wasserstoff
def n_HII_WIM(x0,y0,z0):
	r=np.sqrt(x0**2+y0**2+z0**2)
	P1 = np.exp(-(x**2+(y0-y3)**2)/L3**2)*np.exp(-(z0-z3)**2/H3**2)
	P2 = np.exp(-((r-L2)/(0.5*L2))**2)*sech2(z/H2)
	P3 = np.cos(np.pi*r*0.5/L1)*sech2(z/H1)
	return 8.0*(P1+0.009*P2+0.005*P3)

def n_HII_VHIM(x0,y0,z0):
	e = y0*cos(alphaVH)+z0*sin(alphaVH)
	s = -y0*sin(alphaVH) + z*cos(alphaVH)
	return 0.29*np.exp(-((x0**2+e**2)/LVH**2 + s**2/HVH**2))

def n_HII(x0,y0,z0):
	return n_HII_VHIM(x0,y0,z0) +n_HII_WIM(x0,y0,z0)

def n_HI(x,y,z):
	return n_HI_DISK(x,y,z) + n_HI_CMZ(x,y,z)

def n_H2(x,y,z):
	return n_H2_CMZ(x,y,z) + n_H2_DISK(x,y,z)
