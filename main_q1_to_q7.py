from functionsPy.common import loadFromJSON, saveToJSON, loadConstants
import numpy as np
import os

# Load the variables here
timeInfo = loadFromJSON("inputVariables/time.json")
constants = loadConstants()
SparBuoyData = loadFromJSON("inputVariables/SparBuoyData.json")

g = constants['g']
rhow = constants['rho_water']
rhoa = constants['rho_air']
rhos = SparBuoyData['rho_Steel']
mtu = SparBuoyData['M_Turbine']
zturb = SparBuoyData['z_CM_Turbine']
zhub = SparBuoyData['z_Hub']
fb = SparBuoyData['fb']
draft = SparBuoyData['draft']
Dspar = SparBuoyData['DMonopile']
th = SparBuoyData['Thickness']
Dhp = Dspar
Kmoor = SparBuoyData['K_Moor']
zmoor = SparBuoyData['z_Moor']
Cm = SparBuoyData['Cm']
CM = SparBuoyData['CM']
Cd = SparBuoyData['CD']
mt = SparBuoyData['M_Tower']
zCMt = SparBuoyData['z_CM_Tower']
ICMt = SparBuoyData['I_CM_Tower']
BallastHeight = SparBuoyData['BallastHeightindraft']
BallastCOG = SparBuoyData['Ballast_COG']
mb = SparBuoyData['M_Ballast']


IEA22MWRotor = loadFromJSON('inputVariables/iea22mw.json')

B11 = SparBuoyData["B11"]; 
thrustr = 2*SparBuoyData["MaxThrust"]

##% Preliminary computations
# Calculate center of buoyancy, center of mass, floater inertias
# Location of floater bottom
zbot = -draft;

zballst = BallastCOG; # height of draft is ballast 

# FIXME: Center of buoyancy
# FIXED
zCB = 0.5 * zbot

# FIXME: Displacement volume
# FIXED
Vol = np.pi/4 * Dspar**2 * zbot

# FIXME: Spar length
# FIXED
ls = draft + fb

# FIXME: Spar mass without ballast
# FIXED
ms = SparBuoyData["M_Floater"]

# Floater mass with ballast
mf = ms + mb

# FIXME: Spar center of mass without ballast
# FIXED
zCMs = fb - 0.5 * ls

# FIXME: Floater center of mass with ballast
# FIXED
zCMb = draft * (-1 + 0.025)
zCMf = (ms * zCMs + mb * zCMb)/(mf)

# FIXME: Spar inertia about its Center of Mass without ballast
# FIXED
ICMs = SparBuoyData["I_CM_Floater"]

# FIXME: Floater inertia about floater CM with ballast
# FIXED
ICMf = ICMs + ms * (zCMs - zCMf)**2 + mb * (zCMf - zCMb)**2



##% Q6: System matrices

# FIXME: Total mass
# FIXED
mtot = mf + mt + mtu
print(mtot/10**7) # CHECKS OUT

# FIXME: Total center of mass
# FIXED
zCMtot = (mf * zCMf + mt * zCMt + mtu * SparBuoyData["z_CM_Turbine"]) / mtot
print(zCMtot)
# FIXME: Total inertia about flotation point
# FIXED
IOtot = ICMf + mf * zCMf**2 + ICMt + mt * zCMt**2 + mtu * SparBuoyData["z_CM_Turbine"]**2
print(IOtot/10**11)
# FIXME
# FIXED
M = np.array([[mtot, mtot*zCMtot],[mtot*zCMtot, IOtot]])

A = np.array([[mtot, -1*rhow*np.pi/4*Dspar**2*Cm*(0.5*zbot**2)], [-1*rhow*np.pi/4*Dspar**2*Cm*(0.5*zbot**2), -1*rhow*np.pi/4*Dspar**2*Cm*((1/3)*zbot**3)]])

B = np.array([[B11, 0.],[0. ,0.]])

# FIXME: Water Plane Inertia
# FIXED
IAA = (Dspar**4)*np.pi/64

# FIXME: hydrodynamic stiffness
# FIXED
Chst = np.array([[0., 0.], [0, mtot*g*(zCB-zCMtot)+rhow*g*IAA]])
# Mooring restoring matrix
Cmoor = np.array([[Kmoor, Kmoor*zmoor], [Kmoor*zmoor, Kmoor*(zmoor**2)]])

C = Chst + Cmoor;

SparBuoyData["M"] = M;
SparBuoyData["C"] = C;
SparBuoyData["A"] = A;
SparBuoyData["B"] = B;

print(M)
print(A)
print(C)

##% Natural Frequencies

# FIXME: calculate C over MA
# FIXED
MA = M + A
CoMA = np.dot(np.linalg.inv(MA), C)
eigVal, eigVec = np.linalg.eig(CoMA)

# Natural frequencies
omeganat = np.sqrt(eigVal); 
fnat = omeganat/2/np.pi;

SparBuoyData["fnat"] = fnat;

# Natural periods
Tnat = 1./fnat;

# Display surge and pitch natural periods
print(f'Surge period: {Tnat[0]:.2f} [s]')
print(f'Pitch period: {Tnat[1]:.2f} [s]')

# FIXME: Added mass in heave
# A33 = ...; 
a33 = 0.5 * rhow * 0.64 * 4/3 * np.pi * (Dhp/2)**3

# FIXME: Hydrostatic restoring in heave
c33 = rhow * g * np.pi/4 * Dhp**2

# FIXME: Heave natural period
Theave = 1;

print(f'Heave period: {Theave:.2f} [s]')
os.makedirs("outputVariables", exist_ok=True)
saveToJSON(SparBuoyData, "outputVariables/SparBuoyDataComplete.json")