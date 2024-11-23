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
Cm = SparBuoyData['CM']
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
zCB = 0.5 * zbot

# FIXME: Displacement volume
Vol = np.pi/4 * Dspar**2 * zbot

# FIXME: Spar length
ls = draft + fb

# FIXME: Spar mass without ballast
ms = SparBuoyData["M_Floater"]

# Floater mass with ballast
mf = ms + mb

# FIXME: Spar center of mass without ballast
zCMs = fb - 0.5 * ls

# FIXME: Floater center of mass with ballast
zCMb = draft * (-1 + 0.025)
zCMf = (ms * zCMs + mb * zCMb)/(mf)

# FIXME: Spar inertia about its Center of Mass without ballast
ICMs = SparBuoyData["I_CM_Floater"]

# FIXME: Floater inertia about floater CM with ballast
ICMf = ICMs + ms * (zCMs - zCMf)**2 + mb * (zCMf - zCMb)**2



##% Q6: System matrices

# FIXME: Total mass
mtot = mf + mt + mtu
print(mtot/10**7) # CHECKS OUT

# FIXME: Total center of mass
zCMtot = (mf * zCMf + mt * zCMt + mtu * SparBuoyData["z_CM_Turbine"]) / mtot
print(zCMtot)
# FIXME: Total inertia about flotation point
IOtot = ICMf + mf * zCMf**2 + ICMt + mt * zCMt**2 + mtu * SparBuoyData["z_CM_Turbine"]**2
print(IOtot/10**11)
# FIXME
M = np.array([[mtot, mtot*zCMtot],[mtot*zCMtot, IOtot]])

A = np.array([[mtot, -rhow*np.pi/4*Dspar**2*Cm*(0.5*zbot**2)], [-rhow*np.pi/4*Dspar**2*Cm*(0.5*zbot**2), -rhow*np.pi/4*Dspar**2*Cm*(0.333*zbot**3)]])

B = np.array([[B11, 0.],[0. ,0.]])

# FIXME: Water Plane Inertia
IAA = (Dspar**4)*np.pi/64

# FIXME: hydrodynamic stiffness
Chst = np.array([[0., 0.], [0, -mtot*g*(0.5*zbot-zCMf)-rhow*g*IAA]])

# Mooring restoring matrix
Cmoor = np.array([[Kmoor, Kmoor*zmoor], [Kmoor*zmoor, Kmoor*zmoor**2]])

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
MA = np.dot(M, A)
CoMA = np.dot(C, np.linalg.inv(MA))
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
a33 = 0.5

# FIXME: Hydrostatic restoring in heave
c33 = 1;

# FIXME: Heave natural period
Theave = 1;

print(f'Heave period: {Theave:.2f} [s]')
os.makedirs("outputVariables", exist_ok=True)
saveToJSON(SparBuoyData, "outputVariables/SparBuoyDataComplete.json")