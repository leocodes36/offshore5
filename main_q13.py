from functionsPy.common import *
from functionsPy.waves import *
from functionsPy.wind import *
from functionsPy.regularWaves import *
from functionsPy.integration import ode4
from functionsPy.floaterIntegration import dqdt
from functionsPy.plotting import makeplots
from main_q12 import response as response_waveonly
from main_q12 import wind as no_wind
import numpy as np
import pylab as plt
import os

of = "outputFig"
os.makedirs(of, exist_ok=True)
def ofy(fileName):
    return os.path.join(of, fileName)

plt.close('all')

#%% Q13: response to waves and steady wind

# Load the variables here
timeInfo = loadFromJSON("inputVariables/time.json")
timeInfo["fHighCut"] = 0.2
constants = loadConstants()
SparBuoyData = loadFromJSON("outputVariables/SparBuoyDataComplete.json")
SparBuoyData["CD"] = 0.6; 

z = np.linspace(SparBuoyData["z_Bot"], 0., 100)
SparBuoyData["z"] = z
# Dry decay test
#SparBuoyData["CD"] = 0.

# Wave kinematics - should be zero
waves = loadFromJSON('inputVariables/wave12.json');
waves["z"] = z

with Timer("wind"):
    # Wind speed - should be constant
    # FIXME Modify the input file wind13.json so that it represents a steady wind
    # FIXED
    wind = loadFromJSON('inputVariables/wind13.json');
    wind.update(timeInfo)
    wind["t"] = np.arange(0.,wind["TDur"] ,wind["dt"])
    wind = calculateKaimalSpectrum(wind)
    wind = generateRandomPhases(wind, 1)
    wind = calculateWindTimeSeriesFFT(wind)

# Load the rotor
IEA22MWRotor = loadFromJSON('inputVariables/iea22mw.json')
IEA22MWRotor["ARotor"] = 0.25*np.pi*IEA22MWRotor["DRotor"]**2

# Controller parameter & state of the rotor
IEA22MWRotor["gamma"] = 0.
IEA22MWRotor["active"] = True

with Timer("waves"):
    # calculate the wave kinematics - zero at this stage
    waves.update(timeInfo)
    waves["t"] = np.arange(0.,wind["TDur"] ,wind["dt"])
    waves = calculateRegularWaveFrequencyInformation(waves)
    waves = calculateFreeSurfaceElevationTimeSeries(waves)
    waves = calculateKinematics(waves)

# Integration time array
tode = np.arange(0., timeInfo["TDur"], 2*timeInfo["dt"])

# Response
# Pitch decay with wind
with Timer("Integration"):
    # Initial conditions for surge decay
    q0 = np.array([0,0,0,0,np.nan])
    q = ode4(dqdt, tode, q0, SparBuoyData, IEA22MWRotor, waves, wind)

response = dict()
response["t"] = tode;
response["x1"] = q[:,0]
response["x5"] = q[:,1]

fig13 = makeplots(wind, waves, SparBuoyData, response, timeInfo, 'g');
fig13 = makeplots(no_wind, waves, SparBuoyData, response_waveonly, timeInfo, 'b', ax=fig13)
fig13[0,0].legend(["combined", "wave-only"])
plt.savefig(ofy("fig13.pdf"))

print(f'Q13 Surge Standard deviation [m]: {np.std(q[:,0])}');
print(f'Q13 Surge mean [m]: {np.mean(q[:,0])}');
print(f'Q13 Pitch Standard deviation [deg]: {np.rad2deg(np.std(q[:,1]))}')
print(f'Q13 Pitch mean [deg]: {np.rad2deg(np.mean(q[:,1]))}')