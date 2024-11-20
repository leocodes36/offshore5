from functionsPy.common import loadFromJSON, loadConstants
from functionsPy.waves import calculateFreeSurfaceElevationTimeSeries, calculateKinematics
from functionsPy.regularWaves import calculateRegularWaveFrequencyInformation
from functionsPy.integration import ode4
from functionsPy.floaterIntegration import dqdt
from functionsPy.plotting import makeplots
import numpy as np
import pylab as plt
import os

of = "outputFig"
os.makedirs(of, exist_ok=True)
def ofy(fileName):
    return os.path.join(of, fileName)

plt.close('all')

#%% Q12: Response to regular waves

# Load the variables here
timeInfo = loadFromJSON("inputVariables/time.json")
constants = loadConstants()
SparBuoyData = loadFromJSON("outputVariables/SparBuoyDataComplete.json")

# FIXME set correct CD
SparBuoyData["CD"] = 0.;

z = np.linspace(SparBuoyData["z_Bot"], 0., 100)
SparBuoyData["z"] = z
SparBuoyData["CD"] = 0.

# Wave kinematics 
waves = loadFromJSON('inputVariables/wave12.json');
waves["z"] = z
waves.update(timeInfo)

waves["t"] = np.arange(0.,timeInfo["TDur"] ,timeInfo["dt"])

# FIXME : fix calculateRegularWaveFrequencyInformation to make this work
waves = calculateRegularWaveFrequencyInformation(waves)
waves = calculateFreeSurfaceElevationTimeSeries(waves)
waves = calculateKinematics(waves)

# Wind speed - should be zero
wind = loadFromJSON('inputVariables/nowind.json');
wind.update(timeInfo)
wind["t"] = np.arange(0.,wind["TDur"] ,wind["dt"])
wind["V_hub"] = np.zeros_like(wind["t"])

# Load the rotor
IEA22MWRotor = loadFromJSON('inputVariables/iea22mw.json')
IEA22MWRotor["ARotor"] = 0.25*np.pi*IEA22MWRotor["DRotor"]**2

# Controller parameter & state of the rotor
IEA22MWRotor["gamma"] = 0.
IEA22MWRotor["active"] = False

# Integration time array
tode = np.arange(0., timeInfo["TDur"], 2*timeInfo["dt"])

# FIXME: q0 for pitch decay
q0 = np.array([0,0,0,0,np.nan])
q = ode4(dqdt, tode, q0, SparBuoyData, IEA22MWRotor, waves, wind)

response = dict()
response["t"] = tode;
response["x1"] = q[:,0]
response["x5"] = q[:,1]

fig12 = makeplots(wind, waves, SparBuoyData, response, timeInfo, 'b');
plt.savefig(ofy("fig11.pdf"))

print(f'Q12 Surge Standard deviation [m]: {np.std(q[:,0])}');
print(f'Q12 Pitch Standard deviation [deg]: {np.rad2deg(np.std(q[:,1]))}')