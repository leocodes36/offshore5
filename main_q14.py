from functionsPy.common import *
from functionsPy.waves import *
from functionsPy.wind import *
from functionsPy.monopile import *
from functionsPy.regularWaves import calculateRegularWaveFrequencyInformation
from functionsPy.integration import ode4
from functionsPy.floaterIntegration import dqdt
from functionsPy.plotting import makeplots
import numpy as np
import pylab as plt
import os

# Output folder setup
of = "outputFig"
os.makedirs(of, exist_ok=True)
def ofy(fileName):
    return os.path.join(of, fileName)

# Close any existing plots
plt.close('all')

# %% ----Q14: RESPONSE TO IRREGULAR WAVES--------------------------------------

# Load the variables
timeInfo = loadFromJSON("inputVariables/time.json")
timeInfo["fHighCut"] = 0.2
constants = loadConstants()
SparBuoyData = loadFromJSON("outputVariables/SparBuoyDataComplete.json")

# Load the rotor data
IEA22MWRotor = loadFromJSON('inputVariables/iea22mw.json')
# Controller parameter & state of the rotor
IEA22MWRotor['ARotor'] = 0.25 * np.pi * IEA22MWRotor['DRotor']**2
IEA22MWRotor['gamma'] = 0.0
IEA22MWRotor['active'] = True

# Wind speed - should be constant
# We set up a spectrum with zero turbulence.
# Just an easy way to reuse existing functions
random_seed_wind = 2
wind = loadFromJSON('inputVariables/wind13.json')
wind = update_struct(wind, timeInfo)  # Ensure this function is defined
wind['t'] = ensure_col_vec(np.arange(0, wind['TDur'], wind['dt']))  # Ensure this function is defined

# FIXME: Add WIND calculations here (use 3 functions if applicable)
# FIXED
wind = calculateKaimalSpectrum(wind)
wind = generateRandomPhases(wind, random_seed_wind)
wind = calculateWindTimeSeriesFFT(wind)

# Vertical locations along floater
z = ensure_col_vec(np.linspace(SparBuoyData['z_Bot'], 0, 100))  # Ensure this function is defined
SparBuoyData['z'] = z

# Calculate the wave kinematics from Irregular wave
waves = loadFromJSON('inputVariables/wave14.json')
waves['z'] = z
random_seed_waves = 1
waves = update_struct(waves, timeInfo)  # Ensure this function is defined
waves['t'] = ensure_col_vec(np.arange(0, waves['TDur'], waves['dt']))  # Ensure this function is defined

# FIXME: Add wave calculations here (use 4 functions if applicable)
# FIXED
waves = calculateJONSWAPSpectrum(waves)
waves = generateRandomPhases(waves, random_seed_waves)
waves = calculateFreeSurfaceElevationTimeSeriesFFT(waves)
waves = calculateKinematicsFFT(waves)

# Response
tode = np.arange(0., timeInfo["TDur"], 2 * timeInfo["dt"])
q0 = np.array([0, 0, 0, 0, np.nan])
q = ode4(dqdt, tode, q0, SparBuoyData, IEA22MWRotor, waves, wind)

response = dict()
response["t"] = tode
response["x1"] = q[:, 0]
response["x5"] = q[:, 1]

# Plotting results
fig14 = makeplots(wind, waves, SparBuoyData, response, timeInfo, 'b')
plt.savefig(ofy("fig14.pdf"))

# Calculate and print standard deviations
surge_std = np.std(q[:, 0])
surge_mean = np.mean(q[:, 0])
pitch_std = np.rad2deg(np.std(q[:, 1]))
pitch_mean = np.rad2deg(np.std(q[:, 1]))

print(f'Q14 Surge Standard deviation [m]: {surge_std}')
print(f'Q14 Surge mean [m]: {surge_mean}')
print(f'Q14 Pitch Standard deviation [deg]: {pitch_std}')
print(f'Q14 Pitch mean [deg]: {pitch_mean}')
