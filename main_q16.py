from functionsPy.common import *
from functionsPy.waves import *
from functionsPy.wind import *
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

# %% ----Q16: RESPONSE TO WIND--------------------------------------

# Load the time information
timeInfo = loadFromJSON('inputVariables/time.json')
timeInfo["fHighCut"] = 0.2

# Wind speed - should be constant
# Just a loop to handle both wind cases
winds = ['wind16a', 'wind16b']
colors = ['b', 'g']
labels = ['10 m/s', '16 m/s']

# Initialize the figure outside the loop
fig16 = plt.figure()
#plt.hold(True)  # Keep the plot open for multiple plots

for i in range(len(winds)):
    wind = loadFromJSON(f'inputVariables/{winds[i]}.json')
    wind = update_struct(wind, timeInfo)  # Ensure this function is defined
    wind['t'] = ensure_col_vec(np.arange(0, wind['TDur'], wind['dt']))  # Ensure this function is defined
    wind = calculateKaimalSpectrum(wind)  # Ensure this function is defined
    wind = generateRandomPhases(wind, 1)  # Ensure this function is defined
    wind = calculateWindTimeSeriesFFT(wind)  # Ensure this function is defined

    # Load the rotor and other necessary parameters
    IEA22MWRotor = loadFromJSON('inputVariables/iea22mw.json')
    IEA22MWRotor['ARotor'] = 0.25 * np.pi * IEA22MWRotor['DRotor']**2
    IEA22MWRotor['gamma'] = 0.0
    IEA22MWRotor['active'] = True

    # Disable drag forcing
    SparBuoyData = loadFromJSON("outputVariables/SparBuoyDataComplete.json")
    SparBuoyData['CD'] = 0  # Dry decay test

    # Vertical locations along floater
    z = ensure_col_vec(np.linspace(SparBuoyData['z_Bot'], 0, 100))  # Ensure this function is defined
    SparBuoyData['z'] = z

    # Calculate the wave kinematics from Irregular wave
    waves = loadFromJSON('inputVariables/wave14.json')
    waves['z'] = z
    random_seed_waves = 1
    waves = update_struct(waves, timeInfo)  # Ensure this function is defined
    waves['t'] = ensure_col_vec(np.arange(0, wind['TDur'], wind['dt']))  # Ensure this function is defined
    waves['u'] = np.zeros((len(waves['t']), len(waves['z'])))
    waves['ut'] = np.zeros((len(waves['t']), len(waves['z'])))
    waves['eta'] = ensure_col_vec(np.zeros((len(waves['t']), 1)))

    # Response
    tode = np.arange(0., timeInfo["TDur"], 2 * timeInfo["dt"])
    q0 = np.array([0, 0, 0, 0, np.nan])
    q = ode4(dqdt, tode, q0, SparBuoyData, IEA22MWRotor, waves, wind)

    response = dict()
    response["t"] = tode
    response["x1"] = q[:, 0]
    response["x5"] = q[:, 1]

    # Plot the results on the same figure
    if i == 0:
        fig16 = makeplots(wind, waves, SparBuoyData, response, timeInfo, colors[i])
    else:
        fig16 = makeplots(wind, waves, SparBuoyData, response, timeInfo, colors[i], ax=fig16)

    # Print the standard deviations
    print(f'Q16 Surge Standard deviation [m]: {np.std(q[:,0])}')
    print(f'Q16 Surge mean [m]: {np.mean(q[:,0])}')
    print(f'Q16 Pitch Standard deviation [deg]: {np.rad2deg(np.std(q[:,1]))}')
    print(f'Q16 Pitch mean [deg]: {np.rad2deg(np.mean(q[:,1]))}')

# Add labels and legend
fig16[0,0].legend(labels)  # Use labels from the loop

#plt.hold(False)  # Stop adding to the current plot

# Save the figure
#fig16.set_name('Q16')
plt.savefig(ofy('fig16.pdf'), format='pdf')

# Close the plot
#plt.close(fig16)
