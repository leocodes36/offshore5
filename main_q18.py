import numpy as np
import matplotlib.pyplot as plt
import os
from functionsPy.common import *
from functionsPy.waves import *
from functionsPy.wind import *
from functionsPy.integration import ode4
from functionsPy.floaterIntegration import dqdt
from functionsPy.plotting import makeplots

# Output folder setup
of = "outputFig"
os.makedirs(of, exist_ok=True)

def ofy(fileName):
    return os.path.join(of, fileName)

# Close any existing plots
plt.close('all')

# %% ----Q18: RESPONSE TO GAMMA PARAMETERS--------------------------------------

# Load the time information
timeInfo = loadFromJSON('inputVariables/time.json')
timeInfo["fHighCut"] = 0.2

# Wind speed - should be constant
# Just a loop to handle both wind cases
winds = ['wind16B', 'wind16B', 'wind16B', 'wind16B']
colors = ['b', 'g', 'r', 'm']

# Gamma values for the four cases
gammas = [2., 0.5, 0.2, 0.01]
labels = ['gamma = 2', 'gamma = 0.5', 'gamma = 0.2', 'gamma = 0.01']

# Initialize the figure outside the loop
fig18 = plt.figure()
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
    IEA22MWRotor['gamma'] = gammas[i]
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
    q0 = np.array([0, 0, 0, 0, 0])
    q = ode4(dqdt, tode, q0, SparBuoyData, IEA22MWRotor, waves, wind)

    response = dict()
    response["t"] = tode
    response["x1"] = q[:, 0]
    response["x5"] = q[:, 1]

    # Plot the results on the same figure
    if i == 0:
        fig18 = makeplots(wind, waves, SparBuoyData, response, timeInfo, colors[i])
    if i == 1:
        fig18 = makeplots(wind, waves, SparBuoyData, response, timeInfo, colors[i], ax=fig18)
    if i == 2:
        fig18 = makeplots(wind, waves, SparBuoyData, response, timeInfo, colors[i], ax=fig18)
    if i == 3:
        fig18 = makeplots(wind, waves, SparBuoyData, response, timeInfo, colors[i], ax=fig18)

    # Print the standard deviations
    print(f'Q18 Surge Standard deviation [m]: {np.std(q[:,0])}')
    print(f'Q18 Surge mean [m]: {np.mean(q[:,0])}')
    print(f'Q18 Pitch Standard deviation [deg]: {np.rad2deg(np.std(q[:,1]))}')
    print(f'Q18 Pitch mean [deg]: {np.rad2deg(np.mean(q[:,1]))}')

# Adjust the legend after the loop
# Get handles of all lines in the figure
#allLines = [line for line in fig18.lines]  # Assuming makeplots adds lines sequentially

# Select only the main lines corresponding to the loop iterations
#lineHandles = allLines[-len(winds):]  # Assuming the lines are added sequentially

# Add labels and legend
#plt.legend(lineHandles, labels)  # Use correct labels
#plt.xlabel('Time (s)')
#plt.ylabel('Response')
#plt.hold(False)  # Stop adding to the current plot
fig18[0,0].legend(labels)
# Save the figure
#fig18.set_name('Q18')
plt.savefig(ofy('fig18.pdf'), format='pdf')

# Close the plot
#plt.close(fig18)
