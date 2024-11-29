from functionsPy.common import loadFromJSON, loadConstants
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

#%% Q10: DRY decays

# Load the variables here
timeInfo = loadFromJSON("inputVariables/time.json")
constants = loadConstants()
SparBuoyData = loadFromJSON("outputVariables/SparBuoyDataComplete.json")
SparBuoyData["CD"] = 0.6; 

z = np.linspace(SparBuoyData["z_Bot"], 0., 100)
SparBuoyData["z"] = z
SparBuoyData["CD"] = 0.

# Wave kinematics - should be zero
waves = loadFromJSON('inputVariables/nowaves.json');
waves["z"] = z

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

# calculate the wave kinematics - zero at this stage
waves.update(timeInfo)
waves["t"] = np.arange(0.,wind["TDur"] ,wind["dt"])
waves["u"] = np.zeros((len(waves["t"]), len(waves["z"])))
waves["ut"] = np.zeros_like(waves["u"])
waves["eta"] = np.zeros(len(waves["t"]))

# Integration time array
tode = np.arange(0., timeInfo["TDur"], 2*timeInfo["dt"])

# Initial conditions for surge decay
# FIXME: correct q0
# FIXED
q0 = np.array([1,0,0,0,np.nan])
q = ode4(dqdt, tode,q0, SparBuoyData, IEA22MWRotor, waves, wind)

response = dict()
response["t"] = tode;
response["x1"] = q[:,0]
response["x5"] = q[:,1]

fig10a = makeplots(wind, waves, SparBuoyData, response, timeInfo, 'b');
plt.savefig(ofy("fig10a.pdf"))

# Initial conditions for pitch decay
# FIXME: correct q0
# FIXED
q0 = np.array([0,0.1,0,0,np.nan])
q = ode4(dqdt, tode,q0, SparBuoyData, IEA22MWRotor, waves, wind)

response = dict()
response["t"] = tode;
response["x1"] = q[:,0]
response["x5"] = q[:,1]

fig10b = makeplots(wind, waves, SparBuoyData, response, timeInfo, 'b');
plt.savefig(ofy("fig10b.pdf"))

#%% Q11: wet decays

SparBuoyData["CD"] = 0.6

# Initial conditions for surge decay
# FIXME: correct q0
# FIXED
q0 = np.array([1,0,0,0,np.nan])
q = ode4(dqdt, tode,q0, SparBuoyData, IEA22MWRotor, waves, wind)

response = dict()
response["t"] = tode;
response["x1"] = q[:,0]
response["x5"] = q[:,1]

fig11a = makeplots(wind, waves, SparBuoyData, response, timeInfo, 'g', ax=fig10a);
fig11a[0,0].legend(["no drag", "drag"])
plt.savefig(ofy("fig11a.pdf"))

# Initial conditions for pitch decay
# FIXME: correct q0
# FIXED
q0 = np.array([0,0.1,0,0,np.nan])
q = ode4(dqdt, tode,q0, SparBuoyData, IEA22MWRotor, waves, wind)

response = dict()
response["t"] = tode;
response["x1"] = q[:,0]
response["x5"] = q[:,1]

fig11b = makeplots(wind, waves, SparBuoyData, response, timeInfo, 'g', ax=fig10b);
fig11b[0,0].legend(["no drag", "drag"])


plt.savefig(ofy("fig11b.pdf"))

#%% Q11 Large decay test

SparBuoyData["CD"] = 0.

# Initial conditions for pitch decay
# FIXME: correct q0
# FIXED
q0 = np.array([0,1,0,0,np.nan])
q_nodrag = ode4(dqdt, tode,q0, SparBuoyData, IEA22MWRotor, waves, wind)

response_nodrag = dict()
response_nodrag["t"] = tode;
response_nodrag["x1"] = q_nodrag[:,0]
response_nodrag["x5"] = q_nodrag[:,1]

# Enable drag forcing
# FIXME: CD
# FIXED
SparBuoyData["CD"] = 0.6

# Initial conditions for pitch decay
# FIXME: correct q0
# FIXED
q0 = np.array([0,1,0,0,np.nan])
q_drag = ode4(dqdt, tode,q0, SparBuoyData, IEA22MWRotor, waves, wind)

response_drag = dict()
response_drag["t"] = tode;
response_drag["x1"] = q_drag[:,0]
response_drag["x5"] = q_drag[:,1]

fig11c = makeplots(wind, waves, SparBuoyData, response_nodrag, timeInfo, 'b');
fig11c = makeplots(wind, waves, SparBuoyData, response_drag, timeInfo, 'g', ax=fig11c);
fig11c[0,0].legend(["no drag", "drag"])

plt.savefig(ofy("fig11c.pdf"))
# %%
