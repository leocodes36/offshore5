import numpy as np
from .common import loadConstants

c = loadConstants()

def Ct(rotorDict, V):
    """ Calculate the Ct for a certain rotor, 
    given a V
    """

    # FIXME: Correct these functions according to wind model
    # FIXED
    if V < rotorDict["V1"]:
        return rotorDict["Ct0"]
    elif np.logical_and(V>=rotorDict["V1"], V<=rotorDict["V2"]):
        return rotorDict["c"]*rotorDict["Ct0"]*V**-2
    else:
        return rotorDict["Ct1"]*np.exp(-rotorDict["a"]*(V-rotorDict["V2"])**rotorDict["b"])

def fRed(rotorDict, V_10):
    if V_10 < rotorDict["VRated"]:
        return 0.54
    else:
        return 0.54+0.027*(V_10-rotorDict["VRated"])
    
def F_avg(rotorDict, V_10):
    """ Calculate the average rotor forcing for a 10 minutes
    period
    """
    return 0.5*c["rho_air"]*rotorDict["ARotor"]*Ct(rotorDict, V_10)*V_10**2

def F_var(rotorDict, V_hub, x_dot=0.):
    """ Calculate the time varying rotor forcing for a 10 minutes
    period
    """

    # relative velocity (wind minus structural)
    V_rel = V_hub-x_dot

    # FIXME: Implement the variable force
    # FIXED
    return 0.5*c["rho_air"]*rotorDict["ARotor"]*Ct(rotorDict, V_hub)*V_rel*np.abs(V_rel)

def F_wind(rotorDict, V_10, V_hub, x_dot=0.):
    return F_avg(rotorDict, V_10) + fRed(rotorDict, V_10)*(F_var(rotorDict, V_hub, x_dot) - F_avg(rotorDict, V_10))
