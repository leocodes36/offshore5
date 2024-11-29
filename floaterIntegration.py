from .integration import lookup
import numpy as np
from .loads import forceDistributed
from  .rotor import *
from bisect import bisect_left as lookup
from .floatingRotor import F_wind as F_wind_floating

def dqdt(t, q,
                structure,
                rotor,
                waves,
                wind):
    
    x1 = np.array(q[0:2])
    xdot1 = np.array(q[2:4])
    CT1 = q[4]
    
    rotor["CT"] = CT1

    structure["B"] = np.array(structure["B"])
    structure["C"] = np.array(structure["C"])
    structure["M"] = np.array(structure["M"])
    structure["A"] = np.array(structure["A"])
    
    # Extract time index
    i_ = lookup(waves["t"], t)
    
    # Read wind speed
    V_hub = wind["V_hub"][i_]
    V_10 = wind["V_10"]
    
    # Nacelle speed
    x_dot_rotor = xdot1[0] + structure["zhub"]*xdot1[1];

    # Wind force if the rotor is on
    Thrust, CTVrel = F_wind_floating(rotor, V_10, V_hub, x_dot_rotor)
    Faero = np.array([Thrust, Thrust*structure["zhub"]])
    
    # FIXME for xdot1
    # FIXED
    x_dot_submerged = xdot1[0] + structure["z"]*xdot1[1] # See assignment question 2
    
    u, ut = waves["u"][i_,:], waves["ut"][i_,:]
    df = forceDistributed(structure, u, ut, structure["z"], x_dot_submerged)
    Fhydro = np.array([np.trapz(df, structure["z"]),
                       np.trapz(df*structure["z"], structure["z"])])
    
    output = np.zeros(5)
    output[0:2] = xdot1
    # FIXME for Equation of motion below
    # FIXED
    F = Fhydro + Faero
    output[2:4] = np.dot((F.T - np.dot(structure["B"],xdot1.T) - np.dot(structure["C"], x1.T)), np.linalg.inv(structure["M"] + structure["A"]))
    # FIXME for control equation below
    output[4] = -rotor["gamma"] * (CT1 - CTVrel)
    
    return output