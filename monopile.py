import numpy as np
from .common import loadConstants

g = loadConstants()["g"]
rho_water = loadConstants()["rho_water"]

def forceIntegrate(monopileDict, u, ut, z, x_dot):
    
    h = np.abs(z[0])
    u = u - x_dot
    
    df = forceDistributed(monopileDict, u, ut, z, x_dot)
    
    F = np.trapz(df, z)
    M = np.trapz(df*(z + h), z)
    
    return F, M

def forceDistributed(monopileDict, u, ut, z, x_dot):
    
    u = u - x_dot
    # FIXME:  add back the inertia forces
    # FIXED
    df = 0.5*rho_water*monopileDict["DMonopile"]*monopileDict["CD"]*np.abs(u)*u+rho_water*monopileDict["CM"]*(1/4*np.pi*monopileDict["DMonopile"]**2)*ut
    
    return df

def computeElementwiseQuantities(monopileDict):
    
    outputDict = dict()
    outputDict.update(monopileDict)
    
    # Compute missing element properties
    z = monopileDict["zBeamNodal"]
    dz = np.diff(z)
    outputDict["zBeamElement"] = z[:-1] + dz/2
    outputDict["dz"] = dz
    
    # Compute the phiNodal
    phi = monopileDict["phiNodal"]
    dPhi = np.diff(phi)
    outputDict["phiElement"] = outputDict["phiNodal"][:-1] + dPhi/2
    
    return outputDict
    