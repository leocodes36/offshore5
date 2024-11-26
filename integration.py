import numpy as np
from  .monopile import *
from  .rotor import *
from bisect import bisect_left as lookup

def rk4(f, h, y0, t0, *args, **kwargs):
    """
    h => dt 
    y0 => q0
    t0 => starting point for integration, from which we move with h steps
    args, kwargs => more arguments to the integrand function
    """
    k1 = f(t0, y0, *args, **kwargs)
    k2 = f(t0 + h / 2, y0 + h / 2 * k1, *args, **kwargs)
    k3 = f(t0 + h / 2, y0 + h / 2 * k2,*args, **kwargs)
    k4 = f(t0 + h, y0 + h * k3, *args, **kwargs)
    return y0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6

def ode4(odefun, tspan, y0, *args, **kwargs):
    """
    Python implementation of the rk4 integrator

    Parameters:
    ----------

    odefun: function handle
        Function to be integrated
    tspan: numpy array
        Time span
    y0: numpy array
        Initial conditions
    
    Returns:
    --------
    y : numpy array
        Integrated function

    """
    y = np.zeros((len(tspan),len(y0))) # array to store the solution
    y[0, :] = y0 # save the initial condition at t0
    h = np.diff(tspan) # time step array
    y_ini = y0 # store the initial condition for the first time step

    for i_, t_ in enumerate(tspan[:-1]):
        ti = tspan[i_]
        hi = h[i_]
        y_integrated = rk4(odefun, hi, y_ini, ti, *args, **kwargs)
        y_ini = y_integrated # save the initial condition for the next loop
        y[i_+1, :] = y_integrated # save the result in the solution array

    return y

def dqdt( t, q,
                structure,
                rotor,
                waves,
                wind):

        alpha, alphaDot = q

        GF = GFCalc(t,
                        alphaDot,
                        structure,
                        rotor,
                        waves,
                        wind)
        

        # Monopile modal acceleration
        # FIXME: implement formula for alphaDotDot as a function of GF, GD, GK, GM
        # FIX: alphaDotDot = (GF-structure["GD"]*alphaDot-structure["GK"]*alpha)/(structure["GM"])
        alphaDotDot = (GF-structure["GD"]*alphaDot-structure["GK"]*alpha)/(structure["GM"])
        
        dqdtOut = np.zeros(2)
        dqdtOut[0] = alphaDot
        dqdtOut[1] = alphaDotDot

        return dqdtOut

def GFCalc(t, alphaDot, structure, rotor, waves, wind):    
    
    x_dot = alphaDot*structure["phiNodal"]
    
    # Wave contribution to Generalized Forcing
    wetPart = np.less_equal(structure["zBeamNodal"], 0.)    
    x_dot_submerged = x_dot[wetPart]
    phiNodalSubmerged = structure["phiNodal"][wetPart];
    
    i_ = lookup(waves["t"], t)
    u, ut, z = waves["u"][i_,:], waves["ut"][i_,:], waves["z"]
    df = forceDistributed(structure, u, ut, z, x_dot_submerged)
    
    # FIXME: Add the generalized forcing from the waves
    # FIX: GFWaves = np.trapz(df*phiNodalSubmerged, phiNodalSubmerged)
    GFWaves = np.trapz(df*phiNodalSubmerged, phiNodalSubmerged)
    
    # Wind contribution to Generalized Forcing
    i_ = lookup(wind["t"], t)
    V_10, V_hub = wind["V_10"], wind["V_hub"][i_]
    phiHub = structure["phiNodal"][-1]
    
    # FIXME: Compute the velocity of the hub due to
    # structural deformation.
    # FIX: x_dot_rotor = alphaDot*phiHub
    x_dot_rotor = alphaDot*phiHub
    
    # FIXME: add the generalized forcing from the wind
    # FIX: GFWind = F_aero*phiHub
    F_aero = F_wind(rotor, V_10, V_hub, x_dot_rotor)
    GFWind = F_aero*phiHub
    
    GF = GFWind + GFWaves
    
    return GF