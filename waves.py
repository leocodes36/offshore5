import numpy as np
from .common import dispersion
from .common import pad2

def calculateJONSWAPSpectrum(waveDict):
           
    Hs = waveDict["Hs"]
    Tp = waveDict["Tp"]
    
    # If defined, get the gamma. Ohterwise default to 1.0
    gamma = waveDict.get("gamma", 1.0)
                    
    # Calculate frequency information
    df = waveDict["TDur"]**-1.0
    f = np.arange(df, waveDict["fHighCut"], df)
    
    # Spectral width parameter
    sigma = np.ones(len(f))

    fp = 1./Tp

    sigma[f>fp] = 0.09
    sigma[f<=fp] = 0.07        

    # Calculate the Kaimal spectrum
    # Pierson-Moskowitz spectrum
    Spm = 5/16 * Hs**2 * fp**4 * f**(-5) * \
            np.exp( - 5/4 * (f  / fp)**(-4) ) 
    # Jonswap spectrum
    # FIXME: progarm the correct spectrum, 
    # FIXED
    Spectrum = np.ones_like(f)
    alpha = np.exp(-0.5 * ((f/fp - 1)/(sigma))**2)
    Spectrum = Spm * (1 - 0.287 * np.log(waveDict["gamma"])) * waveDict["gamma"]**alpha
    amplitudeSpectrum = np.sqrt(2*Spectrum*df)
    
    # Store it inside the wind dictionary
    outputDict = dict()
    outputDict.update(waveDict)
    outputDict["Spectrum"] = Spectrum
    outputDict["amplitudeSpectrum"] = amplitudeSpectrum
    outputDict["f"] = f
    
    return outputDict
       
def calculateFreeSurfaceElevationTimeSeries(waveDict):
    
    t = waveDict["t"]
    f = waveDict["f"]
    freeSurfTimeSeries = np.zeros_like(t)
    
    for i_, _ in enumerate(t):
        for j_, _ in enumerate(f):
            freeSurfTimeSeries[i_] += waveDict["amplitudeSpectrum"][j_]*np.cos(2*np.pi*f[j_]*t[i_] + waveDict["randomPhases"][j_])
    
    # Store the result
    
    outputDict = dict()
    outputDict.update(waveDict)    
    outputDict["t"] = t
    outputDict["eta"] = freeSurfTimeSeries
        
    return outputDict

def calculateKinematics(inputDict):
    
    t = inputDict["t"]
    f = inputDict["f"]
    omega = 2*np.pi*f
    
    h = inputDict["h"]
    z = inputDict["z"]
    u = np.zeros((len(t), len(z)))
    ut = np.zeros((len(t), len(z)))
    
    k = dispersion(f, inputDict["h"])
    
    for i_, _ in enumerate(t):
        u[i_, :] += np.sum(inputDict["amplitudeSpectrum"][:, None]*omega[:,None]*np.cosh(k[:,None]*(z[None,:]+h))/np.sinh(k[:,None]*h)*np.cos(omega[:,None]*t[i_] + inputDict["randomPhases"][:,None]), axis=0)
        
        ut[i_, :] += -np.sum(inputDict["amplitudeSpectrum"][:, None]*omega[:,None]**2*np.cosh(k[:,None]*(z[None,:]+h))/np.sinh(k[:,None]*h)*np.sin(omega[:,None]*t[i_] + inputDict["randomPhases"][:,None]), axis=0)
    
    outputDict = dict()
    outputDict.update(inputDict)        
    outputDict["u"] = u
    outputDict["ut"] = ut
    
    return outputDict

def calculateFreeSurfaceElevationTimeSeriesFFT(waveDict):
    
    t = waveDict["t"]
    f = waveDict["f"]
    
    # FIXME: compute the fft kernel and perform the IFFT
    M = len(t)
    # FIXED
    freeSurfTimeSeriesKernel = waveDict["amplitudeSpectrum"]*np.exp(1j*waveDict["randomPhases"])
    freeSurfTimeSeriesKernelPadded = pad2(freeSurfTimeSeriesKernel, M)
    freeSurfTimeSeries = M*np.real(np.fft.ifft(freeSurfTimeSeriesKernelPadded))
    
    # Store the result
    outputDict = dict()
    outputDict.update(waveDict)    
    outputDict["t"] = t
    outputDict["eta"] = freeSurfTimeSeries
        
    return outputDict

def calculateKinematicsFFT(inputDict):
    
    t = inputDict["t"]
    f = inputDict["f"]
    omega = 2*np.pi*f
    
    h = inputDict["h"]
    z = inputDict["z"]
    u = np.zeros((len(t), len(z)))
    ut = np.zeros((len(t), len(z)))
    
    k = dispersion(f, inputDict["h"])
    M = len(t)
    
    for j_, z_ in enumerate(z):
        
        # FIXME: compute the fft kernel and perform the IFFT
        # FIXED
        uKernel = omega*(np.cosh(k*(z_+h)))/(np.sinh(k*h))*inputDict["amplitudeSpectrum"]*np.exp(1j*inputDict["randomPhases"]) # compute the freq. domain kernel and pad to M
        u[:, j_] = M*np.real(np.fft.ifft(pad2(uKernel, M))) # perform the IFFT in this line
        
        utKernel = omega**2*(np.cosh(k*(z_+h)))/(np.sinh(k*h))*inputDict["amplitudeSpectrum"]*np.exp(1j*inputDict["randomPhases"]) # compute the freq. domain kernel and pad to M
        ut[:, j_] = M*np.real(np.fft.ifft(pad2(utKernel, M))) # perform the IFFT in this line
    
    outputDict = dict()
    outputDict.update(inputDict)        
    outputDict["u"] = u
    outputDict["ut"] = ut
    
    return outputDict
