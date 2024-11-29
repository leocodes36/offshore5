import numpy as np
from .common import pad2

def calculateKaimalSpectrum(windDict):
    
    # Store it inside the wind dictionary
    outputDict = dict()
    outputDict.update(windDict)
    
    # Calculate frequency information
    df = windDict["TDur"]**-1
    f = np.arange(df, windDict["fHighCut"], df)

    # Calculate the Kaimal spectrum
    # FIXME: Program spectrum 
    # FIX: Spectrum = (4*windDict["I"]**2*windDict["V_10"]*windDict["l"])/((1+6*f*windDict["l"]/windDict["V_10"])**(5/3))
    Spectrum = (4*windDict["I"]**2*windDict["V_10"]*windDict["l"])/((1+6*f*windDict["l"]/windDict["V_10"])**(5/3))
    amplitudeSpectrum = np.sqrt(2*Spectrum*df)

    outputDict["Spectrum"] = Spectrum
    outputDict["amplitudeSpectrum"] = amplitudeSpectrum
    outputDict["f"] = f
    
    return outputDict

def calculateWindTimeSeries(windDict):
    t = windDict["t"]
    f = windDict["f"]
    windTimeSeries = np.zeros_like(t)
    
    for i_, _ in enumerate(t):
        for j_, _ in enumerate(f):
            # FIXME: add random phases
            # FIX: in the cosine + windDict["randomPhases"][j_], randomPhases need to be added to windDict before calling this function!!
            windTimeSeries[i_] += windDict["amplitudeSpectrum"][j_]*np.cos(2*np.pi*f[j_]*t[i_]+ windDict["randomPhases"][j_])
    
    # Store the result
    outputDict = dict()
    outputDict.update(windDict)
    outputDict["t"] = t
    outputDict["V_hub"] = windTimeSeries + windDict["V_10"]
    
    return outputDict

def calculateWindTimeSeriesFFT(windDict):
    t = windDict["t"]
    f = windDict["f"]

    windTimeSeriesKernel = np.zeros_like(t)
    windTimeSeries = np.zeros_like(t)
    
    M = len(t)
    # FIXME: compute the fft kernel and perform the IFFT
    windTimeSeriesKernel = windDict["amplitudeSpectrum"]*np.exp(1j*windDict["randomPhases"])
    windTimeSeriesKernelPadded = pad2(windTimeSeriesKernel, M)
    windTimeSeries = M*np.real(np.fft.ifft(windTimeSeriesKernelPadded))
    
    # Store the result
    outputDict = dict()
    outputDict.update(windDict)
    outputDict["t"] = t
    outputDict["V_hub"] = windTimeSeries + windDict["V_10"]
    
    return outputDict
