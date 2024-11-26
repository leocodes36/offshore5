import numpy as np

def calculateRegularWaveFrequencyInformation(waveDict):
           
    H = waveDict["Hs"]
    T = waveDict["Tp"]

    
    isRegular = waveDict.get("regular", False)
    if not isRegular:
          raise ValueError("Your input dictionary specifies an irregular sea state, but you have called the regular wave routine.")
                       
    # FIXME: setup a single frequency, and a single amplitude
    # so that you can use the old routine calculateKinematics
    # also to calculate regular waves.
    # Important: f and a should be 1-element arrays, otherwise
    # the function calculateKinematics will fail.
    
    f = np.array([1/T])
    a = np.array([H/2])
    
    # Store it inside the wind dictionary
    outputDict = dict()
    outputDict.update(waveDict)
    outputDict["Spectrum"] = np.nan
    outputDict["amplitudeSpectrum"] = a
    outputDict["f"] = f
    outputDict["randomPhases"] = np.array([0.])
    
    return outputDict