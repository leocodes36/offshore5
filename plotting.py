import pylab as plt
import numpy as np
import scipy.signal as ss

def makeplots(wind, waves, structure, response, timeInfo, color, **kwargs):

    # Check if the axis has been passed otherwise initialize plot
    ax = kwargs.get("ax", np.zeros(2))
    if not ax.any():
        _, ax = plt.subplots(4,2, sharex='col')
    
    transcut = int(0.5 * len(wind["t"]) - 1)
    transcutr = int(0.5 * len(response["t"]) - 1)
    
    ax[0,0].plot(wind["t"], wind["V_hub"], color=color)
    ax[0,0].set_ylabel("V[m/s]")
    
    f, _, S = freqSpectrum(wind["t"][transcut:], wind["V_hub"][transcut:])
    ax[0,1].plot(f, S, color=color)
    ax[0,1].set_ylabel("PSD [m^2/s^2 / Hz]")
    
    ax[1,0].plot(waves["t"], waves["eta"], color=color)
    ax[1,0].set_ylabel("eta[m]")
    
    f, _, S = freqSpectrum(waves["t"][transcut:],  waves["eta"][transcut:])
    ax[1,1].plot(f, S, color=color)
    ax[1,1].set_ylabel("PSD [m^2 / Hz]")    

    ax[2,0].plot(response["t"], response["x1"], color=color)
    ax[2,0].set_ylabel("surge[m]")
    
    f, _, S = freqSpectrum(response["t"][transcutr:], response["x1"][transcutr:])
    ax[2,1].plot(f, S, color=color)
    ax[2,1].set_ylabel("PSD [m^2 / Hz]")    
    ax[2,1].axvline(structure["fnat"][0], color='k', linestyle='--')

    ax[3,0].plot(response["t"], np.rad2deg(response["x5"]), color=color)    
    ax[3,0].set_ylabel("pitch [deg]")
    ax[3,0].set_xlabel("Time [s]")
    
    f, _, S = freqSpectrum(response["t"][transcutr:], np.rad2deg(response["x5"][transcutr:]))
    ax[3,1].plot(f, S, color=color)
    ax[3,1].set_ylabel("PSD [deg^2 / Hz]")
    ax[3,1].set_xlabel("Frequency [Hz]")
    ax[3,1].axvline(structure["fnat"][1], color='k', linestyle='--')
        
    ax[3,1].set_xlim([0., timeInfo["fHighCut"]])
    [ax_.grid(True) for ax_ in ax.ravel()]
    
    plt.tight_layout()
    
    return ax

def freqSpectrum(t, x, flagMean=False):
    
    # Takes a time vector and a signal and returns the one-sided Fourier 
    # coefficients obtained with FFT, as well as the PSD function

    # flagmean=False: f[0]=df, a[0]=mean(x) is removed
    # flagmean=True: f[0]=0,  a[0]=mean(x) is kept
    
    if len(t) == 1:
        raise ValueError("The signal has only one element.")
    
    else:
        Tmax = t[-1] - t[0]
        df = Tmax**-1
        
        if flagMean:
            f = np.arange(len(t))*df
        else:
            f = (np.arange(1,len(t)+1)*df)
            
        a = np.fft.fft(x) / len(t)
        
        dt = t[1] - t[0]        
        f_nyq = 1./2./dt        
        a[f > f_nyq] = 0.
        
        if flagMean:
            a[1:] = 2*a[1:]
        else:
            a = 2*a[1:len(f)]
            a = np.append(a, 0.)
            
        # PSD function
        
        S = np.abs(a)**2 / (2*df)
        
        return f, a, S
        
        