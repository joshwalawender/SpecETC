from telescope import ACF14, ACF14reduced
from spectrograph import Alpy600
from detector import IMX183

import random
import numpy as np
from astropy import units as u
from synphot import SpectralElement, Observation, SourceSpectrum
from synphot.models import BlackBody1D

from matplotlib import pyplot as plt


def SpecETC(source_spectrum, exptime=60*u.second, seeing=2.5*u.arcsec,
        telescope=ACF14reduced,
        spectrograph=Alpy600,
        detector=IMX183):
    spectrograph.generate_binset(detector)
    if not isinstance(seeing, u.Quantity): seeing *= u.arcsec
    seeing = seeing.to(u.arcsec)
    total_efficiency = telescope.efficiency
    total_efficiency *= spectrograph.slit_throughput(seeing, telescope, detector)[0]
    total_efficiency *= spectrograph.efficiency
    total_efficiency *= detector.efficiency
    obs = Observation(source_spectrum, total_efficiency, binset=spectrograph.binset)
    specdata = obs.sample_binned(spectrograph.binset)
    specdata *= spectrograph.AperPix
    specdata *= telescope.area
    signal = (specdata*exptime).value
    noise = ( (specdata*exptime).value + np.ones(specdata.shape)*detector.RN**2)**0.5
    noisy_spectrum = [s+random.gauss(mu=0.0, sigma=noise[i])for i,s in enumerate(signal)]

    mean_signal = np.mean(signal)
    mean_SNR = np.mean(signal/noise)

    plt.figure(figsize=(10,4))

    tstr = (f'Mean Signal={mean_signal:.0f} ct/pix, Mean SNR={mean_SNR:.0f}\n'\
            f'{telescope.name}, {spectrograph.name}, {detector.name}, '\
            f'seeing={seeing:.1f}, exptime={exptime:.0f}')

#     plt.subplot(2,1,1)
    plt.title(tstr)
    plt.plot(spectrograph.binset, signal, 'b-')
    plt.plot(spectrograph.binset, noisy_spectrum, 'k-', alpha=0.3)
    plt.xlabel('Wavelength (A)')
    plt.ylabel('photons/pixel')
    plt.grid()
#     plt.gca().set_xticklabels([])

#     plt.subplot(2,1,2)
#     plt.plot(spectrograph.binset, signal/noise, 'b-')
#     plt.xlabel('Wavelength (A)')
#     plt.ylabel('SNR/pix')
#     plt.grid()

    plt.show()

    return float(mean_signal), float(mean_SNR)
