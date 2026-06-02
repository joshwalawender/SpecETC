from telescope import ACF14, ACF14reduced
from spectrograph import Alpy600
from detector import IMX183

import random
import numpy as np
from astropy import units as u
from astropy import visualization as vis
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
    slit_throughput, trace_profile = spectrograph.slit_throughput(seeing, telescope, detector)
    total_efficiency = telescope.efficiency
    total_efficiency *= slit_throughput
    total_efficiency *= spectrograph.efficiency
    total_efficiency *= detector.efficiency
    obs = Observation(source_spectrum, total_efficiency, binset=spectrograph.binset)
    specdata = obs.sample_binned(spectrograph.binset)
    specdata *= spectrograph.AperPix
    specdata *= telescope.area
    signal = specdata*exptime
    noise = ( (specdata*exptime).value + np.ones(specdata.shape)*detector.RN**2)**0.5
    noisy_spectrum = [s+random.gauss(mu=0.0, sigma=noise[i])for i,s in enumerate(signal.value)]

    # Simulate 2D spectrum
    P = trace_profile/trace_profile.sum()
    image = np.random.normal(0, detector.RN, (len(trace_profile),len(signal)))
    extracted_spectrum = np.zeros(len(signal))
    extracted_variance = np.zeros(len(signal))
    for specpix, specval in enumerate(signal.value):
        image[:,specpix] += P*specval
        # Optimal Extract 2D Spectrum
        V = P*specval + detector.RN**2
        extracted_spectrum[specpix] = np.sum(P*P*specval/V)/np.sum(P**2/V)
        extracted_variance[specpix] = 1/np.sum(P**2/V)

    SNR = extracted_spectrum/extracted_variance**0.5

    mean_signal = np.mean(extracted_spectrum)
    mean_SNR = np.mean(SNR)

    # Make Plots
    plt.figure(figsize=(10,6))
    tstr = (f'Mean Extracted Signal={mean_signal:.0f} ct/pix, Mean SNR={mean_SNR:.0f}\n'\
            f'{telescope.name}, {spectrograph.name}, {detector.name}, '\
            f'seeing={seeing:.1f}, exptime={exptime:.0f}')

    plt.subplot(5,1,1)
    plt.title(tstr)
#     norm = vis.No
    plt.imshow(image, cmap='Grays')
    plt.gca().set_xticks([])
    plt.gca().set_yticks([])
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    plt.subplot(5,1,(2,3))
    plt.plot(spectrograph.binset, extracted_spectrum, 'b-')
#     plt.plot(spectrograph.binset, image.sum(axis=0), 'k-', alpha=0.3)
    plt.xlim(min(spectrograph.binset), max(spectrograph.binset))
    plt.ylim(0,1.1*max(extracted_spectrum))
    plt.ylabel('Extracted (phot/pix)')
    plt.gca().set_xticklabels([])
    plt.grid()

    plt.subplot(5,1,(4,5))
    plt.plot(spectrograph.binset, SNR, 'b-')
    plt.xlim(min(spectrograph.binset), max(spectrograph.binset))
    plt.ylim(0,1.1*max(SNR))
    plt.xlabel('Wavelength (A)')
    plt.ylabel('SNR')
    plt.grid()

    plt.show()

    return extracted_spectrum, extracted_variance, image
