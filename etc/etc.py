from telescope import *
from spectrograph import *
from detector import *
from spectralmodels import get_sky_spectrum

from pathlib import Path
import random
import numpy as np
from astropy import units as u
from astropy import visualization as vis
import synphot

from matplotlib import pyplot as plt


def SpecETC(source_spectrum, exptime=60*u.second, seeing=2.5*u.arcsec,
            telescope=ACF14reduced,
            spectrograph=Alpy600,
            detector=IMX183,
            skymag=20.5,
            plot=True):
    spectrograph.generate_binset(detector)
    if not isinstance(seeing, u.Quantity): seeing *= u.arcsec
    seeing = seeing.to(u.arcsec)
    slit_throughput, trace_profile = spectrograph.slit_throughput(seeing, telescope, detector)
    total_efficiency = telescope.efficiency
    total_efficiency *= slit_throughput
    total_efficiency *= spectrograph.efficiency
    total_efficiency *= detector.efficiency

    # Create synphot.Observation of target
    obs = synphot.Observation(source_spectrum, total_efficiency, binset=spectrograph.binset)
    signal = obs.sample_binned(spectrograph.binset)
    signal *= spectrograph.AperPix
    signal *= telescope.area
    signal *= exptime

    # Create synphot.Observation of sky
    skyspec = get_sky_spectrum(mag=skymag)
    skyobs = synphot.Observation(skyspec, total_efficiency, binset=spectrograph.binset)
    skysignal = skyobs.sample_binned(spectrograph.binset)
    skysignal *= spectrograph.AperPix
    skysignal *= telescope.area
    skysignal *= exptime
    pscale = telescope.pixel_scale(detector.pixel_size)
    skysignal *= (pscale.value)**2

    # Simulate 2D spectrum
    P = trace_profile/trace_profile.sum()
    image = np.random.normal(0, detector.RN, (len(trace_profile),len(signal)))
    extracted_spectrum = np.zeros(len(signal))
    extracted_variance = np.zeros(len(signal))
    for specpix, specval in enumerate(signal.value):
        image[:,specpix] += P*specval
        image[:,specpix] += skysignal.value[specpix]
    for specpix, specval in enumerate(signal.value):
        # Optimal Extract 2D Spectrum
        V = P*specval + detector.RN**2
        extracted_spectrum[specpix] = np.sum(P*P*specval/V)/np.sum(P**2/V)
        extracted_variance[specpix] = 1/np.sum(P**2/V)

    SNR = extracted_spectrum/extracted_variance**0.5

    mean_signal = np.mean(extracted_spectrum)
    mean_SNR = np.mean(SNR)

    # Make Plots
    if plot:
        plt.figure(figsize=(10,5))
        tstr = (f'{telescope.name}, {spectrograph.name}, {detector.name}, '\
                f'seeing={seeing:.1f}, exptime={exptime:.0f}\n'\
                f'Mean Extracted Signal={mean_signal:.0f} ct/pix, Mean SNR={mean_SNR:.1f}')
    
        plt.subplot(5,1,1)
        plt.title(tstr)
        norm = vis.ImageNormalize(image, interval=vis.PercentileInterval(99.9),
                                  stretch=vis.LinearStretch())
        plt.imshow(image, norm=norm, cmap='Grays')
        plt.xlim(0.4*image.shape[1], 0.6*image.shape[1])
        plt.gca().set_xticks([])
        plt.gca().set_yticks([])
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])
    
        plt.subplot(5,1,(2,3))
        plt.plot(spectrograph.binset, extracted_spectrum, 'b-')
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
