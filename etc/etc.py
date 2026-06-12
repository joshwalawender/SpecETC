from telescope import *
from spectrograph import *
from detector import *

from pathlib import Path
import random
import numpy as np
from astropy import units as u
from astropy import visualization as vis
import synphot

from matplotlib import pyplot as plt

def simulate_1D_spectrum(star, sky, telescope, spectrograph, detector, plot=True):
    spectrograph.generate_binset(detector)
    slit_throughput, trace_profile = spectrograph.slit_throughput(sky.seeing, telescope, detector)
    total_efficiency = telescope.efficiency
    total_efficiency *= slit_throughput
    total_efficiency *= spectrograph.efficiency
    total_efficiency *= detector.efficiency

    # Create synphot.Observation of target
    obs = synphot.Observation(star, total_efficiency, binset=spectrograph.binset)
    signal = obs.sample_binned(spectrograph.binset)
    signal *= spectrograph.AperPix
    signal *= telescope.area
    signal *= detector.exptime

    # Create synphot.Observation of sky
    skyobs = synphot.Observation(sky.skyspec, total_efficiency, binset=spectrograph.binset)
    skysignal = skyobs.sample_binned(spectrograph.binset)
    skysignal *= spectrograph.AperPix
    skysignal *= telescope.area
    skysignal *= detector.exptime
    pscale = telescope.pixel_scale(detector.pixel_size)
    skysignal *= (pscale.value)**2

    if plot:
        plt.figure(figsize=(10,4))
        plt.subplot(2,1,1)
        plt.title('Simulated Spectrum')
        plt.plot(spectrograph.binset, signal.value, 'b-')
        plt.xlim(min(spectrograph.binset), max(spectrograph.binset))
        plt.ylim(0,1.1*max(signal.value))
        plt.ylabel('Signal (phot/1Dpix)')
        plt.grid()
        plt.subplot(2,1,2)
        plt.plot(spectrograph.binset, skysignal.value, 'b-')
        plt.xlim(min(spectrograph.binset), max(spectrograph.binset))
        plt.ylim(0,1.1*max(skysignal.value))
        plt.ylabel('Sky Signal (phot/2Dpix)')
        plt.grid()

        plt.show()

    return signal, skysignal, trace_profile


def simulate_2D_spectrum(signal, skysignal, trace_profile,
                         star, sky, telescope, spectrograph, detector, plot=True):
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
        plt.figure(figsize=(10,3))
        tstr = (f'{telescope.name}, {spectrograph.name}, {detector.name}, '\
                f'seeing={sky.seeing:.1f}, exptime={detector.exptime:.0f}\n'\
                f'SNR: min={min(SNR):.1f}, mean={mean_SNR:.1f}, max={max(SNR):.1f}')

        plt.subplot(3,1,1)
        plt.title(tstr)
        norm = vis.ImageNormalize(image, interval=vis.PercentileInterval(99.9),
                                  stretch=vis.LinearStretch())
        plt.imshow(image, norm=norm, cmap='Grays')
#         plt.xlim(0.4*image.shape[1], 0.6*image.shape[1])
#         imwav1 = 5500
#         imwav2 = 6800
#         plt.xlim(min(np.where(Alpy600.binset > imwav1)[0]),
#                  max(np.where(Alpy600.binset < imwav2)[0]))
        plt.gca().set_xticks([])
        plt.gca().set_yticks([])
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])
    
        plt.subplot(3,1,(2,3))
        plt.plot(spectrograph.binset, extracted_spectrum, 'b-')
        plt.plot(spectrograph.binset, image.sum(axis=0), 'k-', alpha=0.2)
        plt.xlim(min(spectrograph.binset), max(spectrograph.binset))
        plt.ylim(0,1.1*max(extracted_spectrum))
        plt.ylabel('Extracted (phot/1Dpix)')
        plt.gca().set_xticklabels([])
        plt.grid()
    
#         plt.subplot(5,1,(4,5))
#         plt.plot(spectrograph.binset, SNR, 'b-')
#         plt.xlim(min(spectrograph.binset), max(spectrograph.binset))
#         plt.ylim(0,1.1*max(SNR))
#         plt.xlabel('Wavelength (A)')
#         plt.ylabel('SNR')
#         plt.grid()
    
        plt.show()

    return image, extracted_spectrum, extracted_variance




def SpecETC(star, sky, telescope, spectrograph, detector, plot=True):
    signal, skysignal, trace_profile = simulate_1D_spectrum(
                        star, sky, telescope, spectrograph, detector, plot=False)
    image, spectrum, variance = simulate_2D_spectrum(signal, skysignal, trace_profile,
                        star, sky, telescope, spectrograph, detector, plot=plot)
    return image, spectrum, variance