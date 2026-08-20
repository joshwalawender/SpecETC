from telescope import *
from spectrograph import *
from detector import *

from pathlib import Path
import random
import numpy as np
from astropy import units as u
from astropy import visualization as vis
from astropy.table import Table
from astropy.modeling.models import Gaussian2D, Moffat2D
import synphot

from matplotlib import pyplot as plt


def calculate_slit_throughput(seeing, tel, slit_size, det, alpha=1, sample_size=1):
    seeing = seeing.to(u.arcsec)
    psf = Moffat2D(amplitude=1, x_0=0, y_0=0,
                   gamma=seeing.value/2, alpha=alpha)
    pscale = tel.pixel_scale(det.pixel_size)
    gstart = -10*seeing.value # arcsec
    gend = -gstart+pscale.value/sample_size # arcsec
    gx = np.arange(gstart, gend, pscale.value/sample_size)
    gy = np.arange(gstart, gend, pscale.value/sample_size)
    xv, yv = np.meshgrid(gx, gy)
    total_psf_flux = np.sum(psf(xv, yv))
    
    # Generate Slit Mask
    slit_width = tel.slit_width(slit_size)
    slit_start = int(np.argmin(abs(gx+slit_width.value/2)))
    slit_end = int(np.argmin(abs(gx-slit_width.value/2)))
    w = (xv > gx[slit_start]) & (xv < gx[slit_end])
    wint = np.array(w, dtype=int)
    
    slit_flux = float(np.sum(psf(xv[w], yv[w])))
    slit_throughput = float(slit_flux/total_psf_flux)
    
    trace_profile = np.sum(psf(xv, yv)*wint, axis=1)
    trace_profile *= slit_throughput/np.sum(trace_profile)

    return slit_throughput, trace_profile


def simulate_1D_spectrum(star, sky, telescope, spectrograph, detector, plot=True):
    spectrograph.set_detector(detector)
    spectrograph.generate_binset()
    spectrograph.get_geometric_efficiency(telescope)
    slit_throughput, trace_profile = calculate_slit_throughput(sky.seeing,
                                                               telescope,
                                                               spectrograph.slit_size,
                                                               detector)
    total_efficiency = telescope.efficiency
    total_efficiency *= slit_throughput
    total_efficiency *= spectrograph.geometric_efficiency
    total_efficiency *= spectrograph.grating.efficiency
    total_efficiency *= detector.efficiency

# #     print(f"Telescope:    {telescope.efficiency:.0%}")
#     print(f"Slit:         {slit_throughput:.0%}")
#     print(f"Grating Beam: {spectrograph.geometric_efficiency:.0%}")
# #     print(f"Spectrograph: {spectrograph.grating_efficiency:.0%}")
# #     print(f"Detector:     {detector.efficiency:.0%}")
# #     print(f"       OTAL = {total_efficiency:.0%}")

    # Create synphot.Observation of target
    obs = synphot.Observation(star, total_efficiency,
                              binset=spectrograph.binset,
                              force='taper')
    signal = obs.sample_binned(spectrograph.binset)
    signal *= spectrograph.AperPix
    signal *= telescope.area
    signal *= detector.exptime

    # Create synphot.Observation of sky
    skyobs = synphot.Observation(sky.skyspec, total_efficiency,
                                 binset=spectrograph.binset,
                                 force='taper')
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
        plt.gca().set_xticks([])
        plt.gca().set_yticks([])
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        plt.subplot(3,1,(2,3))
        plt.plot(spectrograph.binset, extracted_spectrum, 'k-')
        plt.plot(spectrograph.binset, image.sum(axis=0), 'k-', alpha=0.2)

        lines_file = Path(__file__).parent.parent / 'data' / 'spectral_lines.csv'
        lines = Table.read(lines_file, format='ascii.csv')
        wavmin = min(spectrograph.binset)
        wavmax = max(spectrograph.binset)
        for line in lines:
            if line['wav'] > wavmin and line['wav'] < wavmax:
                color = 'b' if line['wav'] < 5500 else 'r'
                alpha = 0.3
                plt.axvline(line['wav'], ymin=0.2, ymax=0.8,
                            color=color, alpha=alpha)
                plt.text(line['wav'], 0.04*max(extracted_spectrum), line['name'],
                         color=color, alpha=alpha)

        plt.xlim(min(spectrograph.binset), max(spectrograph.binset))
        plt.ylim(0,1.1*max(extracted_spectrum))
        plt.ylabel('Extracted (phot/1Dpix)')
        plt.xlabel('Wavelength (A)')
        plt.grid()

        plt.show()

    return image, extracted_spectrum, extracted_variance




def SpecETC(star, sky, telescope, spectrograph, detector, plot=True):
    spectrograph.get_geometric_efficiency(telescope)
    signal, skysignal, trace_profile = simulate_1D_spectrum(
                        star, sky, telescope, spectrograph, detector, plot=False)
    image, spectrum, variance = simulate_2D_spectrum(signal, skysignal, trace_profile,
                        star, sky, telescope, spectrograph, detector, plot=plot)
    return image, spectrum, variance