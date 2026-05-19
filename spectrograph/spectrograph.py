import numpy as np
from astropy import units as u
from astropy.modeling.models import Gaussian2D, Moffat2D
from synphot import SpectralElement
from synphot.models import Box1D


class Spectrograph(object):
    def __init__(self, name, slit_size, dispersion, eff,
                 wav1=3600, wav2=7400):
        self.wav1 = wav1
        self.wav2 = wav2
        self.wavc = (wav1+wav2)/2
        self.wav_width = wav2-wav1
        self.name = name
        self.slit_size = slit_size
        self.dispersion = dispersion
        self.binset = None
        self.AperPix = None
        if isinstance(eff, float):
            self.grating_efficiency = SpectralElement(Box1D, amplitude=eff,
                                              x_0=self.wavc, width=self.wav_width)
        elif type(eff) in [str, Path]:
            efffile = Path(eff).expanduser().absolute()
            if efffile.exists():
                self.grating_efficiency = SpectralElement.from_file(str(efffile))
            else:
                self.grating_efficiency = SpectralElement(Box1D, amplitude=0.25,
                                                  x_0=self.wavc, width=self.wav_width)
        else:
            self.grating_efficiency = SpectralElement(Box1D, amplitude=0.25,
                                              x_0=self.wavc, width=self.wav_width)

    def __str__(self):
        return f"{self.name}"

    def __repr__(self):
        return f"{self.name}"

    def slit_throughput(self, seeing, tel, det, alpha=1, sample_size=1):
        # Generate PSF
        if isinstance(seeing, u.Quantity): seeing = seeing.to(u.arcsec).value
        psf = Moffat2D(amplitude=1, x_0=0, y_0=0,
                       gamma=seeing/2, alpha=alpha)
#         psf = Gaussian2D(amplitude=1, x_mean=0, y_mean=0,
#                          x_stddev=seeing/2.355, y_stddev=seeing/2.355)
        pscale = tel.pixel_scale(det.pixel_size)
        gstart = -5*seeing # arcsec
        gend = -gstart+pscale.value/sample_size # arcsec
        gx = np.arange(gstart, gend, pscale.value/sample_size)
        gy = np.arange(gstart, gend, pscale.value/sample_size)
        xv, yv = np.meshgrid(gx, gy)
        total_psf_flux = np.sum(psf(xv, yv))
    
        # Generate Slit Mask
        slit_width = tel.slit_width(self.slit_size).value
        slit_start = int(np.argmin(abs(gx+slit_width/2)))
        slit_end = int(np.argmin(abs(gx-slit_width/2)))
        w = (xv > gx[slit_start]) & (xv < gx[slit_end])
        wint = np.array(w, dtype=int)
    
        slit_flux = float(np.sum(psf(xv[w], yv[w])))
        slit_throughput = float(slit_flux/total_psf_flux)
    
        trace_profile = np.sum(psf(xv, yv)*wint, axis=1)
        trace_profile *= slit_throughput/np.sum(trace_profile)
        
        self.efficiency = self.grating_efficiency * slit_throughput
    
        return slit_throughput, trace_profile

    def generate_binset(self, det):
        self.AperPix = self.dispersion*det.pixel_size.to(u.mm)
        self.binset = np.arange(self.wav1, self.wav2, self.AperPix.value)
        