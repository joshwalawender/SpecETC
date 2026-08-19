from pathlib import Path
import numpy as np
from astropy import units as u
from synphot import SpectralElement
from synphot.models import Box1D

from matplotlib import pyplot as plt


class Optic(object):
    def __init__(self, name, aperture, fl):
        self.name = name
        self.aperture = aperture.to(u.mm)
        self.fl = fl.to(u.mm)


class Grating(object):
    def __init__(self, alpha, glpmm, diameter,
                 m=1, minus=False, efficiency=0.25):
        self.alpha = alpha.to(u.deg)
        self.glpmm = glpmm
        self.diameter = diameter.to(u.mm)
        self.m = m
        self.minus = minus
        # Grating Efficiency
        if type(efficiency) in [str, Path]:
            efffile = Path(efficiency).expanduser().absolute()
            if efffile.exists():
                print(f'{name}: Reading spectrograph efficiency from {efffile}')
                self.efficiency = SpectralElement.from_file(str(efffile))
            else:
                print('ERROR: {efffile} not found')
        else:
            wav1 = 3000
            wav2 = 10000
            wavc = (wav1+wav2)/2
            wav_width = wav2-wav1
            self.efficiency = SpectralElement(Box1D, amplitude=efficiency,
                                              x_0=wavc, width=wav_width)

    def beta(self, wav):
        '''
        m*wav = d * (sin(alpha) +/- sin(beta))
        m*wav/d = sin(a) +/- sin(b)
        sin(b) = +/- (m*wav/d - sin(a))
        '''
        wav = wav.to(u.mm)
        sinb = self.m*wav*self.glpmm - np.sin(self.alpha.to(u.radian).value)
        if self.minus: sinb *= -1
        b = np.arcsin(sinb).to(u.deg)
        return b

    def wav(self, beta):
        '''
        wav = d*(sin(alpha)+/-sin(beta)) / m
        '''
        d = 1/self.glpmm.to(1/u.Angstrom)
        sign = 1 if self.minus is False else -1
        wav = d*(np.sin(self.alpha.to(u.radian).value) + sign*np.sin(beta.to(u.radian).value))
        return wav.to(u.angstrom)


class Spectrograph(object):
    def __init__(self, name, slit_size, collimator, grating, camera, eff, m=1,
                 cwav=6563*u.Angstrom, a_to_b=None, max_detector_width=None):
        self.name = name
        self.slit_size = slit_size.to(u.micron)
        self.collimator = collimator
        self.grating = grating
        self.camera = camera
        self.cwav = cwav
        self.a_to_b = a_to_b
        self.max_detector_width = max_detector_width
        if a_to_b is not None:
            self.find_alpha()
        self.slit_image = self.slit_size*self.camera.fl/self.collimator.fl
        self.lines = [3968.5, 3933.7, 4063, 4132, 4068, 4076, 5889, 5895,
                      6300, 6363, 6563, 6707, 6717, 6731, 8498, 8542, 8662]
        self.line_names = ['CaH', 'CaK', 'FeI', 'FeI', 'SII', 'SII',
                           'NaI', 'NaI', 'OI', 'OI', 'HII', 'LiI',
                           'SII', 'SII', 'CaII', 'CaII', 'CaII']

    def set_detector(self, det):
        self.det = det
        self.generate_binset()

    def find_alpha(self):
        '''
        sin(b) = +/- (m*wav/d - sin(a))
        sin(a) = +/- (m*wav/d - sin(b))
        '''
        alphas = np.arange(0,85,0.2)
        errs = []
        for alpha in alphas:
            self.grating.alpha = alpha*u.deg
            a_to_b_for_cwav = self.grating.alpha.value - self.grating.beta(self.cwav).to(u.deg).value
            errs.append(abs(a_to_b_for_cwav - self.a_to_b.to(u.deg).value))
        wmin = np.nanargmin(errs)
        self.grating.alpha = alphas[wmin]*u.deg
        if min(errs) > 0.3 or self.grating.alpha.value <= 0 or self.grating.alpha.value >= 85:
            print(f'Issues Finding Alpha for {self.name}')
            print(wmin, alphas[wmin]*u.deg, errs[wmin])
            print([float(e) for e in errs])
        return self.grating.alpha

    def get_geometric_efficiency(self, tel, plot=False):
        # Collimator
        input_beam_size = self.collimator.fl/tel.fratio
        collimator_geometric_eff = min([1, (self.collimator.aperture/input_beam_size)**2])
        if input_beam_size > self.collimator.aperture:
            input_beam_size = self.collimator.aperture
        # Grating
        alpha_rad = self.grating.alpha.to(u.radian).value
        beam_area = input_beam_size**2/4*np.pi
        res = 0.01*u.mm
        grating_projection = self.grating.diameter.to(u.mm)*np.cos(alpha_rad)
        xs = np.arange(0, self.grating.diameter.to(u.mm).value*np.cos(alpha_rad)/2, res.value)
        ys = np.arange(0, self.grating.diameter.to(u.mm).value/2, res.value)
        xvals, yvals = np.meshgrid(xs, ys)
        radius = (xvals**2 + yvals**2)**0.5
        grating_footprint = np.array(radius < input_beam_size.to(u.mm).value/2, dtype=int)
        area = grating_footprint.sum() * res**2 * 4
        self.geometric_efficiency = area/beam_area
        if plot:
            import matplotlib.pyplot as plt
            tlt = [f"Beam = {input_beam_size:.1f}",
                   f"grating = {grating_projection:.1f}",
                   f"efficiency = {self.geometric_efficiency:.0%}"]
            plt.figure(figsize=(6,6))
            plt.title('\n'.join(tlt))
            plt.imshow(grating_footprint)
            plt.show()

    def generate_binset(self):
        beta0 = self.grating.beta(self.cwav)
        cam_ps = 206.265*u.arcsec*self.det.pixel_size.to(u.micron).value/self.camera.fl.to(u.mm).value
        cam_mmscale = 206265*u.arcsec/self.camera.fl.to(u.mm)
        if self.max_detector_width is None:
            beta1 = (beta0-cam_ps*self.det.pixel_shape[0]/2).to(u.deg).value
            beta2 = (beta0+cam_ps*self.det.pixel_shape[0]/2).to(u.deg).value
        elif self.max_detector_width > self.det.size[0]:
            beta1 = (beta0-cam_ps*self.det.pixel_shape[0]/2).to(u.deg).value
            beta2 = (beta0+cam_ps*self.det.pixel_shape[0]/2).to(u.deg).value
        else:
            beta1 = (beta0-cam_mmscale*self.max_detector_width/2).to(u.deg).value
            beta2 = (beta0+cam_mmscale*self.max_detector_width/2).to(u.deg).value
        delta_beta = cam_ps.to(u.deg).value
        betas = np.arange(beta1, beta2, delta_beta)
        self.binset = np.array([self.grating.wav(b*u.degree).to(u.Angstrom).value for b in betas])
        wav1 = self.grating.wav(beta0).to(u.Angstrom)
        wav2 = self.grating.wav(beta0+cam_ps).to(u.Angstrom)
        self.AperPix = wav2-wav1
        self.slit_pix = self.slit_image/self.det.pixel_size.to(u.micron)
        self.slit_span = self.slit_pix*self.AperPix
        self.R = self.binset / self.slit_span.value
        return self.binset

    def plot_coverage(self):
        plt.figure(figsize=(10,3))
        wavmin = min(self.binset)
        wavmax = max(self.binset)
        meanR = np.mean(self.R)
        cwav = self.cwav.to(u.Angstrom).value
        t = (f"{self.name} + {self.det.name}: Cwav = {cwav:.0f} A\n"
             f"{wavmin:.0f} A - {wavmax:.0f} A ({wavmax-wavmin:.0f} A span)\n"
#              f"Slit = {self.slit_pix:.0f} pix = {self.slit_span:.0f} thus "
             f"R ~ {meanR:.0f} ({min(self.R):.0f} - {max(self.R):.0f}) ~ {3e5/meanR:.0f} km/s")
        plt.title(t)
        plt.plot(self.binset, self.R, 'k-')
        for i,line in enumerate(self.lines):
            if line > wavmin and line < wavmax:
                color = 'b' if line < 5500 else 'r'
                alpha = 0.5
                plt.axvline(line, ymax=0.5, color=color, alpha=alpha)
                plt.text(line, max(self.R)*0.98, self.line_names[i], color=color, alpha=alpha)
        plt.xlabel('Wavelength (A)')
        plt.ylabel('R')
        plt.xlim(wavmin, wavmax)
        plt.grid()

        plt.show()

    def __str__(self):
        return f"{self.name}"

    def __repr__(self):
        return f"{self.name}"



class Alpy600Spec(Spectrograph):
    def __init__(self, name, slit_size, dispersion, eff,
                 wav1=None, wav2=None, magnification=1):
        self.name = name
        self.slit_size = slit_size
        self.dispersion = dispersion
        self.binset = None
        self.AperPix = None
        self.magnification = magnification
        self.grating = Grating(0*u.deg, 600, m=1, diameter=25*u.mm, efficiency=eff)

    def get_geometric_efficiency(self, tel):
        if tel.fratio >= 4:
            self.geometric_efficiency = 1
        else:
            self.geometric_efficiency = (tel.fratio/4)**2

    def __str__(self):
        return f"{self.name}"

    def __repr__(self):
        return f"{self.name}"

    def generate_binset(self, det):
        self.AperPix = self.dispersion*det.pixel_size.to(u.mm)
        self.binset = np.arange(self.grating.efficiency.waveset[0].value,
                                self.grating.efficiency.waveset[-1].value,
                                self.AperPix.value)
        