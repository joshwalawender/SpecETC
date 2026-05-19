import numpy as np
from astropy import units as u
from synphot import SpectralElement
from synphot.models import Box1D


class Telescope(object):
    def __init__(self, name, aperture, fl, obstruction=0):
        self.name = name
        self.aperture = aperture
        if not isinstance(self.aperture, u.Quantity): self.aperture *= u.mm
        self.fl = fl
        if not isinstance(self.fl, u.Quantity): self.fl *= u.mm
        self.fratio = self.fl.to(u.mm)/self.aperture.to(u.mm)
        # ALPY 600 Spectrograph only accepts light from f/4 beam or slower
        if self.fratio < 4:
            self.aperture = self.fl/4
        self.obstruction = obstruction
        if not isinstance(self.obstruction, u.Quantity): self.obstruction *= u.mm
        self.area = np.pi*(self.aperture**2-self.obstruction**2).to(u.cm**2)
        self.efficiency = SpectralElement(Box1D, amplitude=0.90**2, x_0=5500, width=3800)


    def slit_width(self, slit_size):
        if not isinstance(slit_size, u.Quantity): slit_size *= u.micron
        return 206265*u.arcsec*slit_size.to(u.mm)/self.fl.to(u.mm)

    def pixel_scale(self, pixel_size):
        if not isinstance(pixel_size, u.Quantity): pixel_size *= u.micron
        return 206265*u.arcsec*pixel_size.to(u.mm)/self.fl.to(u.mm)

    def __str__(self):
        return f"{self.name}"

    def __repr__(self):
        return f"{self.name} (EffD={self.aperture:.0f}mm, FL={self.fl:.0f}mm)"