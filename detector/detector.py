from pathlib import Path
import numpy as np
from astropy import units as u
from synphot import SpectralElement
from synphot.models import Box1D


class Detector(object):
    def __init__(self, name, pixel_size, Nx, Ny, QE, RN):
        self.name = name
        self.pixel_size = pixel_size
        self.pixel_shape = np.array((Nx, Ny))
        self.RN = RN
        self.size = self.pixel_size.to(u.mm)*self.pixel_shape
        if isinstance(QE, float):
            self.efficiency = SpectralElement(Box1D, amplitude=QE,
                                              x_0=6700, width=6000)
        elif type(QE) in [str, Path]:
            QEfile = Path(QE).expanduser().absolute()
            if QEfile.exists():
                self.efficiency = SpectralElement.from_file(str(QEfile))
            else:
                self.efficiency = SpectralElement(Box1D, amplitude=0.5,
                                                  x_0=6700, width=6000)
        else:
            self.efficiency = SpectralElement(Box1D, amplitude=0.5,
                                              x_0=6700, width=6000)


    def __str__(self):
        return f"{self.name}"

    def __repr__(self):
        return f"{self.name}"
