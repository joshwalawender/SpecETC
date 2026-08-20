from pathlib import Path
import numpy as np
from astropy import units as u
from synphot import SpectralElement
from synphot.models import Box1D


class Detector(object):
    def __init__(self, name, pixel_size, Nx, Ny, QE, RN, exptime=300*u.second,
                 wav1=None, wav2=None):
        self.name = name
        self.pixel_size = pixel_size
        self.pixel_shape = np.array((Nx, Ny))
        self.RN = RN
        self.size = self.pixel_size.to(u.mm)*self.pixel_shape
        self.exptime = exptime

        self.efficiency = None
        if type(QE) in [str, Path]:
            QEfile = Path(QE).expanduser().absolute()
            if QEfile.exists():
#                 print(f'{name}: Reading detector efficiency from {QEfile}')
                self.efficiency = SpectralElement.from_file(str(QEfile))
            else:
                print('ERROR: {QEfile} not found')
        else:
            assert wav1 is not None
            assert wav2 is not None
            wavc = (wav1+wav2)/2
            wav_width = wav2-wav1
            if type(QE) in [float, int]:
                self.efficiency = SpectralElement(Box1D, amplitude=QE,
                                                  x_0=wavc, width=wav_width)


    def __str__(self):
        return f"{self.name}"

    def __repr__(self):
        return f"{self.name}"
