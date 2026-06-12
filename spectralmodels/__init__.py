from pathlib import Path
from astropy import units as u
import synphot
from .stellar_models import *


class Sky(object):
    def __init__(self, seeing, mag):
        if not isinstance(seeing, u.Quantity): seeing *= u.arcsec
        self.seeing = seeing.to(u.arcsec)
        self.mag = mag
        path_to_data = Path(__file__).parent.parent / 'data' / 'nonstellar'
        sky_file = path_to_data / 'skybg_50_10_photlam.dat'
        self.skyspec = synphot.SourceSpectrum.from_file(str(sky_file),
                                            flux_unit=synphot.units.PHOTLAM)
    
        self.skyspec *= 10**((20.5-mag)/2.5)
