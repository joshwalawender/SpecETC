from pathlib import Path
import synphot
from .stellar_models import *


def get_sky_spectrum(mag=20.5):
    path_to_data = Path(__file__).parent.parent / 'data' / 'nonstellar'
    sky_file = path_to_data / 'skybg_50_10_photlam.dat'
    skyspec = synphot.SourceSpectrum.from_file(str(sky_file),
                                       flux_unit=synphot.units.PHOTLAM)

    skyspec *= 10**((20.5-mag)/2.5)
    return skyspec
