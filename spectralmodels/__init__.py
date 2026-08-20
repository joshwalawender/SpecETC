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





# Atmospheric Extinction Curves: https://specreduce.readthedocs.io/en/latest/extinction.html

# from specreduce.calibration_data import AtmosphericExtinction
# import numpy as np
# from astropy import units as u
# from astropy.table import Table, Column

# # 1. Fetch the empirical Maunakea Observatory (MKO) atmospheric extinction model
# # This returns extinction coefficients in magnitudes per unit airmass.
# mtham_extinction = AtmosphericExtinction(model="mtham")

# # 2. Extract wavelengths and convert them to Angstroms for synphot compatibility
# waves_angstrom = mtham_extinction.wavelength.to(u.AA).value
# extinction_mags = mtham_extinction.extinction_mag.value

# # 3. Define your target observation airmass
# airmass = 1.5

# # 4. Calculate fractional transmission: T = 10**(-0.4 * extinction * airmass)
# transmission_values = 10**(extinction_mags*airmass/-2.5)

# # Write this to file for future use (like QE)
# data = {'wav': waves_angstrom, 'T': transmission_values}
# t = Table(data)
# t.write('./data/atmosphere/Lick_airmass1p5.csv', format='ascii.csv')