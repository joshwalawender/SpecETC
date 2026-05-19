from astropy import units as u
from .spectrograph import Spectrograph

Alpy600 = Spectrograph('Alpy600', 23*u.micron, 480*u.A/u.mm, 0.3)
