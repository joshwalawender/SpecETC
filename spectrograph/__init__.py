from astropy import units as u
from .spectrograph import Spectrograph

Alpy600 = Spectrograph('Alpy600', 23*u.micron, 480*u.A/u.mm, 0.3, wav1=3600, wav2=7400)
LRISBlue600 = Spectrograph('LRIS-B-600', 725*u.micron, 1.09/15*u.A/u.micron,
                           'data/LRIS-Blue/600_mirr_eff.dat',
                           magnification=1/6.5)
