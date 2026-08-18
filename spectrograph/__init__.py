from astropy import units as u
from .spectrograph import Spectrograph, Alpy600Spec, Grating, Optic


StarExCol = Optic('Collimator', 75/3.5*u.mm, 75*u.mm)
StarExGrating = Grating(72*u.deg, 2400/u.mm, 25*u.mm, m=1, minus=False)
StarExCam = Optic('Camera', 75/3.5*u.mm, 75*u.mm)
StarEx = Spectrograph('StarEx', 32*u.micron, StarExCol, StarExGrating, StarExCam, 0.25)
StarEx600Grating = Grating(28.5*u.deg, 600/u.mm, 25*u.mm, m=1, minus=False)
StarEx600 = Spectrograph('StarEx600', 32*u.micron, StarExCol, StarEx600Grating, StarExCam, 0.25)

# UVEX: alpha-beta = 27 deg
UVEXCol = Optic('Collimator', 25*u.mm, 100*u.mm)
UVEXCam = Optic('Camera', 25*u.mm, 100*u.mm)
UVEX150 = Spectrograph('UVEX150', 23*u.micron, UVEXCol,
                       Grating(15.5*u.deg, 150/u.mm, 25*u.mm, m=1, minus=False),
                       UVEXCam, 0.35, cwav=6500*u.Angstrom, a_to_b=27*u.deg)
UVEX300 = Spectrograph('UVEX300', 23*u.micron, UVEXCol,
                       Grating(15.5*u.deg, 300/u.mm, 25*u.mm, m=1, minus=False),
                       UVEXCam, 0.35, cwav=5400*u.Angstrom, a_to_b=27*u.deg)
UVEX300R = Spectrograph('UVEX300R', 23*u.micron, UVEXCol,
                       Grating(15.5*u.deg, 300/u.mm, 25*u.mm, m=1, minus=False),
                       UVEXCam, 0.35, cwav=7200*u.Angstrom, a_to_b=27*u.deg)
UVEX600 = Spectrograph('UVEX600', 23*u.micron, UVEXCol,
                       Grating(15.5*u.deg, 600/u.mm, 25*u.mm, m=1, minus=False),
                       UVEXCam, 0.35, cwav=5500*u.Angstrom, a_to_b=27*u.deg)
UVEX1200 = Spectrograph('UVEX1200', 23*u.micron, UVEXCol,
                       Grating(15.5*u.deg, 1200/u.mm, 25*u.mm, m=1, minus=False),
                       UVEXCam, 0.35, cwav=6500*u.Angstrom, a_to_b=27*u.deg)
UVEX1800 = Spectrograph('UVEX1800', 23*u.micron, UVEXCol,
                       Grating(15.5*u.deg, 1800/u.mm, 25*u.mm, m=1, minus=False),
                       UVEXCam, 0.35, cwav=6550*u.Angstrom, a_to_b=27*u.deg)

Alpy600 = Alpy600Spec('Alpy600', 23*u.micron, 480*u.A/u.mm, 0.3, wav1=3600, wav2=7400)
