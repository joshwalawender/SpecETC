from astropy import units as u
from .spectrograph import Spectrograph, Alpy600Spec, Grating, Optic


# StarExCol = Optic('Collimator', 75/3.5*u.mm, 75*u.mm)
# StarExGrating = Grating(72*u.deg, 2400/u.mm, 25*u.mm, m=1, minus=False)
# StarExCam = Optic('Camera', 75/3.5*u.mm, 75*u.mm)
# StarEx = Spectrograph('StarEx', 32*u.micron, StarExCol, StarExGrating, StarExCam, 0.25)
# StarEx600Grating = Grating(28.5*u.deg, 600/u.mm, 25*u.mm, m=1, minus=False)
# StarEx600 = Spectrograph('StarEx600', 32*u.micron, StarExCol, StarEx600Grating, StarExCam, 0.25)

# UVEX: alpha-beta = 27 deg
# It offers a useful image length (without degradation of quality) of 20 mm
UVEXCol = Optic('Collimator', 25*u.mm, 100*u.mm)
UVEXCam = Optic('Camera', 25*u.mm, 100*u.mm)
UVEX300 = Spectrograph('UVEX300', 23*u.micron, UVEXCol,
                       Grating(15.5*u.deg, 300/u.mm, 25*u.mm, m=1, minus=False),
                       UVEXCam, 0.35, cwav=6500*u.Angstrom,
                       a_to_b=27*u.deg, max_detector_width=20*u.mm)
UVEX600 = Spectrograph('UVEX600', 23*u.micron, UVEXCol,
                       Grating(15.5*u.deg, 600/u.mm, 25*u.mm, m=1, minus=False),
                       UVEXCam, 0.35, cwav=6500*u.Angstrom,
                       a_to_b=27*u.deg, max_detector_width=20*u.mm)
UVEX1200 = Spectrograph('UVEX1200', 23*u.micron, UVEXCol,
                       Grating(15.5*u.deg, 1200/u.mm, 25*u.mm, m=1, minus=False),
                       UVEXCam, 0.35, cwav=6500*u.Angstrom,
                       a_to_b=27*u.deg, max_detector_width=20*u.mm)

Alpy600 = Alpy600Spec('Alpy600', 23*u.micron, 480*u.A/u.mm, 0.3, wav1=3600, wav2=7400)
