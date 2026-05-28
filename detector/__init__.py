from astropy import units as u
from .detector import Detector

IMX183 = Detector('IMX183', 2.40*u.micron, 5496, 3672, 'data/SonyIMX/QE.csv', 1.2)
IMX533 = Detector('IMX533', 3.76*u.micron, 3003, 3003, 'data/SonyIMX/QE.csv', 1.2)
LRISB = Detector('LRIS-B', 15*u.micron, 4096, 4096, 90, 4.0, wav1=3200, wav2=7000)