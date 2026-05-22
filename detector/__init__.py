from astropy import units as u
from .detector import Detector

IMX183 = Detector('IMX183', 2.40*u.micron, 5496, 3672, 60, 1.2)
IMX533 = Detector('IMX533', 3.76*u.micron, 3003, 3003, 58, 1.2)
