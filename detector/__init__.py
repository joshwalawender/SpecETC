from astropy import units as u
from .detector import Detector

IMX183 = Detector('IMX183', 2.40*u.micron, 5496, 3672, 'data/SonyIMX/QE.csv', 1.2)
IMX533 = Detector('IMX533', 3.76*u.micron, 3003, 3003, 'data/SonyIMX/QE.csv', 1.2)
ICX825 = Detector('ICX825', 6.45*u.micron, 1392, 1040, 'data/SonyIMX/QE.csv', 5)
IMX571 = Detector('IMX571', 3.76*u.micron, 6248, 4176, 'data/SonyIMX/QE.csv', 1.2)
