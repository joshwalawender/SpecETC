from astropy import units as u
from .telescope import Telescope

SVX152 = Telescope('SVX152', 152, 1200)
EdgeHD8 = Telescope('EdgeHD 8', 203, 2032, obstruction=2.7*25.4)
EdgeHD9 = Telescope('EdgeHD 9.25', 235, 2350, obstruction=3.35*25.4)
EdgeHD11 = Telescope('EdgeHD 11', 280, 2800, obstruction=3.75*25.4)
EdgeHD14 = Telescope('EdgeHD 14', 356, 3910, obstruction=4.5*25.4)
ACF14 = Telescope('14ACF', 356*u.mm, 3556*u.mm, obstruction=6.5*25.4*u.mm)
ACF14reduced = Telescope('14ACFx0.63', 356*u.mm, 3556*u.mm*0.63, obstruction=6.5*25.4)
Newtonian13 = Telescope('13inch-f/3', 305, 1052, obstruction=4.0*25.4)
K1 = Telescope('K1', 10*1e3, 15*10*1e3)
