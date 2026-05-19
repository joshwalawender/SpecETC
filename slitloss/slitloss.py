import numpy as np
from astropy.modeling.models import Gaussian2D, Moffat2D
from astropy import units as u


def slit_throughput(seeing, tel, spec, det,
                    amplitude=1, alpha=1,
                    sample_size=10,
                   ):
    # Generate PSF
    if isinstance(seeing, u.Quantity): seeing = seeing.to(u.arcsec).value
    psf = Moffat2D(amplitude=amplitude, x_0=0, y_0=0,
                   gamma=seeing/2, alpha=alpha)
#     psf = Gaussian2D(amplitude=amplitude, x_mean=0, y_mean=0,
#                      x_stddev=seeing/2.355, y_stddev=seeing/2.355)
    pscale = tel.pixel_scale(det.pixel_size)
    gstart = -5*seeing # arcsec
    gend = -gstart+pscale.value/sample_size # arcsec
    gx = np.arange(gstart, gend, pscale.value/sample_size)
    gy = np.arange(gstart, gend, pscale.value/sample_size)
    xv, yv = np.meshgrid(gx, gy)
    total_psf_flux = np.sum(psf(xv, yv))

    # Generate Slit Mask
    slit_width = tel.slit_width(spec.slit_size).value
    slit_start = int(np.argmin(abs(gx+slit_width/2)))
    slit_end = int(np.argmin(abs(gx-slit_width/2)))
    w = (xv > gx[slit_start]) & (xv < gx[slit_end])
    wint = np.array(w, dtype=int)

    slit_flux = float(np.sum(psf(xv[w], yv[w])))
    slit_throughput = float(slit_flux/total_psf_flux)

    trace_profile = np.sum(psf(xv, yv)*wint, axis=1)

    return slit_throughput, trace_profile




#     if plot:
#         plt.subplot(2,Nt,i+1)
#         title = f"{tel.name}\n"
#         title += f"{seeing:.1f}'' Seeing, {slit_width:.2f}'' slit\n"
#         title += f"Eff. = {slit_throughput:.1%} "
#         title += f"Flux = {slit_throughput*tel.area:.0f}"
#         plt.title(title, size=10)
#         plt.imshow(psf(xv, yv)*wint, cmap='Greys', origin='lower', vmin=0, vmax=1e-16)
#         plt.gca().set_xticks([])
#         plt.gca().set_yticks([])
#         plt.axvline(slit_start, color='k')
#         plt.axvline(slit_end, color='k')
# 
#         plt.subplot(2,Nt,i+1+Nt)
#         plt.imshow(psf(xv, yv)*(1-wint), cmap='Reds', origin='lower', vmin=0, vmax=1e-16)
#         plt.gca().set_xticks([])
#         plt.gca().set_yticks([])
#         plt.axvline(slit_start, color='k')
#         plt.axvline(slit_end, color='k')
# if plot: plt.show()
