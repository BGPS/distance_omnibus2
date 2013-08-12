from scipy.ndimage.morphology import binary_dilation
import numpy as np
from astropy.io import fits
from fits_utils import flatten_header, hcongrid

def spec_extract(spec_cube, mask, npix=1, weight=True):
    assert spec_cube.ndim == 3
    assert spec_cube.shape[1:] == mask.shape

    border_mask = binary_dilation(mask, iterations=npix) - mask > 0

    sl = (slice(None),) + np.nonzero(border_mask)
    border_spectrum = spec_cube[sl].squeeze().mean(axis=1)

    sl = (slice(None),) + np.nonzero(mask)

    if weight:
        wts = mask[np.nonzero(mask)]
    else:
        wts = mask.astype('bool')

    spectra = spec_cube[sl].squeeze()
    spectrum = (spectra*wts[np.newaxis,:]).sum(axis=1)/wts.sum()

    return spectrum,border_spectrum

def bgps_to_mask(imagefile, labelmask, target, objectnumber=1):
    imfile = fits.open(imagefile)
    imdata = imfile[0].data
    inhdr = imfile[0].header
    maskdata = fits.getdata(labelmask)

    outhdr = flatten_header(fits.getheader(target))

    mask = (maskdata==objectnumber).astype('float')
    outdata = hcongrid(imdata, inhdr, outhdr)
    outmask = hcongrid(mask, inhdr, outhdr) > 0.5

    return outdata,outmask
