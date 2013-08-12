import numpy as np
try:
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
except ImportError:
    import pyfits
    import pywcs
import scipy.ndimage

def hcongrid(image, header1, header2, **kwargs):
    """
    Interpolate an image from one FITS header onto another

    kwargs will be passed to `scipy.ndimage.map_coordinates`

    Parameters
    ----------
    image : ndarray
        A two-dimensional image 
    header1 : `pyfits.Header` or `pywcs.WCS`
        The header or WCS corresponding to the image
    header2 : `pyfits.Header` or `pywcs.WCS`
        The header or WCS to interpolate onto

    Returns
    -------
    ndarray with shape defined by header2's naxis1/naxis2

    Raises
    ------
    TypeError if either is not a Header or WCS instance
    Exception if image1's shape doesn't match header1's naxis1/naxis2

    Examples
    --------
    (not written with >>> because test.fits/test2.fits do not exist)
    fits1 = pyfits.open('test.fits')
    target_header = pyfits.getheader('test2.fits')
    new_image = hcongrid(fits1[0].data, fits1[0].header, target_header)

    """

    if issubclass(pywcs.WCS, header1.__class__):
        wcs1 = header1
    else:
        try:
            wcs1 = pywcs.WCS(header1)
        except:
            raise TypeError("Header1 must either be a pyfits.Header or pywcs.WCS instance")

    if not (wcs1.naxis1 == image.shape[1] and wcs1.naxis2 == image.shape[0]):
        raise Exception("Image shape must match header shape.")

    if issubclass(pywcs.WCS, header2.__class__):
        wcs2 = header2
    else:
        try:
            wcs2 = pywcs.WCS(header2)
        except:
            raise TypeError("Header2 must either be a pyfits.Header or pywcs.WCS instance")

    if not all([w1==w2 for w1,w2 in zip(wcs1.wcs.ctype,wcs2.wcs.ctype)]):
        # do unit conversions
        raise NotImplementedError("Unit conversions have not yet been implemented.")

    # sigh... why does numpy use matrix convention?  Makes everything so much harder...
    outshape = [wcs2.naxis2,wcs2.naxis1]
    yy2,xx2 = np.indices(outshape)
    lon2,lat2 = wcs2.wcs_pix2sky(xx2, yy2, 0)
    xx1,yy1 = wcs1.wcs_sky2pix(lon2, lat2, 0)
    grid1 = np.array([yy1.reshape(outshape),xx1.reshape(outshape)])

    newimage = scipy.ndimage.map_coordinates(np.nan_to_num(image), grid1, **kwargs)
    
    return newimage

def flatten_header(header):
    """
    Attempt to turn an N-dimensional fits header into a 2-dimensional header
    Turns all CRPIX[>2] etc. into new keywords with suffix 'A'

    header must be a pyfits.Header instance
    """

    # astropy.io.fits != pyfits -> sadness
    #if not hasattr(header,'copy')
    #    raise Exception("flatten_header requires a pyfits.Header instance")

    newheader = header.copy()

    for key in newheader.keys():
        try:
            if int(key[-1]) >= 3 and key[:2] in ['CD','CR','CT','CU','NA']:
                newheader.rename_key(key,'A'+key,force=True)
        except ValueError:
            # if key[-1] is not an int
            pass
        except IndexError:
            # if len(key) < 2
            pass
    newheader.update('NAXIS',2)

    return newheader

def get_velocity_array(header):
    for k in header:
        if header[k] == 'VELOCITY' and k[:-1] == 'CTYPE':
            vaxis = int(k[-1])

    crpix = header['CRPIX%i' % vaxis]
    naxis = header['NAXIS%i' % vaxis]
    crval = header['CRVAL%i' % vaxis]
    cdelt = header['CDELT%i' % vaxis] if 'CDELT%i' % vaxis in header else 'CD%i_%i' % (vaxis,vaxis)

    arr = (np.arange(naxis)+1-crpix)*cdelt + crval

    return arr
