from spectrum_extraction import spec_extract, bgps_to_mask, get_velocity_array
from astropy.io import fits
import pylab as pl
import numpy as np

dirname = 'edmonton2013/'

grsfile = dirname + 'grs-30-cube.fits'
vgpsfile = dirname + 'MOS_029.Tb.fit'
bgpsprefix = dirname + 'v2.0_ds2_l030_13pca'
bgpsfile = bgpsprefix + "_map20.fits"
bgpsmaskfile = bgpsprefix + "_labelmask.fits"

# define velocity arrays
vgpsh = fits.getheader(vgpsfile)
vgpsv = get_velocity_array(vgpsh)
grsh = fits.getheader(grsfile)
grsv = get_velocity_array(grsh)

# cache the data...
if not 'vgps' in locals():
    vgps = fits.getdata(vgpsfile).squeeze()
    grs = fits.getdata(grsfile)

def make_spectra(mask_number, npix=(1,2,3), vrange=[0,100], fignum=1, grsscale=10):

    datav,maskv = bgps_to_mask(bgpsfile,bgpsmaskfile,vgpsfile,mask_number)
    datag,maskg = bgps_to_mask(bgpsfile,bgpsmaskfile,grsfile,mask_number)

    nplots = len(npix)

    pl.figure(fignum)
    pl.clf()

    for ii,n in enumerate(npix):
        ax = pl.subplot(nplots,1,ii+1)
        sv,bv = spec_extract(vgps,maskv*datav,npix=n)
        sg,bg = spec_extract(grs,maskg*datag,npix=n)

        OKv = (vgpsv>vrange[0]*1000)*(vgpsv<vrange[1]*1000)
        OKg = (grsv>vrange[0]*1000)*(grsv<vrange[1]*1000)

        ax.plot(vgpsv[OKv],(sv-bv)[OKv])
        ax.plot(grsv[OKg],(sg-bg)[OKg]*grsscale)
        ax.set_title('Dilation = %i' % n)
        ax.set_ylabel('$T_{MB}$ (K)')
        if ii+1 < nplots:
            ax.set_xticks([])

    ax.set_xlabel("$V_{LSR}$ (m/s)")
    pl.subplots_adjust(hspace=0)

def do_all():
    for mask_number in np.unique(fits.getdata(bgpsmaskfile)):
        if mask_number > 0:
            make_spectra(mask_number)
            pl.savefig("GRS_HI_%i.pdf" % mask_number, bbox_inches='tight')
