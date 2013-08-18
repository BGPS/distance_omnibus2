from astroquery import ukidss
from astroquery.ukidss import Ukidss
from astropy import coordinates as coords
import pylab as pl
import sys
import os
sys.path.append(os.path.split(__file__)[0]+"/../nirex/")
from UKIDSS_extmap import show_contours_on_extinction,make_densitymap,get_contours,get_data
from astropy.io import fits
from astropy.table import Table

example_sources = ['G28.865+0.203', 'G29.841-0.476', 'G29.859-0.490',
                   'G30.397-0.466', 'G28.767-0.186', 'G29.277-0.131',
                   'G30.330+0.116']

coords = {x:coords.GalacticCoordinates(x[1:7],x[7:],unit=('deg','deg')) for x in example_sources}

C = get_contours(fits.open('v2.0_ds2_l029_13pca_map20_crop.fits'))
C30 = get_contours(fits.open('v2.0_ds2_l030_13pca_map20_crop.fits'))
label30 = fits.open('v2.0_ds2_l030_13pca_labelmask.fits')
label30[0].data = np.array((label30[0].data > 0),dtype='float')
C30l = get_contours(label30, contour_level=0.5)

r = 20 # arcmin
for c in coords:

    print "Loading coordinate ",c,'...',

    l,b = coords[c].l.degree,coords[c].b.degree

    fn = "%0.3f%+0.3f_r%0.1f_catalog.fits" % (l,b,r)
    if os.path.exists(fn):
        db = Table.read(fn)
    else:
        db = get_data(l,b, radius=r, get_images=False)

    #db = Ukidss.query_region(coords[c],
    #                         radius = '15 arcmin', programme_id='GPS',system='Galactic')
    #
    #db.write(c+'_jhk.fits')

    print "Making densitymap...",
    layers = []
    ranges = range(11,20)
    for u,l in zip(ranges[1:],ranges[:-1]):
        exm = make_densitymap(db,overwrite=True,kband_lower=l,kband_upper=u)
        layers.append(exm)

    pl.rc('font',size=24)
    pl.close(0)
    pl.figure(0,figsize=(8,12))
    for ii in xrange(8):
        pl.clf()
        pl.title("K=%i to %i" % (11+ii,11+ii+1))
        show_contours_on_extinction(C30l,layers[ii],color='c')
        pl.xlabel("Galactic Longitude")
        pl.ylabel("Galactic Latitude")
        pl.axis( (30.215755981611451, 30.483123039075572, -0.16376485951249242, 0.16176321369869018) )
        pl.savefig("Kband%ito%i_BGPSLabelcontours.png" % (11+ii,11+ii+1),bbox_inches='tight')

        pl.clf()
        pl.title("K=%i to %i" % (11+ii,11+ii+1))
        pl.xlabel("Galactic Longitude")
        pl.ylabel("Galactic Latitude")
        show_contours_on_extinction(C30,layers[ii],color='y')
        pl.axis( (30.215755981611451, 30.483123039075572, -0.16376485951249242, 0.16176321369869018) )
        pl.savefig("Kband%ito%i_BGPScontours.png" % (11+ii,11+ii+1),bbox_inches='tight')

        

    pl.figure()

    if l < 30:
        show_contours_on_extinction(C,exm)
    else:
        show_contours_on_extinction(C30,exm)

    pl.show()
