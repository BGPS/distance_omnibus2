from astroquery import ukidss
from astroquery.ukidss import Ukidss
from astropy import coordinates as coords
import pylab as pl
import sys
import os
sys.path.append(os.path.split(__file__)[0]+"/../nirex/")
from UKIDSS_extmap import show_contours_on_extinction,make_densitymap,get_contours
from astropy.io import fits

example_sources = ['G28.865+0.203', 'G29.841-0.476', 'G29.859-0.490',
                   'G30.397-0.466', 'G28.767-0.186', 'G29.277-0.131',
                   'G30.330+0.116']

coords = {x:coords.GalacticCoordinates(x[1:7],x[7:],unit=('deg','deg')) for x in example_sources}

C = get_contours(fits.open('v2.0_ds2_l029_13pca_map20_crop.fits'))

for c in coords:

    db = Ukidss.query_region(coords[c],
                             radius = '5 arcmin', programme_id='GPS',system='Galactic')

    db.write(c+'_jhk.fits')

    exm = make_densitymap(db,overwrite=True)

    pl.figure()

    show_contours_on_extinction(C,exm)
