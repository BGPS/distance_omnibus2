from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astrodendro.analysis import PPVStatistic
from astrodendro import Dendrogram
import numpy as np
import pylab as pl
import matplotlib as mpl
import itertools
import progressbar
from agpy import fit_a_line

pb = progressbar.ProgressBar()

# 306 = 60 km/s
cropn = 50
grs = fits.getdata('grs-28-cube.fits')[:306,cropn:-cropn,cropn:-cropn]
grsh = fits.getheader('grs-28-cube.fits')

# In [55]: (1*u.K).to(u.Jy,jytok((46*u.arcsec)**2*pi,110*u.GHz))
# Out[55]: <Quantity 58.0861309512 Jy>
metadata = dict(data_unit=u.Jy/58.06, spatial_scale=22.5*u.arcsec, velocity_scale=grsh['CDELT3']*u.m/u.s, vaxis=0, wcs=wcs.WCS(grsh))

if 'grsD' not in locals():
    grsD = Dendrogram.compute(data=grs, min_value=0.25, min_npix=10, min_delta=0.5, verbose=True)

def get_velo_stats(leaf):
    stats = []
    L = leaf
    while hasattr(L,'parent'):
        ps = PPVStatistic(L, metadata=metadata)
        stats.append(ps)
        L = L.parent
    return stats

def get_trunk(x):
    if x.parent is None:
        return x
    return get_trunk(x.parent)

L = grsD.leaves[0]

while hasattr(L,'parent'):
    print L.idx,L.parent,L.get_npix()
    L = L.parent

pb.currval = 0
leaf_trunk = {L.idx:get_trunk(L).idx for L in pb(grsD.leaves)}
pb.currval = 0
trunk_leaves = {}
for k,v in pb(leaf_trunk.iteritems()):
    if v in trunk_leaves:
        trunk_leaves[v].append(k)
    else:
        trunk_leaves[v] = [k]

pb.currval = 0
leaf_stats_full = {L.idx: get_velo_stats(L) for L in pb(grsD.leaves)}
# prune
pb.currval = 0
leaf_stats = {k:v for k,v in leaf_stats_full.iteritems() if v[0].v_cen < 60e3}

pb.currval = 0
leaf_radii = {i:[x.radius for x in leaf_stats[i]] for i in pb(leaf_stats)}
pb.currval = 0
leaf_vrms = {i:[x.v_rms for x in leaf_stats[i]] for i in pb(leaf_stats)}

pb.currval = 0
radii = [[x.value for x in i] for i in pb(leaf_radii.values())]
pb.currval = 0
vrms = [np.array([x.value for x in i])[np.argsort(r)] for i,r in pb(zip(leaf_vrms.values(),radii))]

pl.figure(2)
pl.clf()
for r,v in zip(radii,vrms):
    r = np.array(r)
    pl.plot(r,v,linestyle='-',marker='.',alpha=0.5)
    #plot(r[r<600],v[r<600],alpha=0.5)
    #plot([x.value for x in leaf_radii[i]],
    #     [x.value for x in leaf_vrms[i]],
    #     '-')

markerd = mpl.markers.MarkerStyle.markers.copy()
markerd.pop(None)
markerd.pop('k')
markers = itertools.cycle(markerd)

trunkcolors = {}

pl.figure(3)
pl.clf()
ax = pl.gca()
for trunk in trunk_leaves:
    if trunk_leaves[trunk][0] in leaf_stats:
        color = ax._get_lines.color_cycle.next()
        trunkcolors[trunk] = color
        for leaf in trunk_leaves[trunk]:
            r = np.array([x.value for x in leaf_radii[leaf]])
            v = np.array([x.value for x in leaf_vrms[leaf]])
            pl.plot(r, v, linestyle='-', marker=markers.next(), alpha=0.5, color=color)

trunkslopes = {}
xvals = {}

for trunk in trunk_leaves:
    if trunk_leaves[trunk][0] in leaf_stats:
        x = [j.value for leaf in trunk_leaves[trunk] for j in leaf_radii[leaf] if j.value > 20]
        y = [j.value for leaf in trunk_leaves[trunk] for j,r in zip(leaf_vrms[leaf],leaf_radii[leaf]) if r.value > 20]
        if len(x) <= 3:
            continue
        fitpars = fit_a_line.total_least_squares(np.array(x),np.array(y))
        trunkslopes[trunk] = fitpars
        xvals[trunk] = x

pl.figure(4)
pl.clf()
pl.figure(5)
pl.clf()
ax = pl.gca()
for trunk,fitpars in trunkslopes.iteritems():
    if fitpars[0] > 100:
        pl.figure(4)
    elif fitpars[0] < 2:
        continue
    else:
        pl.figure(5)
    pl.plot(xvals[trunk],np.polyval(fitpars, xvals[trunk]), color=trunkcolors[trunk], linewidth=2, alpha=0.5)

pl.figure(6)
pl.clf()
for jj,trunk in enumerate(xvals):
    if jj > 35:
        continue
    ax = pl.subplot(6,6,1+jj)
    R = [50,50]
    for leaf in trunk_leaves[trunk]:
        r = np.array([x.value for x in leaf_radii[leaf]])
        v = np.array([x.value for x in leaf_vrms[leaf]])
        ax.plot(r, v, linestyle='-', marker=markers.next(), alpha=0.5, color=trunkcolors[trunk], linewidth=0.5)
        R[0] = min([R[0],min(r)])
        R[1] = max([R[1],max(r)])
    ax.plot(R,np.polyval(trunkslopes[trunk], R), color='k', linestyle='--', linewidth=2, alpha=0.5)
    ax.set_yticks([])
    ax.annotate(str(trunkslopes[trunk][0]), (0.7,0.2), xycoords='axes fraction')
    ax.annotate(str(trunk), (0.7,0.1), xycoords='axes fraction')
pl.subplots_adjust(hspace=0,wspace=0,right=1,left=0,top=1,bottom=0)

shallow_slopes = [t for t,(m,b) in trunkslopes.iteritems() if ((m>5) & (m < 25))]
pl.figure(7)
pl.clf()
for jj,trunk in enumerate(shallow_slopes):
    if jj > 35:
        continue
    ax = pl.subplot(6,6,1+jj)
    R = [50,50]
    for leaf in trunk_leaves[trunk]:
        r = np.array([x.value for x in leaf_radii[leaf]])
        v = np.array([x.value for x in leaf_vrms[leaf]])
        ax.plot(r, v, linestyle='-', marker=markers.next(), alpha=0.5, color=trunkcolors[trunk], linewidth=0.5)
        R[0] = min([R[0],min(r)])
        R[1] = max([R[1],max(r)])

    ax.plot(R,np.polyval(trunkslopes[trunk], R), color='k', linestyle='--', linewidth=2, alpha=0.5)
    ax.set_yticks([])
    ax.annotate(str(trunkslopes[trunk][0]), (0.7,0.2), xycoords='axes fraction')
    ax.annotate(str(trunk), (0.7,0.1), xycoords='axes fraction')
pl.subplots_adjust(hspace=0,wspace=0,right=1,left=0,top=1,bottom=0)

viewer = grsD.viewer(fignum=1)

