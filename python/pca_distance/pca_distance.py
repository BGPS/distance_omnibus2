import numpy as np
from astropy.io import fits
import astropy.nddata as nd
import matplotlib.pyplot as p
import scipy.ndimage as ndim
import astropy.wcs as wcs
import scipy.stats as stats
def pca_distance(data, ChannelWidth = 1.0, PixelScale = 1.0, plot = False):
    nv = data.shape[0]
    nspatial = data.shape[1]*data.shape[2]
    nchan = data.shape[0]
    array = np.reshape(data,(data.shape[0],data.shape[1]*data.shape[2]))
    array = array - array[array==array].mean()
    array[np.isnan(array)]=0.0
    covar = np.dot(array,np.transpose(array))
    evals,evec = np.linalg.eig(covar)
# Build reconstructed array
    array = np.dot(np.transpose(evec),array)
    array.shape=(nv,data.shape[1],data.shape[2])
    for index in np.arange(nv):
        WtImage = array[index,:,:]
        AcorImage = nd.convolve_fft(WtImage,WtImage)
        array[index,:,:]=AcorImage

# Build a noise amplitude image
#   Some code goes here       

    xmat,ymat = np.meshgrid(np.arange(data.shape[2]),np.arange(data.shape[1]))
    Dr = np.zeros((nv))+np.nan
    for index in np.arange(nv):
        AcorImage = array[index,:,:]
        MaxPosition = np.unravel_index(np.argmax(AcorImage),AcorImage.shape)
        Distance = ((xmat-MaxPosition[1])**2 + (ymat-MaxPosition[0])**2)**0.5
        Binned = np.floor(Distance)
        Bins = np.unique(Binned)
        BinnedMeans = [AcorImage[Binned == i].mean() \
                           for i in Bins]
# Normalize to Max Value
        if BinnedMeans[0] > 0:
            BinnedMeans = BinnedMeans/BinnedMeans[0]
            BelowIndex = np.min(np.where(BinnedMeans<np.exp(-1)))
            AboveIndex = BelowIndex-1
            Dr[index] = (np.exp(-1)-BinnedMeans[AboveIndex])/\
                (BinnedMeans[AboveIndex]-BinnedMeans[BelowIndex])*\
                (Bins[AboveIndex]-Bins[BelowIndex])+Bins[AboveIndex]

    Dv = np.zeros((nv))+np.nan
    for index in np.arange(nv):
        vec = evec[:,index]
        AcorVec = np.correlate(vec,vec,mode='full')
        pkpix = np.argmax(AcorVec)
        DeltaV = np.abs(np.arange((nv))-pkpix)
        VBinned = np.floor(DeltaV)
        VBins = np.unique(VBinned)
        BinnedVMeans = [AcorVec[VBinned == i].mean() for i in VBins]
        if BinnedVMeans[0] > 0:
            BinnedVMeans/BinnedVMeans[0]
            BelowIndex = np.min(np.where(BinnedVMeans<np.exp(-1)))
            AboveIndex = BelowIndex-1
            Dv[index] = (np.exp(-1)-BinnedVMeans[AboveIndex])/\
                (BinnedVMeans[AboveIndex]-BinnedVMeans[BelowIndex])*\
                (VBins[AboveIndex]-VBins[BelowIndex])+VBins[AboveIndex]
            
#        Dv[index] = VBins[np.min(np.where(BinnedVMeans<np.exp(-1)))]


    UseIndex = np.where((Dv > 2)*(Dr > 2))   
    # Plslope,logIntercept,rval,pval,stderr = \
    #     stats.linregress(np.log(Dr[UseIndex]),np.log(Dv[UseIndex]))
    # Intercept = np.exp(logIntercept)
    Distance = ((Dv*ChannelWidth)/(0.87))**(1/0.67)\
        *(Dr*PixelScale*np.pi/180)**(-1.0)
    Distance = np.exp(np.mean(np.log(Distance[UseIndex])))
    if plot:
        p.loglog((Dr[UseIndex]*PixelScale*np.pi/180*Distance),\
                     ChannelWidth*Dv[UseIndex],marker='.',linestyle='None')
        p.xlabel('$\delta r$ (pc)')
        p.ylabel('$\delta v$ (km/s)')
        p.show()
    return(Distance)

def subcube_extract(glon, glat, vlsr, \
                        border = [30,50,50], threshold = 0.4,\
                        GRSfile='/UBC-O/erosolo/astro/bgps/scutum/edmonton2013/grs-30-cube.fits'):
    GRSCube = fits.getdata(GRSfile)
    GRSheader = fits.getheader(GRSfile)
    w = wcs.WCS(GRSheader)
    pixels = w.wcs_world2pix([[glon,glat,vlsr]],1)
# Make consistent with pixel ordering
    pixels = (pixels[0,::-1])

    x1 = int( np.floor( max([pixels[2]-border[2],0]) ) )
    x2 = int( np.ceil( min([pixels[2]+border[2],GRSCube.shape[2]]) ) )
    y1 = int( np.floor( max([pixels[1]-border[1],0]) ) )
    y2 = int( np.ceil( min([pixels[1]+border[1],GRSCube.shape[1]]) ) )
    v1 = int( np.floor( max([pixels[0]-border[0],0]) ) )
    v2 = int( np.ceil( min([pixels[0]+border[0],GRSCube.shape[0]]) ) )
    subcube = GRSCube[v1:v2,y1:y2,x1:x2]
    label,nFeature = ndim.label(subcube>threshold)
    hit = np.asarray(pixels-[v1,y1,x1],dtype='int8')

    mask = (label==label[hit[0],hit[1],hit[2]])
    
    subcube[np.logical_not(mask)]=np.nan
    return(subcube)

def DPDF_pcadist(distances,glon,glat,vlsr, \
                     GRSfile='/UBC-O/erosolo/astro/bgps/scutum/edmonton2013/grs-30-cube.fits'):
    subcube = subcube_extract(glon,glat,vlsr,GRSfile=GRSfile)
    hd = fits.getheader(GRSfile)
    PixelScale = hd['CDELT2']
    ChannelWidth = hd['CDELT3']/1e3
    Distance = pca_distance(subcube,PixelScale = PixelScale, ChannelWidth = ChannelWidth)
# tAKE WIDTH AS 20% OF DISTANCE DETERMINATION.  WITH CAPS LOCK ON
    dpdf = np.exp(-((distances-Distance)/(0.2*Distance))**2/2)
    return(dpdf)
