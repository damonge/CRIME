import healpy as hp
import numpy as np

fwhmax=10.0
dtor=3.1415926536/180

map_synchrotron=hp.read_map("lambda_haslam408_dsds.fits")
hp.write_map("lambda_haslam408_dsds_ring.fits",map_synchrotron,coord='G',fits_IDL=False)

map_faraday=hp.read_map("faraday.fits",hdu=3)
map_sigma=np.sqrt(hp.smoothing(map_faraday**2,fwhm=fwhmax*dtor))
hp.write_map("faraday_sigmas.fits",map_sigma,coord='G',fits_IDL=False)
