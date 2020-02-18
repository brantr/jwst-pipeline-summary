from astropy.io import fits
from jwst.datamodels import ImageModel
from jwst.assign_wcs.assign_wcs import load_wcs
from astropy import units as u
from astropy.coordinates import SkyCoord

import numpy as np
import sys

#slp dir and filename
slp_dir = "slp/"
slp_file = "udf_cube_rev_F200W_481_289.slp.fits"

#dia dir and filename
dia_dir = "dia/"
dia_file = "udf_cube_rev_F200W_481_289.dia.fits"

#output dir and filename for rate file
output_dir = "rate/"
output_file = "udf_cube_rev_F200W_481_289_rate.fits"

#read in the slp file
slp_fname = slp_dir + slp_file
hdu_slp = fits.open(slp_fname)

header_slp = hdu_slp[0].header
data_slp   = hdu_slp[0].data

#convert from FK5 to ICRS
ra_ref     = header_slp['RA_REF']
dec_ref    = header_slp['DEC_REF']
crval1_ra  = header_slp['CRVAL1']
crval2_dec = header_slp['CRVAL2']
targ_ra    = header_slp['TARG_RA']
targ_dec   = header_slp['TARG_DEC']

c_ref   = SkyCoord(ra=ra_ref*u.degree, dec=dec_ref*u.degree, frame='fk5')
c_crval = SkyCoord(ra=crval1_ra*u.degree, dec=crval2_dec*u.degree, frame='fk5')
c_targ  = SkyCoord(ra=targ_ra*u.degree, dec=targ_dec*u.degree, frame='fk5')
header_slp['RA_REF']  = c_ref.icrs.ra.degree
header_slp['DEC_REF'] = c_ref.icrs.dec.degree
header_slp['CRVAL1']  = c_crval.icrs.ra.degree
header_slp['CRVAL2']  = c_crval.icrs.dec.degree
header_slp['TARG_RA']   = c_targ.icrs.ra.degree
header_slp['TARG_DEC']  = c_targ.icrs.dec.degree
header_slp['RADESYS']  = 'ICRS'

#fix pettern type
header_slp['PATTTYPE']  = 'IMAGING'
del header_slp['MU_RA']
del header_slp['MU_DEC']
header_slp['MU_EPOCH'] = '2000'


#read in the dia file
dia_fname = dia_dir + dia_file
hdu_dia = fits.open(dia_fname)

header_dia = hdu_dia[0].header
data_dia = hdu_dia[0].data


#convert from FK5 to ICRS
ra_ref     = header_dia['RA_REF']
dec_ref    = header_dia['DEC_REF']
crval1_ra  = header_dia['CRVAL1']
crval2_dec = header_dia['CRVAL2']
targ_ra    = header_dia['TARG_RA']
targ_dec   = header_dia['TARG_DEC']

c_ref   = SkyCoord(ra=ra_ref*u.degree, dec=dec_ref*u.degree, frame='fk5')
c_crval = SkyCoord(ra=crval1_ra*u.degree, dec=crval2_dec*u.degree, frame='fk5')
c_targ  = SkyCoord(ra=targ_ra*u.degree, dec=targ_dec*u.degree, frame='fk5')
header_dia['RA_REF']  = c_ref.icrs.ra.degree
header_dia['DEC_REF'] = c_ref.icrs.dec.degree
header_dia['CRVAL1']  = c_crval.icrs.ra.degree
header_dia['CRVAL2']  = c_crval.icrs.dec.degree
header_dia['TARG_RA']   = c_targ.icrs.ra.degree
header_dia['TARG_DEC']  = c_targ.icrs.dec.degree
header_dia['RADESYS']  = 'ICRS'

header_dia['PATTTYPE']  = 'IMAGING'
#del header_dia['MU_RA']
#del header_dia['MU_DEC']
header_dia['MU_EPOCH'] = '2000'
#create the data files for SCI and ERR

data_sci = data_slp[0,:,:]
data_err = data_slp[1,:,:]
data_dq  = data_dia[0,:,:].astype(np.uint32)

#create the primary header
header_rate_primary = header_slp.copy()

#remove NAXIS1, NAXIS2, NAXIS3 key words
del header_rate_primary['NAXIS1']
del header_rate_primary['NAXIS2']
del header_rate_primary['NAXIS3']
header_rate_primary['NAXIS'] = 0

#create the sci header
header_sci = header_slp.copy()

#adjust the axis key words
del header_sci['NAXIS3']
header_sci['NAXIS'] = 2
header_sci['EXTNAME'] = 'SCI'

#create the err header
header_err = header_slp.copy()

#adjust the axis key words
del header_err['NAXIS3']
header_err['NAXIS'] = 2
header_err['EXTNAME'] = 'ERR'


#adjust the DQ header
header_dq = header_dia.copy()
del header_dq['NAXIS3']
header_dq['NAXIS'] = 2
header_dq['EXTNAME'] = 'DQ'

header_vp = header_err.copy()
header_vp['EXTNAME'] = 'VAR_POISSON'

header_vr = header_err.copy()
header_vr['EXTNAME'] = 'VAR_RNOISE'

# create the rate file


# HERE are the components of the rate file
#HDU 	EXTNAME 		HDU Type 	Data Type 	Dimensions
#0		N/A 				primary 	N/A 				N/A
#1		SCI 				IMAGE 		float32 		ncols x nrows
#2 		ERR 				IMAGE 		float32 		ncols x nrows
#3 		DQ 					IMAGE 		uint32 			ncols x nrows
#4 		VAR_POISSON	IMAGE			float32			ncols x nrows
#5 		VAR_RNOISE 	IMAGE 		float32 		ncols x nrows
#6 		ASDF 				BINTABLE 	N/A 				variable

hdul_rate = fits.HDUList()

#create the primary header
hdul_rate.append(fits.PrimaryHDU(header=header_rate_primary))

#create the sci data and header
hdul_rate.append(fits.ImageHDU(header=header_sci, data=data_sci))

#create the err data and header
hdul_rate.append(fits.ImageHDU(header=header_err, data=data_err))

#create the dq data and header
hdul_rate.append(fits.ImageHDU(header=header_dq, data=data_dq))

#create the var_poisson extension
data_poisson = data_err**2

#add the var_poisson data and header
hdul_rate.append(fits.ImageHDU(header=header_vp, data=data_poisson))

#create the readnoise extension
data_rnoise = np.zeros_like(data_err)

#add the var_poisson data and header
hdul_rate.append(fits.ImageHDU(header=header_vr, data=data_rnoise))

hdul_rate.info()

#save the file
hdul_rate.writeto('test_rate.fits',overwrite=True)

from astropy.wcs import WCS
w = WCS(slp_fname)
im = ImageModel(slp_fname)
#ra, dec, lam = im.meta.wcs(x, y)
#print(ra,dec,lam)

im.set_fits_wcs(w)
test = im.find_fits_keyword('V2_REF')
print(test)
test = im.find_fits_keyword('TELESCOP')
print(test)
test = im.find_fits_keyword('NAXIS')
print(test)
im.meta.wcsinfo.v2_ref = header_slp['V2_REF']
print(im.meta.wcsinfo.v2_ref)
im.meta.wcsinfo.v3_ref = header_slp['V3_REF']
print(im.meta.wcsinfo.v3_ref)
im.meta.wcsinfo.pa_v3 = header_slp['PA_V3']
print(im.meta.wcsinfo.pa_v3)
im.meta.wcsinfo.ra_ref = header_slp['RA_REF']
im.meta.wcsinfo.dec_ref = header_slp['DEC_REF']
im.meta.wcsinfo.roll_ref = header_slp['ROLL_REF']



references={"area":"crds/jwst_nircam_area_0017.fits",  "distortion":"crds/jwst_nircam_distortion_0093.asdf",  "drizpars":"crds/jwst_nircam_drizpars_0001.fits",  "flat":"crds/jwst_nircam_flat_0337.fits",  "photom":"crds/jwst_nircam_photom_0074.fits"}
load_wcs(im,references)

im.data = hdul_rate['SCI'].data
im.err  = hdul_rate['ERR'].data
im.dq   = hdul_rate['DQ'].data
im.var_poisson = hdul_rate['VAR_POISSON'].data
im.var_rnoise  = hdul_rate['VAR_RNOISE'].data

im.meta.target.proper_motion_epoch="2000"
im.meta.dither.primary_type = "IMAGING"

ra_ref     = header_slp['RA_REF']
dec_ref    = header_slp['DEC_REF']
crval1_ra  = header_slp['CRVAL1']
crval2_dec = header_slp['CRVAL2']
targ_ra    = header_slp['TARG_RA']
targ_dec   = header_slp['TARG_DEC']

c_ref   = SkyCoord(ra=ra_ref*u.degree, dec=dec_ref*u.degree, frame='fk5')
c_crval = SkyCoord(ra=crval1_ra*u.degree, dec=crval2_dec*u.degree, frame='fk5')
c_targ  = SkyCoord(ra=targ_ra*u.degree, dec=targ_dec*u.degree, frame='fk5')
header_slp['RA_REF']  = c_ref.icrs.ra.degree
header_slp['DEC_REF'] = c_ref.icrs.dec.degree
header_slp['CRVAL1']  = c_crval.icrs.ra.degree
header_slp['CRVAL2']  = c_crval.icrs.dec.degree
header_slp['TARG_RA']   = c_targ.icrs.ra.degree
header_slp['TARG_DEC']  = c_targ.icrs.dec.degree
header_slp['RADESYS']  = 'ICRS'
im.meta.coordinates.reference_frame = 'ICRS'
im.meta.telescope = 'JWST'
im.meta.naxis= 2
im.save('test_im_rate.fits')
