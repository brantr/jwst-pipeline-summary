from astropy.io import fits
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

#read in the dia file
dia_fname = dia_dir + dia_file
hdu_dia = fits.open(dia_fname)

header_dia = hdu_dia[0].header
data_dia = hdu_dia[0].data

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