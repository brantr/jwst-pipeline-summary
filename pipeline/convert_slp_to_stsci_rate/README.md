This directory "v0" contains a first attempt at using the JWST STScI pipeline to reduce a single exposure.

I will attempt to reduce:

lux:/data/groups/comp-astro/brant/jades/Rev_run/F200W/slp/udf_cube_rev_F200W_481_289.slp.fits

I've placed this file in slp/.

First, let's look at the slp data:

python3:

```python
from astropy.io import fits
name = "slp/udf_cube_rev_F200W_481_289.slp.fits"
hdu = fits.open(fname)
hdu.info()
Filename: slp/udf_cube_rev_F200W_481_289.slp.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU     302   (2048, 2048, 2)   float32   

>>> header = hdu[0].header
>>> header
SIMPLE  =                    T / file does conform to FITS standard             
BITPIX  =                  -32 / number of bits per data pixel                  
NAXIS   =                    3 / number of data axes                            
NAXIS1  =                 2048 / length of data axis 1                          
NAXIS2  =                 2048 / length of data axis 2                          
NAXIS3  =                    2 / length of data axis 3                          
EXTEND  =                    T / FITS dataset may contain extensions            
COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H 
COMMENT Signal Ouput image.  Image consists of two planes.                      
COMMENT         Plane 1: Signal in ADU/sec                                      
COMMENT         Plane 2: Signal uncertainty in ADU/sec (noise model)            
COMMENT <----------------------------------------------------------->           
HISTORY Processed by NIRCam analysis software developed at UofA                 
NCDVERS = 'v2.0 [rc111] (16 Dec 2019)' / Version of NCDHAS software used        
NCDHPROC=                    T / Processed with NCDHAS                          
COMMENT Contact for error reporting: /dev/null                                  
COMMENT Contact for fearture requests: /dev/null                                
COMMENT Or misselt@as.arizona.edu might get some response.                      
COMMENT <----------------------------------------------------------->           
COMMENT <----------------------------------------------------------->           
COMMENT                    CALIBRATIONS APPLIED                                 
COMMENT <----------------------------------------------------------->           
DOREFCOR=                    T / Reference correction done                      
DOIPC   =                    T / IPC deconvolution performed                    
DOBPM   =                    T / Bad Pixel Mask applied                         
DOSAT   =                    T / Well depth detection done                      
DODARK  =                    F / Dark correction applied                        
DARKFBYF=                    F / Frame by frame dark correction                 
DOBIAS  =                    T / Bias correction applied                        
DOLIN   =                    T / Linearity correction applied                   
DOFLAT  =                    F / Flat field applied                             
CRDETECT=                    T / Cosmic ray detection done                      
PPIXGAIN=                    T / Gain per pixel used                            
COMMENT <----------------------------------------------------------->           
COMMENT                   CALIBRATION FILES USED                                
COMMENT <----------------------------------------------------------->           
DEF_CONF= 'NIRCam.cfg'         /  Configuration file used for defaults          
SCA_CONF= 'NRCA1_17004_SW_ISIMCV3.cfg' /  SCA Specific configuration file used  
BPMFILE = 'NRCA1_17004_BPM_ISIMCV3_2016-01-21.fits' /  Bad Pixel Mask file used 
IPCKERNL= 'NRCA1_17004_IPCDeconvolutionKernel_2016-03-18.fits' /  IPC kernel use
BIASFILE= 'NRCA1_17004_Bias_ISIMCV3_2016-02-09.fits' /  Bias frame used         
SATFILE = 'NRCA1_17004_WellDepthADU_2016-03-10.fits' /  Well depth file used    
LINFILE = 'NRCA1_17004_LinearityCoeff_2016-03-02.fits' /  Linearity file used   
GAINMAP = 'NRCA1_17004_Gain_ISIMCV3_2016-01-23.fits' /  Gain map used           
COMMENT <----------------------------------------------------------->           
COMMENT               REFERENCE PIXEL PROCESSING KEYS                           
COMMENT <----------------------------------------------------------->           
REFCOR  =                    T / Reference correction applied.                  
SREFCOR =                    F / Side reference pixels correction done.         
HFREFCOR=                    F / High Freq side reference done.                 
RCOLCOR =                    F / Ref pix correction on column basis             
EXCLBOTM=                    F / Exclude bottom reference pixels, 1st frame.    
RSIGREJ =                    T / Sigma rejection is determining ref pix         
RNITER  =                    3 / Number of sigma rejection iterations           
RSIGHI  =                   3. / Upper sigma clip value                         
RSIGLO  =                   3. / Lower sigma clip value                         
COMMENT <----------------------------------------------------------->           
COMMENT                     SIGNAL FITTING KEYS                                 
COMMENT <----------------------------------------------------------->           
ISCDS   =                    F / CDS analysis applied rather than linear fit    
IDRPFRM =                    0 / Number of frames dropped at begining of ramp.  
NFRMFIT =                  -99 / Number of frames to include in fit - override s
FITORDER=                    1 / Order of ramp fit.                             
FITMINIM=                    3 / Minimum number of frames for valid fit.        
COMMENT <----------------------------------------------------------->           
COMMENT              Header copied from raw input ramp                          
COMMENT <----------------------------------------------------------->           
COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H 
DATE    = '2019-12-11T12:35:21.430' / Date this file was created (UTC)          
ORIGIN  = 'AZ_LAB  '           / The organization that created the file         
TIMESYS = 'UTC     '           / Reference time                                 
FILENAME= '/home/marcia/image_create/udf_cube_rev_F' / file name                
FILETYPE= 'UNCALIBRATED'       / type of data in the file                       
TELESCOP= 'Guitarra'           / The overhead-free telescope                    
                                                                                
         Infomation about the coordinates in the file                           
                                                                                
RADESYS = 'FK5     '           / Name of coordinate reference frame             
TITLE   = 'NIRCam mocks'       / Proposal title                                 
PI_NAME = 'Zebigbos'           / Principal investigator name                    
CATEGORY= 'GTO     '           / Program category                               
SUBCAT  = 'NIRCAM  '           / Program sub-category                           
SCICAT  = 'Extragalac'         / Science category assigned by TAC               
CONT_ID =                    0 / continuation of previous program               
                                                                                
         Observation identifiers                                                
                                                                                
DATE-OBS= '2019-12-11'         / UTC date of start of exposure                  
TIME-OBS= '12:35:21.430'       / UTC time of start of exposure                  
DATE-END= '2019-12-11'         / UTC date of end of exposure                    
TIME-END= '12:58:15.735 '      / UTC time of end  of exposure                   
OBS_ID  = ''                   / Programmatic observation identifier            
VISIT_ID= '        '           / Visit identifier                               
PROGRAM = '1180    '           / Program number                                 
OBSERVTN= '7       '           / Observation number                             
VISIT   = ''                   / Visit number                                   
OBSLABEL= '        '           / Proposer label for the observation             
VISITGRP= '        '           / Visit group identifier                         
SEQ_ID  = '        '           / Parallel sequence identifier                   
ACT_ID  = '        '           / Activity identifier                            
EXPOSURE= '1234567 '           / Exposure request number                        
BKGDTARG=                    F / Background Target                              
TEMPLATE= 'NIRCam Imaging'     / Observation template used                      
ENG_QUAL= 'SUSPECT '           / engineering data quality indicator             
                                                                                
         Visit Information                                                      
                                                                                
VISITYPE= 'PRIME_TARGETED_FIXED' / Visit type                                   
VSTSTART= ''                   / UTC visit start                                
VISITSTA= ''                   / Visit status                                   
NEXPOSUR=                    0 / Total number of planned exposures              
INTARGET=                    F / At least one exposure is internal              
TARGOOPP=                    F / Visit scheduled as TOO                         
                                                                                
         Target information                                                     
                                                                                
TARGPROP= 'POINTINGONE-B'      / proposer's name for the target                 
TARGNAME= 'A Really Cool Field' / Standard astronomical catalog name            
TARGTYPE= 'FIXED   '           / Target type (fixed, moving, generic)           
TARG_RA =            53.157833 / Target RA at mid time of exposure (deg)        
TARG_DEC=           -27.807544 / Target Dec at mid time of exposure (deg)       
TARGURA =   1.000000000000E-08 / target RA uncertainty at mid-exposure          
TARGUDEC=   1.000000000000E-08 / target DEC uncertainty at mid-exposure         
MU_RA   =   0.000000000000E+00 / proper motion of the target in RA from A       
MU_DEC  =   0.000000000000E+00 / proper motion of the target in DEC from        
MU_EPOCH=   2.000000000000E+03 / EPOCH of proper motion values                  
PROP_RA =   0.000000000000E+00 / Target proper motion in RA                     
PROP_DEC=   0.000000000000E+00 / Target proper motion in DEC                    
                                                                                
         Instrument configuration                                               
                                                                                
INSTRUME= 'NIRCAM  '           / Instrument used to acquire data                
DETECTOR= 'NRCA1   '           / SCA name                                       
MODULE  = 'A       '           / NIRCam module:  A  or B                        
CHANNEL = 'SHORT   '           / NIRCam channel: short or long                  
FILTER  = 'F200W   '           / Name of filter element used                    
PUPIL   = 'CLEAR   '           / Name of pupil element used                     
PILIN   =                    F / Pupil imaging lens in the optical path         
CORONMSK= 'NONE    '           / coronagraph mask used                          
                                                                                
         Exposure parameters                                                    
                                                                                
EFFEXPTM=            1374.3053 / Effective Exposure time (sec)                  
DURATION=            1374.3053 / Total Duration of Exposure time (sec)          
EXPCOUNT=                    1 / Running count of exposures in visit            
EXPRIPAR= 'PRIME   '           / prime or parallel exposure                     
EXPSTART=    58828.52455358797 / UTC exposure start time (MJD)                  
EXPMID  =    58828.53250674352 / UTC exposure mid   time (MJD)                  
EXPEND  =    58828.54045989908 / UTC exposure end   time (MJD)                  
EXP_TYPE= 'NRC_IMAGE'          /  exposure type                                 
READPATT= 'DEEP8   '           /  detector read-out pattern                     
NINTS   =                    1 / Number of integrations in exposure             
NGROUPS =                    7 / Number groups in an integration                
NFRAMES =                    8 / Number of frames in group                      
GROUPGAP=                   12 / Number of frames skipped                       
NSAMPLES=                    1 / Number of A/D samples per pixel                
TSAMPLE =                  10. / Time  between samples in microsec              
TFRAME  =             10.73676 / Time in seconds between frames                 
TGROUP  =             214.7352 / Delta time between groups                      
EFFINTTM=            1288.4112 / Effective integration time (sec)               
NRSTSTRT=                    1 / Number of resets at start of exposure          
NRESETS =                    0 / Number of resets between integrations          
ZEROFRAM=                    F / Zero frame was downlinked separately           
DATAPROB=                    T / Science telemetry indicated a problem          
SCA_NUM =                  481                                                  
DRPFRMS1=                    0 / Number of frames dropped prior to 1st in       
DRPFRMS3=                    0 / Number of frames dropped between integra       
PNTG_SEQ=                    1 / Pointing sequence number                       
TSOVISIT=                    F / Time  Series Observation visit indicator       
                                                                                
         Subarray parameters                                                    
                                                                                
SUBARRAY= 'FULL    '           / Subarray used                                  
SUBSTRT1=                    1 / Starting pixel in axis 1 direction             
SUBSTRT2=                    1 / Starting pixel in axis 2 direction             
SUBSIZE1=                 2048 / Number of pixels in axis 1 direction           
SUBSIZE2=                 2048 / Number of pixels in axis 2 direction           
FASTAXIS=                    1 / Fast readout axis direction                    
SLOWAXIS=                    1 / slow readout axis direction                    
                                                                                
         NIRCam dither information                                              
                                                                                
PATTTYPE= 'MOSAIC  '           / Primary dither pattern type                    
PATT_NUM=                    1 / Position number in primary dither              
NUMDTHPT=                    1 / Total number of positions in pattern           
SUBPXNUM=                    1 / Subpixel pattern number                        
SUBPXPNS=                    9 / Total number of points in subpixel patte       
XOFFSET =                   0. / x offset from pattern starting position        
YOFFSET =                   0. / y offset from pattern starting position        
                                                                                
         JWST Ephemeris                                                         
                                                                                
REFFRAME= 'EME2000 '           / Ephemeris reference frame (EME2000)            
EPH_TYPE= 'Definitive'         / Definitive or Predicted                        
JWST_X  =        0.0000000E+00 / X spatial coordinate of JWST                   
JWST_Y  =        0.0000000E+00 / Y spatial coordinate of JWST                   
JWST_Z  =        0.0000000E+00 / Z spatial coordinate of JWST                   
JWST_DX =        0.0000000E+00 / X component of JWST velocity                   
JWST_DY =        0.0000000E+00 / Y component of JWST velocity                   
JWST_DZ =        0.0000000E+00 / Z component of JWST velocity                   
                                                                                
         Aperture information                                                   
                                                                                
APERNAME= 'NRCALL_FULL'        / Science Aperture name                          
PA_APER =        2.6000000E+01 / Position angle of aperture used                
PPS_APER= ''                   / original AperName supplied by PPS              
                                                                                
         Dither information                                                     
                                                                                
DVA_RA  =     0.0000000000E+00 / Velocity aberration correction RA offset       
DVA_DEC =     0.0000000000E+00 / Velocity aberration correction Dec offse       
VA_SCALE=     0.0000000000E+00 / Velocity aberration scale factor               
                                                                                
         TIME   information                                                     
                                                                                
BARTDELT=                   0. / Barycentric time correction (sec)              
BSTRTIME=        58828.5245536 / Barycentric exposure start time                
BENDTIME=        58828.5404599 / Barycentric exposure end time                  
BMIDTIME=        58828.5325067 / Barycentric exposure mid time                  
HELIDELT=                   0. / Heliocentric time correction (sec)             
HSTRTIME=        58828.5245536 / Heliocentric exposure start time               
HENDTIME=        58828.5404599 / Heliocentric exposure end time                 
HMIDTIME=        58828.5325067 / Heliocentric exposure mid time                 
                                                                                
         Photometry parameters                                                  
                                                                                
PHOTMJSR=     2.4696923105E+05 / Flux density (MJy/steradian) producing 1       
PHOTUJA2=     5.8048711366E+06 / Flux density (uJy/arcsec2) producing 1 c       
PIXAR_SR=     2.3619367204E-14 / Nominal pixel area in steradians               
PIXAR_A2=     1.0048900000E-03 / Nominal pixel area in steradians               
PHOTPLAM=     1.9886482775E+04 / Pivot wavelength Angstroms                     
PHOTFLAM=              26.3451 / Flux for 1e s-1 in erg cm**-2 s-1 A-1          
STMAG   =               30.798 / STMAG zeropoint                                
ABMAG   =              27.9973 / ABMAG zeropoint                                
                                                                                
                                                                                
PA_V3   =         6.95258E-310 / Position angle of V3-axis of JWST              
V2_REF  =              -0.3174 / V2 coordinate of ref point (arcsec)            
V3_REF  =            -492.5913 / V3 coordinate of ref point (arcsec)            
                                                                                
         WCS information                                                        
                                                                                
WCSAXES =                    2                                                  
CRPIX1  =               1024.5 / Axis 1 coordinate of reference pixel           
CRPIX2  =               1024.5 / Axis 2 coordinate of reference pixel           
CRPIX3  =                   0. / Axis 3 coordinate of reference pixel           
CRVAL1  =        53.1871366141 / RA at reference pixel (degrees)                
CRVAL2  =        -27.830983877 / DEC at reference pixel (degrees)               
CRVAL3  =           -171.78816 / T at reference pixel (seconds)                 
CTYPE1  = 'RA---TAN'           / Projection type                                
CTYPE2  = 'DEC--TAN'           / Projection type                                
CTYPE3  = '        '           / Projection type                                
CDELT1  =   8.805555555556E-06 / First axis increment per pixel                 
CDELT2  =   8.805555555556E-06 / Second axis increment per pixel                
CDELT3  =   2.147352000000E+02 / Third axis increment per pixel                 
CUNIT1  = 'deg     '           / First axis units                               
CUNIT2  = 'deg     '           / Second axis units                              
CUNIT3  = 'sec     '           / Third axis units                               
PC1_1   =     8.9303010626E-01                                                  
PC1_2   =     4.2403192369E-01                                                  
PC2_1   =    -4.2364855349E-01                                                  
PC2_2   =     8.9383823160E-01                                                  
PC3_1   =     0.0000000000E+00                                                  
PC3_2   =     0.0000000000E+00                                                  
RA_REF  =            53.157833 / RA  of the reference point (deg)               
DEC_REF =           -27.807544 / Dec of the reference point (deg)               
ROLL_REF=                  26. / Telescope roll angle of V3 at ref point        
VPARITY =                   -1 / Relative sense of rotation between Ideal       
V3I_YANG=           0.0000E+00 / Angle from V3 axis to Ideal y axis (deg)       
EQUINOX =                2000.                                                  
                                                                                
         Simulator parameters                                                   
                                                                                
INC_KTC =                    1 / include KTC F(0) T(1)                          
INC_RON =                    1 / include readnoise F(0) T(1)                    
INC_BKG =                    1 / include background F(0) T(1)                   
INC_CR  =                    1 / include Cosmic Rays F(0) T(1)                  
INC_DARK=                    1 / include darks F(0) T(1)                        
INC_LAT =                    0 / include latents F(0) T(1)                      
INC_NLIN=                    1 / include non-linearity F(0) T(1)                
NOISELES=                    F / NOISELES (T or F)                              
KTC     =               38.295 / KTC value (e-)                                 
BIAS    =               12213. / BIAS value (e-)                                
RDNOISE =                 11.3 / Readout noise (e-)                             
BKG     =            0.1539636 / background (e-/sec/pixel)                      
NSAMPLE =                    1 / background (e-/sec/pixel)                      
NINT    =                    1 / background (e-/sec/pixel)                      
HISTORY Science data not written by FITSWriter                                  
                                                                                
         Keywords required for DHAS                                             
                                                                                
SUBARRAY=                    F                                                  
COLCORNR=                    1                                                  
ROWCORNR=                    1                                                  
BREFROW =                    4 / bottom reference pixel rows                    
TREFROW =                    4 / top reference pixel rows                       
DRPFRMS1=                    0 / Number of frames skipped prior to first        
NGROUP  =                    7 / Number groups for ncdhas                       
NFRAME  =                    8 / Number of frames for ncdhas                    
SCA_ID  =                  481   

```

We need to convert this to the JWST image data format.  Let's write a python script to do that.

