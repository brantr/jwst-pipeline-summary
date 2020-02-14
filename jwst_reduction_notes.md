# JWST Pipeline Reduction Notes

### Required software

* astropy -- available via pip
* glue -- available via pip
* ginga -- available via pip
* photutils -- available via pip
* astroimtools -- available via pip
* imexam -- available via pip
* asdf -- available via pip
* gwcs -- available via pip
* synphot -- available via pip
* jwst -- see below

#### ASDF depends on
* numpy -- available via pip
* jsonschema -- available via pip
* pyyaml -- available via pip
* six -- available via pip

### JWST Info

* the jwst package on pip is a dummy package that indicates no public release is available
* The STScI github site has the package available as [jwst](https://github.com/spacetelescope/jwst).  You can 

``` bash
git clone https://github.com/spacetelescope/jwst.git
```
Then, to make the package available via python, make a link to the $jwst/jwst$ directory in $lib/python3.6/site-packages/$.

* jwst depends on crds -- which is available via pip
* There is a (jwst pipeline introduction)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/introduction.html] on the STScI website.
* The strun script is in jwst/scripts.  Note that "python" make need to be changed to "python3".

### GCWS Info

* [STScI notebook on GWCS](https://github.com/spacetelescope/JWSTUserTraining2016/blob/master/Workshop_Notebooks/GWCS/GWCS.ipynb)
* The "forward" direction is from detector to sky.
* The WCS is 0-based.

### ASDF Info

* (The ASDF Read the Docs)[https://asdf-standard.readthedocs.io/en/latest/asdf_in_fits.html]

### Pipeline Stages

* (STScI Pipeline Stages)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/main.html#pipelines]
* (Stage 2 Imaging Processing)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_image2.html#calwebb-image2]

Stage 2 involves the (background)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/background/index.html#background-step], (assign_wcs)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/assign_wcs/index.html#assign-wcs-step], (flat_field)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/flatfield/index.html#flatfield-step], (photom)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/photom/index.html#photom-step], and (resample)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/index.html#resample-step] steps.

The input is a countrate exposure, either "\_rate" or "\_rateints" files, or an ASN file with multiple inputs. The outputs are "\_bsub" (if save_bsub==True) the data output from the background step, "\_cal" a fully calibrated but unrectified expsoure, and "\_i2d" a resampled data that is not permanent and is not passed along to Stage 3.

* (Stage 3 Imaging Processing)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_image3.html#calwebb-image3]

Stage 3 involves the combination of calibrated dat from multiple exposures into a single, rectified (distortion corrected) image.  The input exposures are corrected for astrometric alignment, background matching, and outlier rejection.  The steps are (tweakreg)[]https://jwst-pipeline.readthedocs.io/en/latest/jwst/tweakreg/index.html#tweakreg-step, (skymatch)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/skymatch/index.html#skymatch-step], (outlier_detection)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/outlier_detection/index.html#outlier-detection-step], (resample)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/index.html#resample-step], and (source_catalog)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/source_catalog/index.html#source-catalog-step] steps.

The input are the calwebb_image2 image products ("\_cal") or an ASN file.  If a single "\_cal" file is used as input, only the resample and source_catalog steps are performed.

The outputs are cosmic ray-flagged exposures ("\_crf"), which has an updated DQ array.  A resampled and combined image "\_i2d" from the resample step (note the Stage 2 "\_i2d" are usually not used).  And a source catalog in an ASCII file "ecsv" format ("\_cat").


## Stage 2

* (background)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/background/description.html] Subtracts background signal. Input is a single target exposure and one or more background exposures.  For multiple background exposures, averaged and then iterative sigma clipped.  The steps are clip SCI arrays of all background exposures, compute mean of unclipped sci values, sum the ERR arrays in quadratures for the bkg exposures clipping same input from the SCI arrays, convert to uncertainty in the mean, then combine the DQ arrays using bitwise OR.  The average bkg is then subtracted from the targets.  The bkg SCI is subtracted from the target SCI.  The target ERR is unchanged (full error propagation not yet implemented). DQ arrays of bkg and target expourse are combined with bitwise OR. The "S_BKDSUB" keyword in the output will read as COMPLETE.

* (assign_wcs)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/assign_wcs/main.html] Associates a WCS object with each science exposure.  There may be intermediate coordinate frames depending on the instrument. the WCS is saved in the ASDF extension of the FITS file. "It can be accessed as an attribute of the meta object when the fits file is opened as a data model."  Forward transforms is detector to world, and input positions or 0-based.  Expects the basic WCS keywords in the SCI header.  Distortion and spectral models are stored in the ASDF format.  For each mode in "EXP_TYPE" keyword, assign_wcs gets the reference files from CRDS and applies transforms from detector to v2v3, then the basic WCS keywords transforms from v2v3 to world.  The header keywords from v2v3 to world are (RA_REF, DEC_REF, V2_REF, V3_REF, ROLL_REF, RADESYS).  If the FITS is opened as a DataModel, then WCS can be accessed via model.meta.wcs:

```python
>>> from jwst.datamodels import ImageModel
>>> exp = ImageModel('miri_fixedslit_assign_wcs.fits')
>>> ra, dec, lam = exp.meta.wcs(x, y)
>>> print(ra, dec, lam)
    (329.97260532549336, 372.0242999250267, 5.4176100046836675)
```

* (flat_field)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/flatfield/main.html] The SCI array from the flat-field reference is divided into the SCI and ERR arrays, and the flat-field DQ array is combined with SCI DQ with bitwise OR.  For Imaging datasets, find pixels with NaN or zero in FLAT reference SCI array, and set DQ to "NO_FLAT_FIELD".  Reset the values of pixels in the flat with DQ="No_FLAT_FIELD" to 1.0.  Apply the flat by dividing into SCI and ERR. Update SCI DQ with FLAT DQ.  The total ERR array is updated as the square root of the quadratice sum of "VAR_POISSON", "VAR_RNOISE", "VAR_FLAT". 

* (photom)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/photom/main.html] Photometric calibration step.  Converts from units of countrate to surface brightness.  The calibration info is from a photometric reference file.  The assign_wcs must be applied before photom.  The "S_PHOTOM" keyword will read as "COMPLETE".  The "BUNIT" will be updated to reflect the change in units.  The PHOTOM reference files include a scalar flux conversion constant depending on detector, filter, pupil, and grating. The scalar conversion constant is copied to the header PHOTMJSR, gives conversion from DN/s to megaJy/steradian.  Stores PHOTUJA2 as well. The photom step also retrieves a pixel area map and copies into "AREA" extension in the science product. Also populates PIXAR_SR and PIXAR_A2.  These are averages.

* (resample)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/main.html] NOTE -- these data products at Stage 2 are not used in Stage 3. Resamples images based on WCS and distortion info, and produces a single undistorted product.  Either a single image or an association table.  The parameters for drizzle get provided by DRIZPARS reference from CRDS. Output data uses WCS of all inputs, the output WCS defines a field-of-view encompasseds all input images with same orientation and plate scale as the FIRST LISTED INPUT.  Uses the "C-based" cdriz routine. A full description drizzle is at (DrizzlePac Handbook)[http://www.stsci.edu/scientific-community/software/drizzlepac.html].

## Stage 3

* (tweakreg)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/tweakreg/README.html] Creates an image catalog of point-like sources, uses their centroids to compute corrections to WCS such that sky catalogs from input images will align.  Stars are found using photutils "daofind".  Linear coordinate transformations that align the star catalogs are defined. Tangent-plane corrections applied to GWCS pipeline are determined.  The correction is implemented using "TPCorr", and inserted into GWCS pipeline of image's WCS.  More info on tweakreg (here)[https://tweakwcs.readthedocs.io/en/latest/].

* (skymatch)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/skymatch/README.html] Computes the sky values in a collection of images that have both sky and signal.  Can compute image sky values separately, or in a way to minimize sky differences between images. Compares total signal levels in the overlap regions and computes signal offsets that will minimize the residuals across the whole. Performed directly without resampling. Assumes assign_wcs has been performed.  There are three versions of background computations.  "localmin" computes sky statistics for each image, incorporating DQ flags. "globalamin" computes the minimum sky value across all input images.  A single sky value is the background in all images. "match" algorithm computes a constant (within image) value to be applied to each input such that the mismatch in the computed background between all pairs of images is least squares minimized.  The globalmin+match uses globalmin to find baseline sky  and the match to equalize sky values among images.  NOTE: I am guessing that match is the correct algorithm.  DEFAULT is all non-zero DQ pixels are bad.  Performed on flattend but not distortion corrected(!) images.  Mask point like sources!

* (outlier_detection)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/outlier_detection/main.html] IDs bad pixels or cosmic rays in input images. Builds a stack of inmput data with same WCS, create median image, create blotted data for each original input, pixel-by-pixel stat comparison for difference from mean by some number of sigma.  Flag bad pixels in the DQ array.

* (resample)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/main.html] Resamples images based on WCS and distortion info, and produces a single undistorted product.  Either a single image or an association table.  The parameters for drizzle get provided by DRIZPARS reference from CRDS. Output data uses WCS of all inputs, the output WCS defines a field-of-view encompasseds all input images with same orientation and plate scale as the FIRST LISTED INPUT.  Uses the "C-based" cdriz routine. A full description drizzle is at (DrizzlePac Handbook)[http://www.stsci.edu/scientific-community/software/drizzlepac.html].

* (source_catalog)[https://jwst-pipeline.readthedocs.io/en/latest/jwst/source_catalog/main.html] Produces source photometry and morphologies.  Image segmentation is from photutils source extraction, based on threshold method. Deblending uses photutils (deblender)[https://photutils.readthedocs.io/en/latest/segmentation.html#source-deblending]. Properties computed for each object are fluxes and errors, ab mags and errors, area, semi major and semiminor, orientation, and sky coordinates at bounding box corners.