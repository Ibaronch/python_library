# python_library

Useful python 2.7 functions and classes for astronomical uses

Ivano Baronchelli 2018


# List of classes:
-------------------------------------------------------------------
 Deg2sex    Given RA and DEC in input, they will be returnet in
            output converted from 
            decimal to sexagesimal or vice-versa (depending on
            intype and outype variables set).
-------------------------------------------------------------------
 Match      matches a series of (integer) input IDs (like IDL match)
-------------------------------------------------------------------
 Match_cat  Matches coordinates [deg], in two lists of coordinates,
            using an user specyfied searching radius dt[arcsec]
            For each source in the first catalog, the closest
            counterpart in the second catalog is associated.
-------------------------------------------------------------------
 Readcol  Read specified columns in an ascii file
          and returns them in the specified format
-------------------------------------------------------------------


#  List of functions:
-------------------------------------------------------------------
 gaussfit        Simple mono-dimensional gaussian fit
-------------------------------------------------------------------
 image_depth     measure depth (1sigma) of an image.fits
                 using the random apertures method. 
-------------------------------------------------------------------
 match           matches a series of (integer) IDs
-------------------------------------------------------------------
 mk_regionfile   creates simple ds9 region files. 
-------------------------------------------------------------------
 image_shift     shift the RA, DEC center of an image
-------------------------------------------------------------------
