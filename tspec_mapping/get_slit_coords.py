"""
TripleSpec's slit is centered at coordinates
250,416.4 
REVISED:
    250.11, 416.53
    angle 0.66 degrees
    length 89.3 pix, height 3.24 pix
in the slitviewer (proc coordinates, not raw)
Its endpoints are approximately:
    207,416  # REVISED 205.6 415.7
    293,417  # REVISED 294.9 416.7

These are manually determined but it would be nice to use the mask to fit it
directly...
"""

import numpy as np
try:
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
except ImportError:
    import pywcs
    import pyfits

def get_slit_center(finder):
    header = pyfits.getheader(finder)
    wcs = pywcs.WCS(header)
    ((ra,dec),) = wcs.all_pix2world([[250.1,416.5]],1)
    return ra,dec

def get_slit_coordinates(finder, slitlength):
    x = np.linspace(205.6,294.9,slitlength)
    y = np.linspace(415.7,416.7,slitlength)
    header = pyfits.getheader(finder)
    wcs = pywcs.WCS(header)
    (ra,dec) = wcs.all_pix2world(np.array([x,y]).T,1).T
    return ra,dec
