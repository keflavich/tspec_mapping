"""
Given a region file, create finder images 5x5' of the field in K-band
"""
import pyregion
import subprocess
import coords
import shutil
import os

def make_finders(regfile,usetext=False, band='K'):
    reg = pyregion.open(regfile)
    for R in reg:
        if R.coord_format == 'galactic':
            ra,dec = coords.Position(R.coord_list[:2],system='galactic').j2000()
        else:
            ra,dec = R.coord_list[:2]

        if usetext and 'text' in R.attr[1]: 
            filename = R.attr[1]['text'].replace(" ","_")+'_%s_finder.fits' % (band)
        else:
            filename = "%06.3f%+06.3f_%s_finder.fits" % (ra,dec,band)

        print "ra,dec: %06.3f %+06.3f" % (ra,dec)," coords:", R.coord_list

        make_finder(ra,dec, outfilename=filename, band=band)

def make_finder(ra, dec, outfilename=None, band='K', size_arcmin=7):
    """
    Use montage to create a 2MASS finder chart from the specified ra,dec
    coordinates.  

    Parameters
    ----------
    ra,dec : float,float
        The RA, Dec coordinates to center the finder on
    outfilename : string
        A FITS file name
    band : [ 'J', 'H', 'K' ]
        The 2MASS band to use
    size_arcmin : float
        The size of the field to observe in arcminutes (square)
    """

    make_header(ra, dec, size_arcmin=size_arcmin)

    prefix = "%06.3f%+06.3f" % (ra,dec)
    command = "mExec -d 1 -k -o %s -f %s.hdr 2MASS %s %s/" % (outfilename,prefix,band,prefix)

    if os.path.exists(prefix):
        shutil.rmtree(prefix)

    print command
    sp = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)

    sp.wait()
    print sp.returncode, "".join(sp.stdout.readlines())
    #shutil.move("%s/mosaic.fits" % prefix,"%s_%s_finder.fits" % (prefix,band))

    shutil.rmtree(prefix)
    #    raise IOError("Failed to create %s from %s" % (outfilename,prefix))

def make_header(ra, dec, size_arcmin=7):
    prefix = "%06.3f%+06.3f" % (ra,dec)

    cdelt = 0.0002
    size_pix = size_arcmin / 60. / cdelt

    outf = open(prefix+".hdr",'w')

    print >>outf,"SIMPLE  =                    T"
    print >>outf,"BITPIX  =                  -64"
    print >>outf,"NAXIS   =                    2   / # of Axes"
    print >>outf,"NAXIS1  =                  %i" % (size_pix)
    print >>outf,"NAXIS2  =                  %i" % (size_pix)
    print >>outf,"CTYPE1  = 'RA---TAN'           / Orthographic Projection"
    print >>outf,"CTYPE2  = 'DEC--TAN'           / Orthographic Projection"
    print >>outf,"CRPIX1  =                  %i /   Axis 1 Reference Pixel" % (size_pix/2)
    print >>outf,"CRPIX2  =                  %i /   Axis 2 Reference Pixel" % (size_pix/2)
    print >>outf,"CRVAL1  =         %f /   RA  at Frame Center, J2000 (deg)" % (ra)
    print >>outf,"CRVAL2  =         %f /   Dec at Frame Center, J2000 (deg)" % (dec)
    print >>outf,"CROTA2  =        0.00000000000 /   Image Twist +AXIS2 W of N, J2000 (deg)"
    print >>outf,"CDELT1  =     -0.0002000000000 /   Axis 1 Pixel Size (degs)"
    print >>outf,"CDELT2  =      0.0002000000000 /   Axis 2 Pixel Size (degs)"
    
    outf.close()

    print "Created header file %s.hdr" % prefix

if __name__ == "__main__":
    import optparse
    parser=optparse.OptionParser()
    parser.add_option("--band",default='K',help="2MASS band name")
    parser.add_option("--usetext",default=False,help="Use the text attribute of the region file to name the finder?",action='store_true')
    options,args = parser.parse_args()

    make_finders(args[0], usetext=options.usetext, band=options.band)
