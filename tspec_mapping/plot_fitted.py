"""
Orion-specific test code...
"""
import aplpy
import os
import pyfits
import matplotlib.pyplot as pl
pl.ioff()

verbose_global = False

# orion-specific
ROrion = pyfits.getdata('/Volumes/disk5/Users/adam/observations/Q4CU05/UT121105/Robberto2010Orion.fits')
UOrion = pyfits.getdata('/Volumes/disk5/Users/adam/observations/Q4CU05/UT121105/UKIDSS_Orion.fits')
MOrion = pyfits.getdata('/Volumes/disk5/Users/adam/observations/tcam/UT121105/2mass_15arcmin.fits')

def plot_fitted(fn, plow=1, phi=99, verbose=verbose_global, 
        zoom_center=None, zoom_radius=1/60.):
    if verbose:
        print "Plotting %s" % fn
    prefix = fn.replace('.fits','')
    if not os.path.exists(prefix+'.corr.fits'):
        print "SKIP %s" % fn
        return
    F = aplpy.FITSFigure(prefix+'.corr.fits')
    F.show_grayscale(pmin=plow,pmax=phi)
    F.save(prefix+".corr.gray.png")

    # orion-specific...
    F.show_markers(ROrion['RAJ2000'], ROrion['DEJ2000'], marker='.', color=(0,1,0),
            facecolor=(0,1,0), edgecolor=(0,1,0), s=2)
    F.show_markers(UOrion['RA'], UOrion['Dec'], marker='.', color=(0,0,1),
            facecolor=(0,0,1), edgecolor=(0,0,1), s=2)
    F.show_markers(MOrion['RA'], MOrion['DEC'], marker='.', color=(0,1,1),
            facecolor=(0,1,1), edgecolor=(0,1,1), s=2)

    if os.path.exists(prefix+'.rdls'):
        rdls = pyfits.getdata(prefix+'.rdls')
        F.show_markers(rdls['RA'],rdls['DEC'],marker='x',s=5)
        F.save(prefix+".corr.rdls.png")
    else:
        F.save(prefix+".corr.no-rdls.png")

    if zoom_center:
        try:
            F.recenter(float(zoom_center[0]), float(zoom_center[1]),
                    float(zoom_radius))
            F.save(prefix+".corr.zoom.png")
        except:
            pass
    F.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose',default=False,action='store_true')
    parser.add_argument('--zoom-center',default=None,nargs=2)
    parser.add_argument('--zoom-radius',default=0.5/60.,nargs=1)
    parser.add_argument('files',nargs='*')
    args=parser.parse_args()
    print args

    verbose_global=args.verbose

    for fn in args.files:
        plot_fitted(fn,verbose=args.verbose, zoom_center=args.zoom_center,
                zoom_radius=args.zoom_radius)
