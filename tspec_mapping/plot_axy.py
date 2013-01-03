"""
Essentially a replacement of the plot_axy code in astrometry.net because mine
didn't work
"""
import matplotlib.pyplot as pl
pl.ioff()
import numpy as np
import pyfits

verbose_global=False

def plot_axy(fn,plow=1,phi=99,verbose=verbose_global):
    if verbose:
        print "Plotting %s" % fn
    prefix = fn.replace('.fits','')
    axy = pyfits.getdata(prefix+".axy")
    data = pyfits.getdata(fn)
    vmin = np.percentile(data,plow)
    vmax = np.percentile(data,phi)
    pl.clf()
    pl.imshow(data,cmap='gray',vmin=vmin,vmax=vmax)
    pl.axis([0,data.shape[1],0,data.shape[0]])
    pl.savefig(prefix+"_axy_nomarkers.png",bbox_inches='tight')
    pl.plot(axy['X']-1,axy['Y']-1,'gx',markeredgewidth=1)
    try:
        xyls = pyfits.getdata(prefix+"-indx.xyls")
        pl.plot(xyls['X']-1,xyls['Y']-1,'r+')
    except IOError:
        print "No XYLS for file %s" % prefix
    pl.axis([0,data.shape[1],0,data.shape[0]])
    pl.savefig(prefix+"_axy.png",bbox_inches='tight')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose',default=False,action='store_true')
    parser.add_argument('files',nargs='*')
    args=parser.parse_args()

    verbose_global=args.verbose

    for fn in args.files:
        plot_axy(fn,verbose=args.verbose)
