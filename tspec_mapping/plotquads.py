"""
Wrapper of plotquads to produce a halfway-reasonable .png image with quads overlaid
"""
import PIL.Image
import numpy as np
import pyfits
import astrometry
import os

verbose_global=False

def plot_quads(fn,plow=1,phi=99,verbose=verbose_global, color='red'):
    if verbose:
        print "Plotting %s" % fn
    prefix = fn.replace('.fits','')
    data = pyfits.getdata(fn)
    vmin = np.percentile(data,plow)
    vmax = np.percentile(data,phi)
    newimg = data
    newimg[data<vmin] = vmin
    newimg[data>vmax] = vmax
    newimg -= vmin
    newimg /= newimg.max() 
    newimg *= 255
    png = PIL.Image.fromarray(newimg.astype('int32'))
    lapng = png.convert('LA')
    lapng.save(prefix+"_grey.png")

    cmd = "plotquad -I {prefix}_grey.png -C {color} -c -m {prefix}.match > {prefix}.quads.png"
    cmd = cmd.format(prefix=prefix,color=color)
    if verbose:
        print cmd

    return os.popen(cmd).read()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose',default=False,action='store_true')
    parser.add_argument('--color',default='red')
    parser.add_argument('files',nargs='*')
    args=parser.parse_args()

    verbose_global=args.verbose

    for fn in args.files:
        stdout = plot_quads(fn,verbose=args.verbose,color=args.color)
        if args.verbose:
            print stdout

