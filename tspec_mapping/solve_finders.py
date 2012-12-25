# Step 1: Grab 2MASS catalog
#
"""
import astroquery.irsa
ORION = astroquery.irsa.query_gator_box('pt_src_cat','83.808 -5.391',900)
import atpy
import astropy
T = atpy.Table()
[T.add_column(n,ORION[n]) for n in ORION.columns]
T.write('2MASS_orion_cat.fits')
T2 = atpy.Table()
T2.add_column('RA',T.ra)
T2.add_column('Dec',T.dec)
T2.add_column('Kmag',T.k_m)
T2.write('2MASS_orion_cat_simple.fits')
"""

#    query = {}
#    query["-source"] = "II/314/gcs8" # UKIDSS DR8
#    query["-out"] = ["RAJ2000","DEJ2000","Kmag"]
#    query["-c"] = "05 35 14 -05 22.4"
#    query["-c.u"] = 'arcmin'
#    query["-c.r"] = "25"
#    query["-c.geom"] = 'r'
#    table1 = vizquery(query)
    


# Step 2: build indices
#for P in xrange(-4,3):
#    #cmd = 'build-index -S Kmag -E -v  -P {0} -i 2MASS_orion_cat_simple.fits -o 2MASS_orion_simple_index_P{0}'.format(P)
#    #cmd = 'build-index -S Kmag -E -v  -P {0} -i Robberto2010Orion_fixed.fits -o Robberto_orion_index_P{0}'.format(P)
#    cmd = 'build-index -S Kmag -E -v  -P {0} -i UKIDSS_Orion.fits -o UKIDSS_Orion_index_P{0}'.format(P)
#    print cmd
#    #os.system(cmd)

import multiprocessing
import glob
import os
import subprocess
import shlex
import shutil
import pyfits
import numpy as np
import time

mask = pyfits.getdata('mask0001.fits')

verbose_global = False
resolve_global = False

def solvefield(fn,verbose=verbose_global,resolve=resolve_global):
    prefix = fn.replace('.fits','')
    if os.path.exists(prefix+'.solved') and not resolve:
        if verbose:
            print "Solved %s" % prefix
    else:
        cmd = 'solve-field -g --crpix-center --no-plots --overwrite --continue %s.fits' % prefix
        if not os.path.exists("%s.fits" % prefix):
            raise IOError("CANT FIND %s.fits" % prefix)
        if verbose:
            print cmd
        p = subprocess.Popen(shlex.split(cmd), shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        #while p.poll() is None:
        #    time.sleep(0.1)
        #    print ".",
        stdout,stderr = p.communicate()
        #if stderr:
        #    raise Exception(stderr)
        cmd2 = "solve-field --verify {0}.wcs --no-plots --continue {0}.fits".format(prefix)
        p = subprocess.Popen(shlex.split(cmd2), shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        stdout2,stderr2 = p.communicate()
        return stdout.strip(),stderr.strip(),stdout2.strip(),stderr2.strip()

import plot_axy
# def plot_axy(fn):
#     #prefix = fn.replace('.fits','')
#     cmd2 = "python plot_axy.py %s" % fn
#     p = subprocess.Popen(shlex.split(cmd2), shell=False,
#             stdout=subprocess.PIPE,
#             stderr=subprocess.PIPE)
#     stdout2,stderr2 = p.communicate()
#     return stdout2.strip(),stderr2.strip()

def mask_field(fn,verbose=verbose_global):
    prefix = fn.replace('.fits','')
    # mask file
    if os.path.exists(prefix+'.new'):
        if verbose:
            print "Masking %s" % prefix
        try:
            f = pyfits.open(prefix+'.new')
            f[0].data[mask<0.75] = np.nan
            f.writeto(prefix+'.new', clobber=True, output_verify='fix')
        except Exception as ex:
            print ex
    else:
        print "FAILED TO SOLVE & MASK %s" % prefix
        return 'FAILURE'

def deproject(fn,verbose=verbose_global):
    prefix = fn.replace('.fits','')
    if verbose:
        print "Deproject ",fn
    # deproject file
    #cp $fn $fn.corr
    shutil.copy(prefix+'.wcs',prefix+'.wcs.corr')
    # modhead $fn.corr CTYPE1 RA---TAN
    cmd = "modhead %s.wcs.corr CTYPE1 RA---TAN" % prefix
    p = subprocess.Popen(shlex.split(cmd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    #p.check_call()
    # stderr = p.stderr.read()
    # stdout = p.stderr.read()
    # modhead $fn.corr CTYPE2 DEC--TAN
    stdout,stderr = p.communicate()
    cmd = "modhead %s.wcs.corr CTYPE2 DEC--TAN" % prefix
    p = subprocess.Popen(shlex.split(cmd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    stdout,stderr = p.communicate()
    # p.check_call()
    # stderr = p.stderr.read()
    # stdout = p.stderr.read()
    # wcs-resample -w $fn ${fn%wcs}new $fn.corr ${fn%wcs}corr.fits' sh
    cmd = "wcs-resample -w {0}.wcs {0}.new {0}.wcs.corr {0}.corr.fits".format(prefix)
    if verbose:
        print cmd
    p = subprocess.Popen(shlex.split(cmd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    stdout,stderr = p.communicate()
    return stdout

def fix_header_time(fn,verbose=verbose_global):
    if verbose:
        print "Fix header ",fn
    prefix = fn.replace('.fits','')
    header = pyfits.getheader(fn)
    f = pyfits.open(prefix+".corr.fits")
    if not 'DATE-OBS' in header:
        raise IOError("%s is apparently invalid" % fn)
    else:
        DO = header['DATE-OBS']
    try:
        f[0].header.add_blank('DATE-OBS')
        f[0].header['DATE-OBS'] = DO
    except Exception as ex:
        print "pyfits failure",ex
    f.writeto(prefix+".corr.fits",clobber=True,output_verify='fix')
    # Fix header TIME
    # cmd = 'sethead {0}.corr.fits DATE-OBS=`gethead DATE-OBS {0}.fits`'.format(prefix)
    # if verbose:
    #     print cmd
    # p = subprocess.check_call(shlex.split(cmd),
    #         stdout=subprocess.PIPE,
    #         stderr=subprocess.PIPE)
    return 

def process_file(fn, verbose=verbose_global, resolve=resolve_global):
    #prefix = fn.replace('.fits','')
    print os.getpid(),fn

    try:
        sf = solvefield(fn, resolve=resolve, verbose=verbose)
        if verbose:
            print sf
        #P = multiprocessing.Process(target=plot_axy.plot_axy,args=(fn,))
        #P.start()
        mask_field(fn,verbose=verbose)
        deproject(fn,verbose=verbose)
        fix_header_time(fn,verbose=verbose)
    except Exception as ex:
        print "FILE %s FAILED" % fn
        print ex

    print "SUCCESS? ",fn,os.getpid()
    return "SUCCESS"
     
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--nprocs',metavar='n',type=int,default=None)
    parser.add_argument('--verbose',default=False,action='store_true')
    parser.add_argument('--resolve',default=False,action='store_true')
    parser.add_argument('--glob',metavar='g',type=str,default='proc-t[0-9][0-9][0-9][0-9].fits')
    args=parser.parse_args()
    nprocs = int(args.nprocs)

    verbose_global = args.verbose
    resolve_global = args.resolve


    print nprocs, args.glob

    if nprocs is None:
        p = multiprocessing.Pool()
    elif nprocs > 1:
        p = multiprocessing.Pool(nprocs)
    else:
        p = None

    files = glob.glob(args.glob)

    if args.resolve:
        for fn in files:
            solved = fn.replace(".fits",".solved")
            if "solved" in solved:
                if os.path.exists(solved):
                    print "Removing %s" % solved
                    os.remove(solved)
            else:
                raise IOError("file name was wrong: %s" % fn)



    import time
    t0 = time.time()

    x = p.map_async(process_file,files)
    print "After map_async",time.time()-t0
    print x.get()

    # for fn in files:
    #     if args.verbose:
    #         print "Applying to %s" % (fn)
    #     if p is not None:

    #         RETURN = p.apply(process_file,args=(fn,),kwds={'verbose':args.verbose})
    #         print "returned " ,RETURN," dt = %0.1f" % (time.time()-t0)
    #     else:
    #         RETURN = process_file(fn,verbose=args.verbose)
    #         print "returned " ,RETURN," dt = %0.1f" % (time.time()-t0)

    #p.map(process_file,files,verbose=args.verbose)

# Step 3: Solve fields

# BASH:
# for fn in `ls proc-t[0-9][0-9][0-9][0-9].fits`; do solve-field -g --crpix-center $fn; done
# find ./ -name 'proc-t[0-9][0-9][0-9][0-9].fits' -print0 | xargs -0 -n 1 -P 6 sh -c 'fn=$1; if [ ! -e ${fn%fits}solved ]; then solve-field -g --crpix-center --overwrite --continue $fn; else echo ${fn%fits}solved; fi' sh
# 

# Step 4: mask fields
# python mask.py  NOT THIS: python -c "import glob,pyfits; fl = glob.glob('proc-t[0-9][0-9][0-9][0-9].new'); for fn in fl:
#
# Step 5: Deproject fields (optional)
# parallel:
# find ./ -name '*.wcs' -print0 | xargs -0 -n 1 -P 12 sh -c 'fn=$1; cp $fn $fn.corr; modhead $fn.corr CTYPE1 RA---TAN; modhead $fn.corr CTYPE2 DEC--TAN; wcs-resample -w $fn ${fn%wcs}new $fn.corr ${fn%wcs}corr.fits' sh
# single thread:
# for fn in `ls *.wcs`; do cp $fn $fn.corr; modhead $fn.corr CTYPE1 RA---TAN; modhead $fn.corr CTYPE2 DEC--TAN; wcs-resample -w $fn ${fn%wcs}new $fn.corr ${fn%wcs}corr.fits; done
# 
#
# Step 6: montage fields
# montage proc-t*corr.fits --outfile=orion_finder_mosaic.fits --header=header.hdr --combine=median --tmpdir=tmp &
#
# 
