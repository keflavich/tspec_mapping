import subprocess
import shlex
import os
import numpy as np

def _runcmd(cmd):
    """
    Execute a generic command string
    (wrapper for subprocess.Popen or whatever you want)
    """

    p = subprocess.Popen(shlex.split(cmd), shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

    stdout,stderr = p.communicate()

    return stdout.strip(),stderr.strip()

# found with solve-field | perl -ne 'm/--([a-z-]*)/; print "$1\n"' | uniq
def _solve_field_defaults():
    return {
    "help":None, "verbose":None, "dir":None, "out":None, "backend-config":None,
    "backend-batch":None, "files-on-stdin":None, "no-plots":None,
    "plot-scale":None, "plot-bg":None, "use-wget":None, "overwrite":None,
    "continue":None, "skip-solved":None, "new-fits":None, "kmz":None,
    "scamp":None, "scamp-config":None, "index-xyls":None, "just-augment":None,
    "no-delete-temp":None, "scale-low":None, "scale-high":None,
    "scale-units":None, "parity":None, "code-tolerance":None,
    "pixel-error":None, "quad-size-min":None, "quad-size-max":None,
    "odds-to-tune-up":None, "odds-to-solve":None, "odds-to-reject":None,
    "odds-to-stop-looking":None, "use-sextractor":None,
    "sextractor-config":None, "use-sextractor":None, "sextractor-path":None,
    "use-sextractor":None, "ra":None, "dec":None, "radius":None, "depth":None,
    "objs":None, "cpulimit":None, "resort":None, "extension":None,
    "no-fits":None, "invert":None, "downsample":None,
    "no-background-subtraction":None, "sigma":None, "no-remove-lines":None,
    "uniformize":None, "no-verify-uniformize":None, "no-verify-dedup":None,
    "no-fix-sdss":None, "cancel":None, "solved":None, "solved-in":None,
    "match":None, "rdls":None, "sort-rdls":None, "tag":None, "tag-all":None,
    "scamp-ref":None, "corr":None, "wcs":None, "pnm":None, "keep-xylist":None,
    "dont-augment":None, "verify":None, "no-verify":None, "guess-scale":None,
    "crpix-center":None, "crpix-x":None, "crpix-y":None, "no-tweak":None,
    "tweak-order":None, "temp-dir":None, "fields":None, "width":None,
    "height":None, "x-column":None, "y-column":None, "sort-column":None,
    "sort-ascending":None,         }

def solve_field(filename, **kwargs):
    """
    Solve a field.  Specify command line arguments with kwargs.  Short versions
    of command-line args are not supported. If an option has a dash in it, use
    underscore instead, e.g. instead of plot-scale=True (which is not valid
    python) use plot_scale=True.
    
    Program docs:
This program is part of the Astrometry.net suite.
For details, visit  http://astrometry.net .
Subversion URL svn+ssh://astrometry.net/svn/tags/tarball-0.38/astrometry/util/
Revision 16745, date 2010-11-19 20:47:53 +0000 (Fri, 19 Nov 2010).

Usage:   solve-field [options]  [<image-file-1> <image-file-2> ...] [<xyls-file-1> <xyls-file-2> ...]

You can specify http:// or ftp:// URLs instead of filenames.  The "wget" or "curl" program will be used to retrieve the URL.

Options include:
  -h / --help: print this help message
  -v / --verbose: be more chatty -- repeat for even more verboseness
  -D / --dir <directory>: place all output files in the specified directory
  -o / --out <base-filename>: name the output files with this base name
  -b / --backend-config <filename>: use this config file for the "backend"
          program
  --backend-batch: run backend once, rather than once per input file
  -f / --files-on-stdin: read filenames to solve on stdin, one per line
  -p / --no-plots: don't create any plots of the results
  --plot-scale <scale>: scale the plots by this factor (eg, 0.25)
  --plot-bg <filename (JPEG)>: set the background image to use for plots
  -G / --use-wget: use wget instead of curl
  -O / --overwrite: overwrite output files if they already exist
  -K / --continue: don't overwrite output files if they already exist; continue
          a previous run
  -J / --skip-solved: skip input files for which the 'solved' output file
          already exists; NOTE: this assumes single-field input files
  -N / --new-fits <filename>: output filename of the new FITS file containing
          the WCS header; "none" to not create this file
  -Z / --kmz <filename>: create KMZ file for Google Sky.  (requires wcs2kml)
  -i / --scamp <filename>: create image object catalog for SCAMP
  -n / --scamp-config <filename>: create SCAMP config file snippet
  -U / --index-xyls <filename>: output filename for xylist containing the image
          coordinate of stars from the index
  --just-augment: just write the augmented xylist files; don't run backend.
  -7 / --no-delete-temp: don't delete temp files (for debugging)

  -L / --scale-low <scale>: lower bound of image scale estimate
  -H / --scale-high <scale>: upper bound of image scale estimate
  -u / --scale-units <units>: in what units are the lower and upper bounds?
     choices:  "degwidth", "degw", "dw"   : width of the image, in degrees (default)
               "arcminwidth", "amw", "aw" : width of the image, in arcminutes
               "arcsecperpix", "app": arcseconds per pixel
  -8 / --parity <pos/neg>: only check for matches with positive/negative parity
          (default: try both)
  -c / --code-tolerance <distance>: matching distance for quads (default: 0.01)
  -E / --pixel-error <pixels>: for verification, size of pixel positional error
          (default: 1)
  -q / --quad-size-min <fraction>: minimum size of quads to try, as a fraction
          of the smaller image dimension, default: 0.1
  -Q / --quad-size-max <fraction>: maximum size of quads to try, as a fraction
          of the image hypotenuse, default 1.0
  --odds-to-tune-up <odds>: odds ratio at which to try tuning up a match that
          isn't good enough to solve (default: 1e6)
  --odds-to-solve <odds>: odds ratio at which to consider a field solved
          (default: 1e9)
  --odds-to-reject <odds>: odds ratio at which to reject a hypothesis (default:
          1e-100)
  --odds-to-stop-looking <odds>: odds ratio at which to stop adding stars when
          evaluating a hypothesis (default: HUGE_VAL)
  --use-sextractor: use SExtractor rather than built-in image2xy to find sources
  --sextractor-config <filename>: use the given SExtractor config file (default:
          etc/sextractor.conf).  Note that CATALOG_NAME and CATALOG_TYPE values
          will be over-ridden by command-line values.  This option implies
          --use-sextractor.
  --sextractor-path <filename>: use the given path to the SExtractor executable.
          Default: just 'sex', assumed to be in your PATH.  Note that you can
          give command-line args here too (but put them in quotes), eg:
          --sextractor-path 'sex -DETECT_TYPE CCD'.  This option implies
          --use-sextractor.
  -3 / --ra <degrees or hh:mm:ss>: only search in indexes within 'radius' of the
          field center given by 'ra' and 'dec'
  -4 / --dec <degrees or [+-]dd:mm:ss>: only search in indexes within 'radius'
          of the field center given by 'ra' and 'dec'
  -5 / --radius <degrees>: only search in indexes within 'radius' of the field
          center given by ('ra', 'dec')
  -d / --depth <number or range>: number of field objects to look at, or range
          of numbers; 1 is the brightest star, so "-d 10" or "-d 1-10" mean look
          at the top ten brightest stars only.
  --objs <int>: cut the source list to have this many items (after sorting, if
          applicable).
  -l / --cpulimit <seconds>: give up solving after the specified number of
          seconds of CPU time
  -r / --resort: sort the star brightnesses by background-subtracted flux; the
          default is to sort using acompromise between background-subtracted and
          non-background-subtracted flux
  -6 / --extension <int>: FITS extension to read image from.
  -2 / --no-fits2fits: don't sanitize FITS files; assume they're already valid
  --invert: invert the image (for black-on-white images)
  -z / --downsample <int>: downsample the image by factor <int> before running
          source extraction
  --no-background-subtraction: don't try to estimate a smoothly-varying sky
          background during source extraction.
  --sigma <float>: set the noise level in the image
  -9 / --no-remove-lines: don't remove horizontal and vertical overdensities of
          sources.
  --uniformize <int>: select sources uniformly using roughly this many boxes
          (0=disable; default 10)
  --no-verify-uniformize: don't uniformize the field stars during verification
  --no-verify-dedup: don't deduplicate the field stars during verification
  -0 / --no-fix-sdss: don't try to fix SDSS idR files.
  -C / --cancel <filename>: filename whose creation signals the process to stop
  -S / --solved <filename>: output file to mark that the solver succeeded
  -I / --solved-in <filename>: input filename for solved file
  -M / --match <filename>: output filename for match file
  -R / --rdls <filename>: output filename for RDLS file
  --sort-rdls <column>: sort the RDLS file by this column; default is ascending;
          use "-column" to sort "column" in descending order instead.
  --tag <column>: grab tag-along column from index into RDLS file
  --tag-all: grab all tag-along columns from index into RDLS file
  -j / --scamp-ref <filename>: output filename for SCAMP reference catalog
  -B / --corr <filename>: output filename for correspondences
  -W / --wcs <filename>: output filename for WCS file
  -P / --pnm <filename>: save the PNM file as <filename>
  -k / --keep-xylist <filename>: save the (unaugmented) xylist to <filename>
  -A / --dont-augment: quit after writing the unaugmented xylist
  -V / --verify <filename>: try to verify an existing WCS file
  -y / --no-verify: ignore existing WCS headers in FITS input images
  -g / --guess-scale: try to guess the image scale from the FITS headers
  --crpix-center: set the WCS reference point to the image center
  --crpix-x <pix>: set the WCS reference point to the given position
  --crpix-y <pix>: set the WCS reference point to the given position
  -T / --no-tweak: don't fine-tune WCS by computing a SIP polynomial
  -t / --tweak-order <int>: polynomial order of SIP WCS corrections
  -m / --temp-dir <dir>: where to put temp files, default /tmp
The following options are valid for xylist inputs only:
  -F / --fields <number or range>: the FITS extension(s) to solve, inclusive
  -w / --width <pixels>: specify the field width
  -e / --height <pixels>: specify the field height
  -X / --x-column <column-name>: the FITS column containing the X coordinate of
          the sources
  -Y / --y-column <column-name>: the FITS column containing the Y coordinate of
          the sources
  -s / --sort-column <column-name>: the FITS column that should be used to sort
          the sources
  -a / --sort-ascending: sort in ascending order (smallest first); default is
          descending order

Note that most output files can be disabled by setting the filename to "none".
 (If you have a sick sense of humour and you really want to name your output
  file "none", you can use "./none" instead.)
    """
    
    solve_field_args = _solve_field_defaults()
    for key,val in kwargs.iteritems():
        if key in solve_field_args.iteritems():
            solve_field_args[key] = val
        elif key.replace("_","-") in solve_field_args:
            solve_field_args[key.replace("_","-")] = val

    cmd = "solve_field "
    for key,val in solve_field_args:
        if val is True:
            cmd += "--%s " % key
        elif val is not None:
            cmd += "--%s %s " % (key,val)

    return _runcmd(cmd)


def _build_index_args(**kwargs):
    argkeys = {'sort_column':'S', 'scale_number':'P', 'nside':'N',
            'min_quad_size':'l', 'max_quad_size':'u',
            'reverse_sortorder':'f', 'healpix_nside':'U', 'big_healpix':'H',
            'big_healpix_nside':'s', 'margin':'m', 'sweeps':'n',
            'dedup_radius':'r', 'jitter_arcsec':'j', 'dimquads':'d',
            'passes':'p', 'reuse_times':'R', 'max_reuses':'L',
            'scan_catalog':'E', 'unique_id':'I', 'in_memory':'M',
            'keep_temp':'T', 'verbose':'v', }
    defaults = dict([(k,None) for k in argkeys.keys()])
    defaults['reverse_sortorder'] = False
    defaults['big_healpix_nside'] = 1
    defaults['sweeps'] = 10
    defaults['margin'] = 0
    defaults['jitter_arcsec'] = 1
    defaults['dimquads'] = 4
    defaults['passes'] = 16
    defaults['reuse_times'] = 8
    defaults['in_memory'] = False
    defaults['keep_temp'] = False
    defaults['verbose'] = False

    args = {}

    for k,v in kwargs.iteritems():
        if k not in argkeys.keys():
            raise ValueError("%s is not a valid key" % k)
        else:
            charkey = argkeys[k]
            if k in defaults.keys():
                # only set arg if it's not the default
                if defaults[k] != v:
                    args[charkey] = v
            else: # those with no defaults get set 
                args[charkey] = v

    return args



def build_index(infile, outfile=None, index_dir=None, debug=False, **kwargs):
    """
Usage: build-index
      (
         -i <input-FITS-catalog>  input: source RA,DEC, etc
    OR,
         -1 <input-index>         to share another index's stars
      )
      -o <output-index>        output filename for index
      (
         -P <scale-number>: use 'preset' values for '-N', '-l', and '-u'
               (the scale-number is the last two digits of the pre-cooked
                index filename -- eg, index-205 is  "-P 5".
                -P 0  should be good for images about 6 arcmin in size
                    and it goes in steps of sqrt(2), so:
                -P 2  should work for images about 12 arcmin across
                -P 4  should work for images about 24 arcmin across
                -P 6  should work for images about 1 degree across
                -P 8  should work for images about 2 degree across
                -P 10 should work for images about 4 degree across
                 etc... up to -P 19
  OR,
         -N <nside>            healpix Nside for quad-building
         -l <min-quad-size>    minimum quad size (arcminutes)
         -u <max-quad-size>    maximum quad size (arcminutes)
      )
      [-S]: sort column (default: assume the input file is already sorted)
      [-f]: sort in descending order (eg, for FLUX); default ascending (eg, for MAG)
      [-U]: healpix Nside for uniformization (default: same as -n)
      [-H <big healpix>]; default is all-sky
      [-s <big healpix Nside>]; default is 1
      [-m <margin>]: add a margin of <margin> healpixels; default 0
      [-n <sweeps>]    (ie, number of stars per fine healpix grid cell); default 10
      [-r <dedup-radius>]: deduplication radius in arcseconds; default no deduplication
      [-j <jitter-arcsec>]: positional error of stars in the reference catalog (in arcsec; default 1)
      [-d <dimquads>] number of stars in a "quad" (default 4).
      [-p <passes>]   number of rounds of quad-building (ie, # quads per healpix cell, default 16)
      [-R <reuse-times>] number of times a star can be used (default: 8)
      [-L <max-reuses>] make extra passes through the healpixes, increasing the "-r" reuse
                     limit each time, up to "max-reuses".
      [-E]: scan through the catalog, checking which healpixes are occupied.
      [-I <unique-id>] set the unique ID of this index
      [-M]: in-memory (don't use temp files)
      [-T]: don't delete temp files
      [-v]: add verbosity.
    """        
    args = _build_index_args(**kwargs)

    build_index = os.popen('which build-index').read()
    cmd = "%s -i %s " % (build_index,infile)

    outdir = _get_index_dir()
    if outfile is not None:
        cmd += "-o %s " % (outdir+outfile)
    else:
        if "P" in args:
            outfile = infile+"_P%s.index" % args["P"]
        elif "N" in args:
            outfile = infile+"_N%s.index" % args["N"]
        else:
            outfile = infile+".index"
        cmd += "-o %s " % (outdir+outfile)

    for key,val in args.iteritems():
        if val is True:
            # booleans just get the -A arg added...
            cmd += "-%s " % key
        elif val is not None:
            cmd += "-%s %s " % (key,val)

    if debug:
        print cmd

    return _runcmd(cmd)

def _get_index_dir():
    """
    Use the path to the build-index executable to determine the path to the
    astrometry.net build indices
    """
    build_index_path = os.path.split(os.popen('which build-index').read().strip())[0]
    data_path = os.path.split(build_index_path)[0]+"/data/"
    if not os.access(data_path,os.W_OK):
        raise IOError("Permissions in output directory %s do not allow writing" % data_path)
    return data_path


def get_closest_preset(fieldsize):
    """
    Given a field size in arcminutes, return the preset with the closest field
    size
    """
    # from build-index-main.c
    scales = [0.35, 0.5, 0.7, 1., 1.4,
              2., 2.8, 4., 5.6, 8., 11., 16., 22., 30., 42., 60., 85.,
              120., 170., 240., 340., 480., 680., 1000., 1400., 2000. ]
    scale_numbers = [x-5 for x in range(len(scales))]

    closest_fieldsize = np.argmin(np.abs(np.array(scales)-fieldsize))
    closest_scale = scale_numbers[closest_fieldsize]
    return closest_scale

