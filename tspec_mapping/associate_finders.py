import datetime
import os
import numpy as np
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import glob
try: 
    import progressbar
    widgets = [progressbar.FormatLabel('Processed: %(value)d spectra in %(elapsed)s)'), progressbar.Percentage()]
    progress = progressbar.ProgressBar(widgets=widgets)
except ImportError:
    def progress(x):
        yield x

def find_nearest_finder(time, findertimes):
    """
    Given a time and a list of times, returns the index of the time closest to
    the requested time
    """

    deltat = [(t1-time).seconds for t1 in findertimes] 
    #print [(t1,time,(t1-time).seconds) for t1 in findertimes]

    nearest = np.argmin(np.abs(deltat))

    return nearest

def find_finders_between(starttime, endtime, findertimes):
    """
    return indices of finders with acceptable times
    (I don't think numpy arrays can do time comparison...)
    """
    bools =  np.array([((T>starttime) and (T<endtime)) for T in findertimes])
    indices = np.nonzero(bools)[0]
    return indices


def get_times(filelist):
    """
    Extract a list of datetime.datetime objects from a series of fits files
    with DATE-OBS header keywords
    """
    times = []
    for fn in filelist:
        hdr = pyfits.getheader(fn)
        timestr = hdr.get('DATE-OBS')
        times.append(datetime.datetime.strptime(timestr[:-4],'%Y-%m-%dT%H:%M:%S'))
    return times

def find_finders(spectrum_list, finder_list, do_average=True):

    finder_times = get_times(finder_list)

    finder_dict = {}

    for spec in progress(spectrum_list):
        hdr = pyfits.getheader(spec)
        timestr = hdr.get('DATE-OBS')
        starttime = datetime.datetime.strptime(timestr[:-4],'%Y-%m-%dT%H:%M:%S')
        endtime = starttime + datetime.timedelta(seconds=hdr['EXPTIME'])
        if do_average:
            finder_indices = find_finders_between(starttime,endtime,finder_times)
            finders = [finder_list[ii] for ii in finder_indices]
            if len(finders) == 0:
                finder_dict[spec] = finder_list[find_nearest_finder(starttime, finder_times)]
            else:
                combined_finder = spec.replace(".fits","_finder.fits")
                average_images(finders, combined_finder)
                finder_dict[spec] = combined_finder
        else:
            finder = finder_list[find_nearest_finder(starttime, finder_times)]
            finder_dict[spec] = finder

    return finder_dict

def average_images(image_list, outfilename, combine=np.median, clobber=True):
    header = pyfits.getheader(image_list[0])
    data_list = [pyfits.getdata(fn) for fn in image_list]
    new_data = combine(np.array(data_list), axis=0)
    newHDU = pyfits.PrimaryHDU(data=new_data, header=header)
    newHDU.writeto(outfilename,clobber=clobber)



def find_finders_dirs(spec_dir, finder_dir, includeglob="*", specsuffix=".fits", findersuffix='.fits'):

    spglob = spec_dir+includeglob+".[0-9][0-9][0-9][0-9]"+specsuffix
    print spglob

    spectra = glob.glob(spglob)

    finderglob = finder_dir+'proc-t[0-9][0-9][0-9][0-9]'+findersuffix
    print finderglob

    finders = glob.glob(finderglob)

    if len(spectra) == 0:
        raise ValueError("No spectra found")
    if len(finders) == 0:
        raise ValueError("No finders found")

    return find_finders(spectra, finders)


if __name__ == "__main__":
    import optparse
    parser=optparse.OptionParser()
    #parser.add_option("--band",default='K',help="2MASS band name")
    #parser.add_option("--usetext",default=False,help="Use the text attribute of the region file to name the finder?",action='store_true')
    parser.set_description("""Associate each spectrum with the first finder taken after the spectrum exposure was started.
    Arguments are spectrum directory, guider directory, e.g.:
    associate_finders('/path/to/Q2UV01/UT123456','/path/to/tcam/UT123456')

    An optional 3rd argument is an 'includeglob', defaulting to *
            """)
    options,args = parser.parse_args()

    matches = find_finders_dirs(*args)

    longest_key = 0
    longest_value = 0
    for key,value in matches.iteritems():
        longest_key = len(key) if len(key) > longest_key else longest_key
        if hasattr(value,'__len__'):
            longest_value = len(value) if len(value) > longest_value else longest_value
        printstr = "%%-%is : %%-%is" % (longest_key, longest_value)
    for key,value in matches.iteritems():
        print printstr % (key,value)
        f = pyfits.open(key)
        if 'FINDER' in f[0].header: del f[0].header['FINDER']
        if 'FINDER_1' in f[0].header: del f[0].header['FINDER_1']
        if 'FINDER_2' in f[0].header: del f[0].header['FINDER_2']
        f[0].header['FINDER'] = str(value)
        f.writeto(key,clobber=True,output_verify='fix')
        #os.system('sethead %s FINDER=%s' % (key,value))
