import numpy as np 
try:
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
except ImportError:
    import pyfits
    import pywcs
from astropy import coordinates
import coords
import AG_image_tools
import agpy
import montage
import numpy as np
from agpy import cross_correlation
np.seterr(all='ignore')
import sys
try:
    sys.path.append('/Users/adam/repos/aposoftware/reduction/')
    from get_slit_coords import get_slit_coordinates,get_slit_center
    from associate_finders import find_finders
except Exception as ex:
    print ex
    #matches=

def make_cube(filelist, outfilename, clobber=True, brgfile=None,
        stellarfile=True, h2file=True, brdfile=True, brefile=True,
        reference_image=None, reference_radec=None, reference_pix=None,
        pabfile=True, rot90=False, cd1=0.5, cd2=0.5,):
    """
    Create a data cube from a list of reduced TripleSpec spectra

    Parameters
    ----------
    filelist : list
        A list of the files to include.  Generate by, e.g.,
        filelist = glob.glob("*JHK_cal.fits")
    outfilename : string
        The output FITS file name
    clobber : bool
        Overwrite the output file if it exists?
    brgfile : string
        Create a Brackett Gamma integrated image if this is set to a FITS file
        name
    stellarfile : string
        Create a stellar continuum image from a relatively clean part of the
        spectrum if this is set to a FITS file name
    h2file : string
        Create an H2 2.12um integrated image if this is set to a FITS file name
    reference_image : string
        A FITS file to use as the pointing standard for cross-correlation based
        pointing correction.  Makes use of agpy's cross_correlation_shift tool
    reference_radec : (float,float)
        The RA/Dec (CRVAL1/CRVAL2) of the image reference position.  If not
        specified, will be extracted from one of the input spectra (likely
        incorrectly)
    reference_pix : (float,float)
        The X,Y pixel (CRPIX1/CRPIX2) reference in the image.  If not
        specified, will default to the center.
    rot90 : bool
        Should the spectra be rotated by 90 degrees before stacking?
    cd1 : float
    cd2 : float
        The CD1_1 and CD2_2 / CDELT1 & CDELT2 parameters for the FITS header
    """

    f0 = pyfits.open(filelist[0])

    cubeJHK = np.zeros([f0[0].header['NAXIS2']+20,f0[0].header['NAXIS1'],len(filelist)+20])
    print 'cubeshape (no adds): ',[f0[0].header['NAXIS2'],f0[0].header['NAXIS1'],len(filelist)]
    print 'cubeshape (with adds): ',[f0[0].header['NAXIS2']+20,f0[0].header['NAXIS1'],len(filelist)+20]


    fout = pyfits.open(filelist[0])
    print "cube shape: ",cubeJHK.shape
    fout[0].data = cubeJHK = cubeJHK.swapaxes(0,1)
    if rot90:
        fout[0].data = cubeJHK = fout[0].data.swapaxes(1,2)[:,::-1,:]
    fout[0].header.update('NAXIS3',cubeJHK.shape[0])
    fout[0].header.update('NAXIS',3)
    print "fout.data shape: ",fout[0].data.shape

    if reference_radec is not None:
        ra,dec = reference_radec
    else:
        fcen = pyfits.open(filelist[len(filelist)/2])
        ra = fcen[0].header.get('RA')
        dec = fcen[0].header.get('DEC')
        #ra,dec = coords.Position(ra+" "+dec).j2000()
        C = coordinates.ICRS(ra,dec,unit=('deg','deg'))
        ra,dec = C.ra.deg,C.dec.deg

    if reference_pix is not None:
        xpix,ypix = reference_pix
    else:
        xpix = fout[0].data.shape[2]/2.+1
        ypix = fout[0].data.shape[1]

    fout[0].header.update('CRVAL3', fout[0].header['CRVAL1'])
    fout[0].header.update('CRPIX3', fout[0].header['CRPIX1'])
    fout[0].header.update('CD3_3', fout[0].header['CDELT1'])
    fout[0].header.update('CRVAL1', ra)
    fout[0].header.update('CRVAL2', dec)
    fout[0].header.update('CRPIX1', xpix)
    fout[0].header.update('CRPIX2', ypix)
    fout[0].header.update('CD1_1', -cd1/3600.)
    fout[0].header.update('CD2_2',  cd2/3600.)
    fout[0].header.update('CTYPE1', 'RA---TAN')
    fout[0].header.update('CTYPE2', 'DEC--TAN')
    fout[0].header.update('CTYPE3', 'LINEAR')
    del fout[0].header['CDELT1']
    del fout[0].header['CDELT2']

    outhd = fout[0].header
    flathead = agpy.cubes.flatten_header(outhd)
    outwcs = pywcs.WCS(flathead)

    wtcube = cubeJHK*0

    for ii,fn in enumerate(filelist): 
        d = pyfits.getdata(fn)
        h = pyfits.getheader(fn)
        findername = h.get('FINDER')
        if findername is None or findername == '':
            raise IOError("No finder found!")
        slitcoords = get_slit_coordinates(findername,d.shape[0]) 
        pixcoords = outwcs.wcs_world2pix(np.array(zip(*slitcoords)),0)     
        print ii,fn,len(slitcoords),len(pixcoords),len(pixcoords[0])
        for x,y,di in zip(pixcoords[:,0],pixcoords[:,1],d):                  
            xrd = round(x)
            yrd = round(y)
            xrem = x-xrd
            yrem = y-yrd
            weights = ((0.5+(abs(xrem)-0.5)*np.sign(xrem)) * (0.5+(abs(yrem)-0.5)*np.sign(yrem)),
                       (0.5+(abs(xrem)-0.5)*np.sign(xrem)) * (0.5-(abs(yrem)-0.5)*np.sign(yrem)),
                       (0.5-(abs(xrem)-0.5)*np.sign(xrem)) * (0.5-(abs(yrem)-0.5)*np.sign(yrem)),
                       (0.5-(abs(xrem)-0.5)*np.sign(xrem)) * (0.5+(abs(yrem)-0.5)*np.sign(yrem)))
            if np.any(np.array(weights) < 0):
                print weights
                raise ValueError('negative weight')
            wtsum = sum(weights)
            if (wtsum - 1) > 1e-5:
                raise ValueError("wtsum = 1+%g" % (wtsum-1))
            cubeJHK[:,xrd,yrd]     += di * weights[0]
            cubeJHK[:,xrd,yrd-1]   += di * weights[1]
            cubeJHK[:,xrd-1,yrd-1] += di * weights[3]
            cubeJHK[:,xrd-1,yrd]   += di * weights[2]
            wtcube[:,xrd,yrd]     += weights[0]
            wtcube[:,xrd,yrd-1]   += weights[1]
            wtcube[:,xrd-1,yrd-1] += weights[3]
            wtcube[:,xrd-1,yrd]   += weights[2]
            #wtcube[:,x,y] += (1-xrem) * (1-yrem)

    print "wtcube shape: ",wtcube.shape,"  wtcube.sum(): ",wtcube.sum()," single plane: ",wtcube[0,:,:].sum()
    print "max wt: ",wtcube[0,:,:].max()
    print "min wt: ",wtcube[0,:,:].min()

    fout[0].data = cubeJHK/wtcube

    fout.writeto(outfilename,clobber=clobber,output_verify='fix')
    #fout[0].data = wtcube
    #fout.writeto(outfilename.replace('cube','weightcube'),clobber=clobber,output_verify='fix')

    cubeJHK = cubeJHK/wtcube

    if brgfile is not None:
        if brgfile is True:
            if 'cube.fits' in outfilename:
                brgfile = outfilename.replace("cube.fits","brg.fits")
            else:
                raise 
        brg = cubeJHK[3582:3590,:,:].sum(axis=0) #- np.median(cubeJHK[3592:3605,:,:],axis=0)*8
        fout[0].data=brg
        del fout[0].header['CTYPE3']
        del fout[0].header['CRVAL3']
        del fout[0].header['CRPIX3']
        del fout[0].header['CD3_3']
        fout.writeto(brgfile,clobber=clobber,output_verify='fix')

        if reference_image is not None:
            offsets = AG_image_tools.cross_correlation_shifts_FITS(brgfile, reference_image)
            fout[0].header['CRPIX1'] -= offsets[0][0]
            fout[0].header['CRPIX2'] -= offsets[1][0]
            fout.writeto(brgfile,clobber=clobber,output_verify='fix')

            foutcube = pyfits.open(outfilename)
            foutcube[0].header['CRPIX1'] -= offsets[0][0]
            foutcube[0].header['CRPIX2'] -= offsets[1][0]
            foutcube.writeto(outfilename,clobber=clobber,output_verify='fix')

    if stellarfile is not None:
        if stellarfile is True:
            if 'cube.fits' in outfilename:
                stellarfile = outfilename.replace("cube.fits","brg_cont.fits")
            else:
                raise 
        stellar = (np.sum(cubeJHK[3570:3575,:,:],axis=0) + np.sum(cubeJHK[3592:3594,:,:],axis=0) + np.sum(cubeJHK[3598:3601,:,:],axis=0)) * 8./10.
        fout[0].data=stellar
        fout.writeto(stellarfile,clobber=clobber,output_verify='fix')

    if h2file is not None:
        if h2file is True:
            if 'cube.fits' in outfilename:
                h2file = outfilename.replace("cube.fits","h2.fits")
            else:
                raise 
        h2 = np.sum(cubeJHK[3431:3435,:,:],axis=0)# - np.median(cubeJHK[3435:3445,:,:],axis=0)*4
        fout[0].data=h2
        fout.writeto(h2file,clobber=clobber,output_verify='fix')

        if h2file is True:
            if 'cube.fits' in outfilename:
                h2cfile = outfilename.replace("cube.fits","h2_cont.fits")
            else:
                raise 
        else:
            if 'h2' in h2file:
                h2cfile = h2file.replace("h2","h2_cont")
            else:
                raise
        h2c = (np.sum(cubeJHK[3446:3452,:,:],axis=0))*4./6.
        fout[0].data=h2c
        fout.writeto(h2cfile,clobber=clobber,output_verify='fix')

    if brdfile is not None:
        if brdfile is True:
            if 'cube.fits' in outfilename:
                brdfile = outfilename.replace("cube.fits","brd.fits")
            else:
                raise 
        brd = np.sum(cubeJHK[2819:2824,:,:],axis=0) #- np.median(cubeJHK[2825:2835,:,:],axis=0)*5
        fout[0].data=brd
        fout.writeto(brdfile,clobber=clobber,output_verify='fix')

        if brdfile is True:
            if 'cube.fits' in outfilename:
                brdcfile = outfilename.replace("cube.fits","brd_cont.fits")
            else:
                raise 
        else:
            if 'brd' in brdfile:
                brdcfile = brdfile.replace("brd","brd_cont")
            else:
                raise
        brdc = (np.sum(cubeJHK[2810:2815,:,:],axis=0) + np.sum(cubeJHK[2824:2829,:,:],axis=0))/2.
        fout[0].data=brdc
        fout.writeto(brdcfile,clobber=clobber,output_verify='fix')

    if brefile is not None:
        if brefile is True:
            if 'cube.fits' in outfilename:
                brefile = outfilename.replace("cube.fits","bre.fits")
            else:
                raise 
        bre = np.sum(cubeJHK[2377:2382,:,:],axis=0) # - np.median(cubeJHK[2347:2353,:,:],axis=0)*5
        fout[0].data=bre
        fout.writeto(brefile,clobber=clobber,output_verify='fix')

        if brefile is True:
            if 'cube.fits' in outfilename:
                brecfile = outfilename.replace("cube.fits","bre_cont.fits")
            else:
                raise 
        else:
            if 'bre' in brefile:
                brecfile = brefile.replace("bre","bre_cont")
            else:
                raise
        brec = (np.sum(cubeJHK[2370:2375,:,:],axis=0) + np.sum(cubeJHK[2382:2387,:,:],axis=0))/2.
        fout[0].data=brec
        fout.writeto(brecfile,clobber=clobber,output_verify='fix')

    if pabfile is not None:
        if pabfile is True:
            if 'cube.fits' in outfilename:
                pabfile = outfilename.replace("cube.fits","pab.fits")
            else:
                raise 
        pab = np.sum(cubeJHK[519:523,:,:],axis=0) # - np.median(cubeJHK[2347:2353,:,:],axis=0)*5
        fout[0].data=pab
        fout.writeto(pabfile,clobber=clobber,output_verify='fix')

        if pabfile is True:
            if 'cube.fits' in outfilename:
                pabcfile = outfilename.replace("cube.fits","pab_cont.fits")
            else:
                raise 
        else:
            if 'pab' in pabfile:
                pabcfile = pabfile.replace("pab","pab_cont")
            else:
                raise
        pabc = (cubeJHK[519,:,:]+cubeJHK[523,:,:])*2.
        fout[0].data=pabc
        fout.writeto(pabcfile,clobber=clobber,output_verify='fix')

    for line in line_ranges:
        make_integral(cubeJHK, filename=True, outfilename=outfilename, line=line, fout=fout, clobber=clobber)

    kcont = np.sum([cubeJHK[a:b,:,:].sum(axis=0) for (a,b) in Kcont],axis=0)
    fout[0].data=kcont
    fout.writeto(outfilename.replace('cube','kcont'),clobber=clobber,output_verify='fix')

Kcont = [(3274, 3381), (3597,3742), (3891,4200)]

line_ranges = {
        'h2s3': [2859, 2863],
        'h2s2': [3122, 3127],
        'h2s1': [3427,3433],
        'h2s0': [3779,3785],
        'brg': [3582, 3586],
        'h2q1':[4415,4420],
        'h2q2':[4439,4444],
        'h2q3':[4474,4480],
        'h2q4':[4523,4528],
        'h2q5':[4580,4588],
        'h2s1_2-1':[3864,3870],
        'h2s2_2-1':[3539, 3545],
        'h2s3_2-1':[3258, 3266],
        'h2s4_2-1':[3009, 3011],
        'h2s5_2-1':[2816, 2819],
        'he2.058':[3209, 3214],
        'pab':[517,520],
        'hemaybe':[505,510], # 1.275 microns
        '1.70':[1968,1971],
        'Br10':[2092,2097],
        'Br11':[1899,1904],
        'Br12':[1761, 1766],
        'Br13':[1658, 1662],
        'Br14':[1578, 1582],
        'Fe1.644':[1770,1776]}

def make_integral(cubeJHK, filename, outfilename, line, fout, clobber=True):

    irange = line_ranges[line]

    if filename is True:
        if 'cube.fits' in outfilename:
            filename = outfilename.replace("cube.fits",line+".fits")
        else:
            raise 
    lineim = np.sum(cubeJHK[irange[0]:irange[1],:,:],axis=0)# - np.median(cubeJHK[3435:3445,:,:],axis=0)*4
    fout[0].data=lineim
    fout.writeto(filename,clobber=clobber,output_verify='fix')

def average_cubes(filenames, outfilename, imsize=(4647,300,300), cd1=0.5, cd2=0.5):
    """
    Average a bunch of data cubes together... if you want it to work right, you
    need to make sure the imsize you specify includes all of the cube data.
    """
    #raise
    # this won't work
    image = np.zeros(imsize)
    nhits = np.zeros(imsize[1:])
    outheader = pyfits.Header()
    outheader.fromTxtFile(headerfile)
    outwcs = pywcs.WCS(outheader)
    for file in filenames:
        data = pyfits.getdata(file)
  
        header = pyfits.getheader(file)
        yinds,xinds = np.indices(data.shape[1:])
        raarr = header['CRVAL1'] + (xinds-header['CRPIX1']+1) * header['CD1_1']
        decarr = header['CRVAL2'] + (yinds-header['CRPIX2']+1) * header['CD2_2']

        flathead = agpy.cubes.flatten_header(header)
        flathead['CRPIX1'] = imsize[2]/2.
        flathead['CRPIX2'] = imsize[1]/2.
        flathead['CD1_1'] = -cd1/3600.
        flathead['CD2_2'] =  cd2/3600.
        #wcs = pywcs.WCS(flathead)
  
        xx,yy = outwcs.wcs_sky2pix(raarr,decarr,0)
  
  
        for kk,(ii,jj) in enumerate(zip(xinds.ravel(),yinds.ravel())):
  
            #x,y = wcs.wcs_sky2pix(raarr[jj,ii],decarr[jj,ii],0)
            x,y = xx[kk], yy[kk]
            #print x,y,ii,jj
            #image[int(np.round(x)),int(np.round(y))] += data.TSYS[ii]
            #nhits[int(np.round(x)),int(np.round(y))] += 1
            datavect = data[:,jj,ii]
            OK = (datavect == datavect)
            image[OK,int(np.round(y)),int(np.round(x))]  += datavect[OK]
            nhits[int(np.round(y)),int(np.round(x))]     += 1
  
        imav = image/nhits
  
        for k in ['CRPIX3','CRVAL3','CDELT3','CUNIT3']:
            outheader.update(k,header.get(k))
        HDU = pyfits.PrimaryHDU(data=imav,header=outheader)
        HDU.writeto(outfilename,clobber=True,output_verify='fix')
        # oldwayHDU = pyfits.PrimaryHDU(data=imav,header=header)
        # oldwayHDU.writeto(outfilename,clobber=True,output_verify='fix')

if __name__ == "__main__":
    import glob
    # ref stars = s1 = [83.807522, -5.368632]  dra =0.004632 ddec = 0.008919  dx=34.2 dy=63.8
    #             s2 = [83.812154, -5.377551]  cd = y0.503, x0.487
    make_cube(glob.glob("OMC_BN1_30s_map.*JHK_cal.fits"),"OMC_BN1_30s_map_cal_cube.fits",brgfile=True,
            reference_radec=(83.813977,-5.3711047),reference_pix=(10.9,59.8), rot90=True, cd2=0.503, cd1=0.487)
elif False: #if __name__ == "__main__":
    import glob
    make_cube(glob.glob("BBRW51a_30s_series.*JHK_cal.fits"),"BBRW51a_series1_JHK_cal_cube.fits",brgfile=True,reference_radec=(290.92990,14.51399),reference_pix=(41.5,29.5), rot90=True, cd2=0.97, cd1=0.472)
    make_cube(glob.glob("BBRW51a_30s_series.*JHK_cal_backsub.fits"),"BBRW51a_series1_JHK_cal_backsub_cube.fits",brgfile=True,reference_radec=(290.92990,14.51399),reference_pix=(41.5,29.5), rot90=True, cd2=0.97, cd1=0.472)
    make_cube(glob.glob("BBRW51a_rot_30s_series2.*JHK_cal.fits"),"BBRW51a_rot_series2_JHK_cal_cube.fits",brgfile=True,reference_radec=(290.92679,14.51138),reference_pix=(13.37,44.96), cd1=0.59, cd2=0.457)
    make_cube(glob.glob("BBRW51a_rot_30s_series2.*JHK_cal_backsub.fits"),"BBRW51a_rot_series2_JHK_cal_backsub_cube.fits",brgfile=True,reference_radec=(290.92679,14.51138),reference_pix=(13.37,44.96), cd1=0.59, cd2=0.457)
    make_cube(glob.glob("BBRW51a_rot_30s_series.*JHK_cal.fits"),"BBRW51a_rot_series_JHK_cal_cube.fits",brgfile=True,reference_radec=(290.9308725,14.514113),reference_pix=(33.92,66.10),cd1=0.526,cd2=0.476)
    make_cube(glob.glob("BBRW51a_rot_30s_series.*JHK_cal_backsub.fits"),"BBRW51a_rot_series_JHK_cal_backsub_cube.fits",brgfile=True,reference_radec=(290.9308725,14.514113),reference_pix=(33.92,66.10),cd1=0.526,cd2=0.476)
    make_cube(glob.glob("W51N_30s_series.02[012][0-9]_JHK_cal.fits")+glob.glob("W51N_30s_series.023[0-6]_JHK_cal.fits"),"W51N_series_JHK_cal_cube.fits",brgfile=True,rot90=True,reference_pix=(46.83,7.39),reference_radec=(290.91738,14.51877),cd1=0.5,cd2=0.5)
    make_cube(glob.glob("W51N_30s_series.02[012][0-9]_JHK_cal_backsub.fits")+glob.glob("W51N_30s_series.023[0-6]_JHK_cal_backsub.fits"),"W51N_series_JHK_cal_backsub_cube.fits",brgfile=True,rot90=True,reference_pix=(46.83,7.39),reference_radec=(290.91738,14.51877),cd1=0.5,cd2=0.5)
    average_cubes(glob.glob("*_cal_cube.fits"),"averaged_cube.fits")

    import agpy.agpy_montage
    import os
    os.mkdir('tmpdrive')
    agpy.agpy_montage.wrapper(["*_backsub_brg.fits"], outfile="BrG_backsub_montage_tspec.fits", header='header.hdr', tmpdrive=os.getcwd()+"/tmpdrive/")
    agpy.agpy_montage.wrapper(["*_backsub_brd.fits"], outfile="BrD_backsub_montage_tspec.fits", header='header.hdr', tmpdrive=os.getcwd()+"/tmpdrive/")
    agpy.agpy_montage.wrapper(["*_backsub_bre.fits"], outfile="BrE_backsub_montage_tspec.fits", header='header.hdr', tmpdrive=os.getcwd()+"/tmpdrive/")
    agpy.agpy_montage.wrapper(["*_cal_brg.fits"], outfile="BrG_cal_montage_tspec.fits", header='header.hdr', tmpdrive=os.getcwd()+"/tmpdrive/")
    agpy.agpy_montage.wrapper(["*_cal_brd.fits"], outfile="BrD_cal_montage_tspec.fits", header='header.hdr', tmpdrive=os.getcwd()+"/tmpdrive/")
    agpy.agpy_montage.wrapper(["*_cal_bre.fits"], outfile="BrE_cal_montage_tspec.fits", header='header.hdr', tmpdrive=os.getcwd()+"/tmpdrive/")
    agpy.agpy_montage.wrapper(["*_cal_brg_cont.fits"], outfile="BrG_cont_cal_montage_tspec.fits", header='header.hdr', tmpdrive=os.getcwd()+"/tmpdrive/")
    agpy.agpy_montage.wrapper(["*_cal_brd_cont.fits"], outfile="BrD_cont_cal_montage_tspec.fits", header='header.hdr', tmpdrive=os.getcwd()+"/tmpdrive/")
    agpy.agpy_montage.wrapper(["*_cal_bre_cont.fits"], outfile="BrE_cont_cal_montage_tspec.fits", header='header.hdr', tmpdrive=os.getcwd()+"/tmpdrive/")

    import pyfits
    brg = pyfits.getdata('BrG_cal_montage_tspec.fits')
    bre = pyfits.getdata('BrE_cal_montage_tspec.fits')
    brd = pyfits.getdata('BrD_cal_montage_tspec.fits')

    brgc = pyfits.getdata('BrG_cont_cal_montage_tspec.fits')
    brec = pyfits.getdata('BrE_cont_cal_montage_tspec.fits')
    brdc = pyfits.getdata('BrD_cont_cal_montage_tspec.fits')

    brgcs = brg-brgc
    brecs = bre-brec
    brdcs = brd-brdc


