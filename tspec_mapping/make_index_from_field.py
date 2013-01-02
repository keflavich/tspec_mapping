import astroquery.irsa
import astrometry
import astropy.io.fits
import numpy as np
import atpy
import os

def make_index_from_table(table,fieldname,fov=None,clobber=False,**kwargs):
    """
    Given a table with RA and Dec columns (case-sensitive!), build an astrometry.net
    quad index

    Parameters
    ----------
    table : Table
        astropy.io.table or atpy.table instance (recarrays require different cleaning operations)
    fieldname : str
        Seed name for the output catalog file and index file
    clobber : bool
        Overwrite existing fits files / indices?
    kwargs : 
        Are passed to astrometry.build_index
    """

    #fitstable = astropy.io.fits.BinTableHDU(data=table)
    newtable = atpy.Table(name=fieldname)
    for colname in table.dtype.names:
        newtable.add_column(colname, table[colname])

    # sanitize fieldname
    fieldname = fieldname.replace(" ","_") # what other chars should I be careful of? 

    #fitstable.writeto(fieldname+".fits",clobber=clobber)
    newtable.write(fieldname+".fits",overwrite=clobber)

    if fov is None:
        # guess the FOV... sort of
        rarange = (np.max(table['RA']) - np.min(table['RA']))*3600
        decrange = (np.max(table['Dec'])-np.min(table['Dec']))*3600
        fov = (rarange+decrange)/2.

    return make_index_from_fitstable(fieldname+'.fits',fieldname,fov=fov,**kwargs)

def make_index_from_fitstable(fitstablename, fieldname=None, fov=None, preset_list=None, **kwargs):
    """
    Build an index from a FITS table already on disk (very thin wrapper of build_index)

    Parameters
    ----------
    fitstablename : str
        Full path to a .fits table with the 2nd header being a BinTableHDU for
        astrometry's build-index to parse
    preset_list : list
        List of presets, in the range -5 to 21, to build indices for
    fov : int
        field of view in arcseconds
    fieldname : str
        output prefix for the index file.  If not specified, will use the root string
        of the fitsfilename
    """
    
    if fov is None and 'scale_number' not in kwargs and preset_list is None:
        raise ValueError("Must specify a preset or a FOV")
    elif 'scale_number' in kwargs:
        presets = [kwargs.pop('scale_number')]
    elif preset_list is not None:
        presets = preset_list
    else:
        # determine appropriate "presets" to use
        preset = astrometry.get_closest_preset(fov/60.)
        if preset > -4:
            presets = [preset-2, preset-1,preset,preset+1]
        elif preset > -5:
            presets = [preset-1,preset,preset+1]
        else:
            presets = [preset,preset+1]
    
    if fieldname is None:
        fieldname = os.path.split( os.path.splitext(fitstablename)[0] )[1]

    stdout,stderr = "",""
    for preset in presets:
        _stdout,_stderr = astrometry.build_index(fieldname+".fits",scale_number=preset,**kwargs)
        stdout += _stdout
        stderr += _stderr

    return stdout,stderr

def make_index_from_field_2MASS(coords,fieldname,fov=900,clobber=False,**kwargs):
    """
    Create an index file.  The input should be IRSA-parseable coordinates, e.g.
    a name, ra/dec, or glon/glat coords

    Example
    -------
    >>> make_index_from_field_2MASS('Sgr C','Sgr C',300,scan_catalog=True,clobber=True)
    >>> make_index_from_field_2MASS('266.1512 -29.4703','Sgr C',300,scan_catalog=True,clobber=True)
    >>> make_index_from_field_2MASS('359.4288 -00.0898 gal','Sgr C',300,scan_catalog=True,clobber=True)
    
    """

    fieldname = fieldname.replace(" ","_") # what other chars should I be careful of? 

    table = astroquery.irsa.query_gator_box('pt_src_cat',coords,fov)
    table.rename_column('ra','RA')
    table.rename_column('dec','Dec')
    cleantable = _clean_table(table)

    return make_index_from_table(cleantable,fieldname,fov=fov,clobber=clobber,**kwargs)

def make_index_from_field_UKIDSS(glon,glat,fieldname,catalog='GPS',fov=900,clobber=False,**kwargs):
    """
    Create an index file.  The input should be UKIDSS-parseable coordinates, e.g.
    glon,glat (so far, only a galactic lon/lat query tool is implemented

    Example
    -------
    >>> make_index_from_field_UKIDSS(359.4288,-00.0898,'Sgr C',fov=300,scan_catalog=True,clobber=True)
    
    """

    fieldname = fieldname.replace(" ","_") # what other chars should I be careful of? 

    ukquery = astroquery.ukidss.UKIDSSQuery()
    ukquery.programmeID = catalog
    uktable = ukquery.get_catalog_gal(glon,glat,radius=fov/60.)[0]
    uktable.writeto(fieldname+".fits",clobber=clobber)
    #bintab = table[0][1]
    #bintab.data = bintab.data.astype(newtype)
    #table.rename_column('ra','RA')
    #table.rename_column('dec','Dec')
    #cleantable = _clean_table(table)

    return make_index_from_fitstable(fieldname+".fits",fieldname=fieldname,fov=fov,**kwargs)


def _clean_table(table):
    """
    Hack to convert a table to a FITS-friendly numpy ndarray;
    this will become obsolete when astropy's table includes a FITS writer
    """
    float_types = [np.float, np.float128, np.float16, np.float32, np.float64, np.float_, np.floating]
    new_fields = [(k,np.dtype('S8')) if v[0].type == np.object_ else 
                  (k,np.float64) if v[0].type in float_types else (k,v[0])  
                  for (k,v) in table._data.dtype.fields.iteritems()]
    new_array = np.array(table._data, dtype=new_fields)

    return new_array


