import astroquery.irsa
import astrometry
import astropy.io.fits
import numpy as np
import atpy

def make_index_from_field(coords,fieldname,fov=900,clobber=False,**kwargs):
    """
    Input coords... arcsecs...

    Example
    -------
    make_index_from_field('Sgr C','Sgr C',300,scan_catalog=True,scale_number=-1,clobber=True)
    tbl = astroquery.irsa.query_gator_box('pt_src_cat','83.808 -5.391',300)
    
    """

    table = astroquery.irsa.query_gator_box('pt_src_cat',coords,fov)
    table.rename_column('ra','RA')
    table.rename_column('dec','Dec')
    cleantable = _clean_table(table)
    fitstable = astropy.io.fits.BinTableHDU(data=cleantable)
    newtable = atpy.Table()
    for colname in table.columns:
        newtable.add_column(colname, table[colname])

    # sanitize fieldname
    fieldname = fieldname.replace(" ","_") # what other chars should I be careful of? 

    #fitstable.writeto(fieldname+".fits",clobber=clobber)
    newtable.write(fieldname+".fits",overwrite=clobber)

    status = astrometry.build_index(fieldname+".fits",**kwargs)

    return status

    #T = atpy.Table()
    #[T.add_column(n,ORION[n]) for n in ORION.columns]
    #T.write('2MASS_orion_cat.fits')
    #T2 = atpy.Table()
    #T2.add_column('RA',T.ra)
    #T2.add_column('Dec',T.dec)
    #T2.add_column('Kmag',T.k_m)
    #T2.write('2MASS_orion_cat_simple.fits')

def _clean_table(table):
    """
    Hack to convert a table to a FITS-friendly numpy ndarray;
    this will become obsolete when astropy's table includes a FITS writer
    """
    new_fields = [(k,np.dtype('S%i' % v[1])) if v[0].type == np.object_ else (k,v)  
                  for (k,v) in table._data.dtype.fields.iteritems()]
    new_array = np.array(table._data, dtype=new_fields)

    return new_array


