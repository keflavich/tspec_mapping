import astroquery.irsa

def make_index_from_field(coords,fieldname,fov=900):
    """
    Input coords... arcsecs...
    """

    table = astroquery.irsa.query_gator_box('pt_src_cat',coords,fov)
    table.write(fieldname+".fits")

    #T = atpy.Table()
    #[T.add_column(n,ORION[n]) for n in ORION.columns]
    #T.write('2MASS_orion_cat.fits')
    #T2 = atpy.Table()
    #T2.add_column('RA',T.ra)
    #T2.add_column('Dec',T.dec)
    #T2.add_column('Kmag',T.k_m)
    #T2.write('2MASS_orion_cat_simple.fits')
