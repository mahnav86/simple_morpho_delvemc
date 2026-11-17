#!/usr/bin/env python
"""
Generic python script.

Modified for DELVE-MC by 
Joanna Sakowska
"""
__author__ = "Sidney Mau"

import time

import numpy as np
from dl import queryClient as qc


def _query_retry(sql, attempts=3):
    """Data Lab sometimes stalls; retry before giving up."""
    for attempt in range(1, attempts + 1):
        try:
            return qc.query(sql=sql, fmt='structarray')
        except Exception as exc:
            if attempt == attempts:
                raise
            print(f"Data Lab query failed ({type(exc).__name__}), retrying ({attempt}/{attempts})...")
            time.sleep(10 * attempt)

#-------------------------------------------------------------------------------

# DES Redenning coefficients

#R_u = 3.9631
R_g = 3.1863
R_r = 2.1401
R_i = 1.5690
R_z = 1.196
R_y = 1.048 


#--------------------------------------------------------------------------------

def query(profile, ra, dec, radius, gmax, stars=True, galaxies=False):
    """Return data queried from datalab
    Parameters
    ----------
    profile : Profile for data lab query [str]
    ra      : Right Ascension [deg]
    dec     : Declination [deg]
    radius  : radius around (ra, dec) [deg]

    Returns
    -------
    data : numpy recarray of data
    """


    qc.set_profile(profile)


    sql_stars = f'''
        SELECT ra,
               dec,
               gmag,
               gmag-{R_g}*ebv AS gmag0,
               gerr,
               rmag,
               rmag-{R_r}*ebv AS rmag0,
               rerr,
               imag,
               imag-{R_i}*ebv AS imag0,
               ierr,
               ebv,
               ndet,
               ndetr,
               chi,
               prob,
               sharp,
               brickuniq,
               nalfdetiter
        FROM delve_dr3.mc_object 
        WHERE q3c_radial_query(ra,dec,{ra},{dec},{radius})
              AND brickuniq = TRUE
              AND rmag < 23.5
              AND ndet > 5
              AND abs(gmag-{R_g}*ebv - rmag-{R_r}*ebv) < 1.2 
              AND prob >0.6
              AND abs(sharp) < 0.8
              AND chi < 3
    '''
            #  delvemc_y4t2.object 


    sql_galaxies = f'''
        SELECT ra,
               dec,
               gmag,
               gmag-{R_g}*ebv AS gmag0,
               gerr,
               rmag,
               rmag-{R_r}*ebv AS rmag0,
               rerr,
               imag,
               imag-{R_i}*ebv AS imag0,
               ierr,
               ebv,
               ndet,
               chi,
               prob,
               sharp,
               brickuniq
        FROM delve_dr3.mc_object
        WHERE q3c_radial_query(ra,dec,{ra},{dec},{radius})
              AND brickuniq = TRUE
              AND rmag < 23.5
              AND ndet > 5
              AND abs(gmag-{R_g}*ebv - rmag-{R_r}*ebv) < 1.2
              AND prob < 0.4
              AND abs(sharp)> 0.9             
    '''
              #AND chi > 3 
                #   AND brickuniq = 1


    if stars:
        query_stars = _query_retry(sql_stars)
    if galaxies:
        query_galaxies = _query_retry(sql_galaxies)


    if stars and not galaxies:
        return(query_stars)
    elif not stars and galaxies:
        return(query_galaxies)
    elif stars and galaxies:
        return(query_stars, query_galaxies)
    else:
        return(None)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile',type=str,required=False,
                        help='Profile for data lab query [str]')
    parser.add_argument('--ra',type=float,required=True,
                        help='Right Ascension of target position [deg]')
    parser.add_argument('--dec',type=float,required=True,
                        help='Declination of target position [deg]')
    parser.add_argument('--radius',type=float,default=1.0,
                        help='Radius around target position [deg]')
    parser.add_argument('--gmax',type=float,default=24,
                        help='Maximum g-band magnitude [mag]')
    args = parser.parse_args()
    data = query(args.profile, args.ra, args.dec, args.radius, args.gmax)
