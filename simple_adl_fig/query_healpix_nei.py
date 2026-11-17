#!/usr/bin/env python
"""
Generic python script.

Modified for DELVE-MC by 
Mahdieh Navabi
"""
__author__ = "Sidney Mau"
import pandas as pd
import numpy as np
import healpy as hp
from astropy.table import Table
from dl import queryClient as qc

#-------------------------------------------------------------------------------

# DES Redenning coefficients

#R_u = 3.9631
R_g = 3.1863
R_r = 2.1401
R_i = 1.5690
R_z = 1.196
R_y = 1.048 

# Define the calibration function
def apply_calib(list, healpix_index):
    number = np.where(list['# id'] == healpix_index)[0][0]
    # Define calibration values
    calib_r = list['r'][number]
    calib_g = list['g'][number]
    calib_i = list['i'][number]
    return calib_g, calib_r, calib_i

def nside32toring256(pixel_32_array):
    """
    Find all HEALPix pixels at nside=256 related to a given nside=32 pixel and
    construct and execute a SQL query to fetch data for these pixels.

    Parameters:
    pixel_32 (int): The HEALPix pixel number at nside=32.

    Returns:
    result: all HEALPix pixels at nside=256
    """
    all=[]
    for pixel_32 in pixel_32_array:
        # Define the nside values
        nside_32 = 32
        nside_256 = 256

        # Get the boundary vertices of the pixel at nside=32
        vertices_theta_phi = hp.boundaries(nside_32, pixel_32, nest=False)

        # Convert vertices to a list of coordinates
        vertices = np.transpose(vertices_theta_phi)

        # Find all nside=256 pixels within the boundary using query_polygon
        pixels_256 = hp.query_polygon(nside_256, vertices, nest=False)

        # Convert the list of pixels to a comma-separated string
        pixels_256_str = ','.join(map(str, pixels_256))

        # Print the resulting string
        # print("Formatted pixel list for SQL IN clause:", pixels_256_str)

        # Define the calibration values
        calib_g, calib_r, calib_i = apply_calib(list_calib_values, pixel_32)

        all.append((pixels_256_str, calib_g, calib_r, calib_i))
    return all

# Define the main path
list_calib_values = pd.read_csv('/user/HS501/mn00962/Documents/delve_calibration/table_calibration.csv')

#--------------------------------------------------------------------------------

def query(profile, healpix, gmax, stars=True, galaxies=False):
    """Return data queried from datalab

    Parameters
    ----------
    profile : Profile for data lab query [str]
    ra      : Right Ascension [deg]
    dec     : Declination [deg]
    radius  : Radius around (ra, dec) [deg]

    Returns
    -------
    data : numpy recarray of data
    """
    
    # Define the healpix pixels in ring256    
    all_results = nside32toring256(healpix)
    # Unpack results from nside32toring256
    all_pixels_256 = ','.join([result[0] for result in all_results])
    calib_g = all_results[0][1]
    calib_r = all_results[0][2]
    calib_i = all_results[0][3]
    # import pdb; pdb.set_trace()
    qc.set_profile(profile)

    sql_stars = f'''
        SELECT ra,
               dec,
               gmag,
               (gmag+{calib_g}-{R_g}*ebv) AS gmag0,
               gerr,
               rmag,
               (rmag+{calib_r}-{R_r}*ebv) AS rmag0,
               rerr,
               imag,
               (imag+{calib_i}-{R_i}*ebv) AS imag0,
               ierr,
               ebv,
               ndetg,
               ndetr,
               chi,
               prob,
               sharp,
               brickuniq,
               nalfdetiter,
               ring256,
               brick
        FROM delvemc_y4t2.object 
        WHERE ring256 IN ({all_pixels_256})
              AND brickuniq = 1
              AND gmag < 23
              AND ndetg > 2
              AND ndetr > 2
              AND abs(gmag-{R_g}*ebv - rmag-{R_r}*ebv) < 1.2 
              AND prob > 0.35
              AND abs(sharp) < 0.8
    '''
            #   AND chi < 2.5
 
    sql_galaxies = f'''
        SELECT ra,
               dec,
               gmag,
               (gmag+{calib_g}-{R_g}*ebv) AS gmag0,
               gerr,
               rmag,
               (rmag+{calib_r}-{R_r}*ebv) AS rmag0,
               rerr,
               imag,
               (imag+{calib_i}-{R_i}*ebv) AS imag0,
               ierr,
               ebv,
               ndetg,
               ndetr,
               chi,
               prob,
               sharp,
               brickuniq,
               ring256,
               brick
        FROM delvemc_y4t2.object
        WHERE ring256 IN ({all_pixels_256})
              AND brickuniq = 1
              AND gmag < 23
              AND ndetg > 2
              AND ndetr > 2
              AND abs(gmag-{R_g}*ebv - rmag-{R_r}*ebv) < 1.2
              AND prob < 0.35
              AND abs(sharp)> 0.8
    '''
                #  AND chi>5

    if stars:
        query_stars = qc.query(sql=sql_stars, fmt='structarray')
    if galaxies:
        query_galaxies = qc.query(sql=sql_galaxies, fmt='structarray')

    if stars and not galaxies:
        return query_stars
    elif not stars and galaxies:
        return query_galaxies
    elif stars and galaxies:
        return query_stars, query_galaxies
    else:
        return None

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', type=str, required=False,
                        help='Profile for data lab query [str]')
    parser.add_argument('--ra', type=float, required=True,
                        help='Right Ascension of target position [deg]')
    parser.add_argument('--dec', type=float, required=True,
                        help='Declination of target position [deg]')
    parser.add_argument('--radius', type=float, default=1.0,
                        help='Radius around target position [deg]')
    parser.add_argument('--gmax', type=float, default=24,
                        help='Maximum g-band magnitude [mag]')
    args = parser.parse_args()
    data = query(args.profile, args.ra, args.dec, args.radius, args.gmax)
