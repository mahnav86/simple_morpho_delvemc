
#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Sid Mau"

# Python libraries
import os
from sre_constants import AT_END_STRING
import yaml
import numpy as np
import healpy as hp
import scipy.ndimage
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import simple_adl_fig.survey
import simple_adl_fig.isochrone
import simple_adl_fig.coordinate_tools
from simple_adl_fig.search import cut_isochrone_path

#-------------------------------------------------------------------------------

props = dict(facecolor='white', edgecolor='black', linewidth=1) #dict(facecolor='white', edgecolor='none', alpha=0.7)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=False,default='config.yaml',
                        help='Config file [.yaml]')
    parser.add_argument('--outfile',type=str,required=False,default='out.png',
                        help='Output file [.png/.pdf]')
    parser.add_argument('--ra',type=float,required=True,
                        help='Right Ascension of target position [deg]')
    parser.add_argument('--dec',type=float,required=True,
                        help='Declination of target position [deg]')
    parser.add_argument('--mod',type=float,required=False,
                        help='Distance modulus [mag]')
    args = parser.parse_args()

    with open(args.config, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
        survey = simple_adl_fig.survey.Survey(cfg)

    #---------------------------------------------------------------------------

    region = simple_adl_fig.survey.Region(survey, args.ra, args.dec)
    print('Plot coordinates: (RA, Dec) = ({:0.2f}, {:0.2f})'.format(region.ra, region.dec))

    # Apply filters

    data = region.load_data(stars=True, galaxies=False) # Data is stars as default
    #import pdb; pdb.set_trace()

    stars = data

    galaxies = region.load_data(stars=False, galaxies=True)
    

    print('Found {} objects'.format(len(stars)))
    #if (len(stars) == 0):
    #    print('Ending search.')
    #    exit()



    #---------------------------------------------------------------------------

    iso = simple_adl_fig.isochrone.Isochrone(survey=survey.isochrone['survey'],
                              band_1=survey.band_1.lower(),
                              band_2=survey.band_2.lower(),
                              age=12.0, #survey.isochrone['age'],
                              metallicity=0.00010) #survey.isochrone['metallicity']
    if args.mod:
        iso.distance_modulus = args.mod
    #iso_sep = iso.separation(data[mag_1], data[mag_2])

    
    # Isochrone for stars


    iso_filter_stars = cut_isochrone_path(
        stars[survey.mag_dered_1], 
        stars[survey.mag_dered_2], 
        stars[survey.mag_err_1], 
        stars[survey.mag_err_2], 
        iso, 
        mag_max=survey.catalog['mag_max'], 
        radius=0.1, 
        return_all=False)

        
    
    # projection of image
    proj = region.proj
    x_stars, y_stars = proj.sphereToImage(stars[survey.catalog['basis_1']], stars[survey.catalog['basis_2']])
    x_galaxies, y_galaxies = proj.sphereToImage(galaxies[survey.catalog['basis_1']], galaxies[survey.catalog['basis_2']])

    # hess
    mag = stars[survey.mag_dered_1]
    color = stars[survey.mag_dered_1] - stars[survey.mag_dered_2]
    
    # Filters
    
    # OLD ESTIMATE

    #extension = 0.05 # SM setting

    #if args.radius is not None:
    #    r_estimate = args.radius
    #else:
    #    r_estimate = extension

    #r0 = 3.0 * r_estimate # 3.0 # set as g-radius 
    #r1 = 5.0 * r_estimate # 5.0 
    #r2 = np.sqrt(r0**2 + r1**2) 

    # NEW ESTIMATE
    #  
    extension = 0.025 #JS: Minimised the extension
    r0 = 1.5 * extension # 3.0 # set as g-radius 
    r1 = 3.0 * extension # 5.0 rnear? 
    r2 = np.sqrt(r0**2 + r1**2) # rfar?
    angsep_stars = simple_adl_fig.coordinate_tools.angsep(args.ra, args.dec, stars[survey.catalog['basis_1']], stars[survey.catalog['basis_2']])
    inner = (angsep_stars < r0)
    outer = ((angsep_stars > r1) & (angsep_stars < r2))
    background = (angsep_stars > r2)

    #---------------------------------------------------------------------------

    fig, axs = plt.subplots(nrows=2,  ncols=3, figsize=(16, 8), dpi=300)
    fig.subplots_adjust(wspace=0.5)

#     #---------------------------------------------------------------------------
    
#     # Stellar histogram
#     ax = axs[0,0]
#     plt.sca(ax)

#     #r0 = g_radius
    
#     bound = 2
#     steps = 100
#     bins = np.linspace(-bound, bound, steps)
#     signal = np.histogram2d(x_stars[iso_filter_stars], y_stars[iso_filter_stars], bins=[bins, bins])[0]
#     sigma = 0.01 * (0.25 * np.arctan(0.25 * r0 * 60. - 1.5) + 1.3)
#     convolution = scipy.ndimage.filters.gaussian_filter(signal, sigma/(bound/steps)).T
#     pc = ax.pcolormesh(bins, bins, convolution, cmap='Greys', rasterized=True)
    
#     # search kernel
#     #x, y = ra_proj, dec_proj
#     #delta_x = 0.01
#     #area = delta_x**2
#     #smoothing = 2. / 60. # Was 3 arcmin
#     #bins = np.arange(-8., 8. + 1.e-10, delta_x)
#     #centers = 0.5 * (bins[0: -1] + bins[1:])
#     #yy, xx = np.meshgrid(centers, centers)
#     #h = np.histogram2d(x, y, bins=[bins, bins])[0]
#     #h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)
#     #pc = ax.pcolormesh(bins, bins, h_g.T, cmap='Greys', rasterized=True)
    
#     ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top', bbox=props)
    
#     ax.set_xlim(0.5, -0.5)
#     ax.set_ylim(-0.5, 0.5)
#     ax.set_xlabel(r'$\Delta {\rm RA}$ (deg)')
#     ax.set_ylabel(r'$\Delta {\rm Dec}$ (deg)')
    
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes('right', size='5%', pad=0)
#     cb = fig.colorbar(pc, cax=cax)
#     #cb.set_label('Counts')

#     #---------------------------------------------------------------------------
    
#     # Galactic histogram
#     ax = axs[0,1]
#     plt.sca(ax)

#     bound = 0.5
#     steps = 100
#     bins = np.linspace(-bound, bound, steps)
#     signal = np.histogram2d(x_galaxies, y_galaxies, bins=[bins, bins])[0]
#     sigma = 0.01 * (0.25 * np.arctan(0.25 * r0 * 60. - 1.5) + 1.3)
#     convolution = scipy.ndimage.filters.gaussian_filter(signal, sigma/(bound/steps)).T
#     pc = ax.pcolormesh(bins, bins, convolution, cmap='Greys', rasterized=True, vmin=vmin, vmax=vmax)
    
#     ## search kernel
#     ##x, y = ra_proj, dec_proj
#     ##delta_x = 0.01
#     ##area = delta_x**2
#     ##smoothing = 2. / 60. # Was 3 arcmin
#     ##bins = np.arange(-8., 8. + 1.e-10, delta_x)
#     ##centers = 0.5 * (bins[0: -1] + bins[1:])
#     ##yy, xx = np.meshgrid(centers, centers)
#     ##h = np.histogram2d(x, y, bins=[bins, bins])[0]
#     ##h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)
#     ##pc = ax.pcolormesh(bins, bins, h_g.T, cmap='Greys', rasterized=True)
    
#     ax.text(0.05, 0.95, 'Galaxies', transform=ax.transAxes, verticalalignment='top', bbox=props)
    
#     ax.set_xlim(0.25, -0.25)
#     ax.set_ylim(-0.25, 0.25)
#     ax.set_xlabel(r'$\Delta \alpha_{2000}$ (deg)')
#     ax.set_ylabel(r'$\Delta \delta_{2000}$ (deg)')
#     # Get the min and max for the color scale


#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes('right', size='5%', pad=0)
#     cb = fig.colorbar(pc, cax=cax)
#     #cb.set_label('Counts')
# #---------------------------------------------------------------------------
   
    # Calculate axis limits and color scale based on both histograms
    
    common_xlim = (0.5, -0.5)
    common_ylim = (-0.5, 0.5)

    # Calculate the global color scale based on the two datasets
    bound = 0.5
    steps = 100
    bins = np.linspace(-bound, bound, steps)
    signal_stars = np.histogram2d(x_stars[iso_filter_stars], y_stars[iso_filter_stars], bins=[bins, bins])[0]
    signal_galaxies = np.histogram2d(x_galaxies, y_galaxies, bins=[bins, bins])[0]

    sigma_stars = 0.01 * (0.25 * np.arctan(0.25 * r0 * 60. - 1.5) + 1.3)
    sigma_galaxies = 0.01 * (0.25 * np.arctan(0.25 * r0 * 60. - 1.5) + 1.3)

    convolution_stars = scipy.ndimage.filters.gaussian_filter(signal_stars, sigma_stars/(bound/steps)).T
    convolution_galaxies = scipy.ndimage.filters.gaussian_filter(signal_galaxies, sigma_galaxies/(bound/steps)).T

    # Get the min and max for the color scale
    vmin = convolution_stars.min()
    vmax = convolution_stars.max()

    # Stellar histogram
    ax = axs[0, 0]
    plt.sca(ax)

    pc = ax.pcolormesh(bins, bins, convolution_stars, cmap='Greys', rasterized=True, vmin=vmin, vmax=vmax)
    ax.set_xlim(common_xlim)
    ax.set_ylim(common_ylim)
    ax.set_xlabel(r'$\Delta {\rm RA}$ (deg)')
    ax.set_ylabel(r'$\Delta {\rm Dec}$ (deg)')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0)
    cb = fig.colorbar(pc, cax=cax)

    # Galactic histogram
    ax = axs[0, 1]
    plt.sca(ax)

    vmin = convolution_galaxies.min()
    vmax = convolution_galaxies.max()
    
    pc = ax.pcolormesh(bins, bins, convolution_galaxies, cmap='Greys', rasterized=True, vmin=vmin, vmax=vmax)
    ax.set_xlim(common_xlim)
    ax.set_ylim(common_ylim)
    ax.set_xlabel(r'$\Delta \alpha_{2000}$ (deg)')
    ax.set_ylabel(r'$\Delta \delta_{2000}$ (deg)')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0)
    cb = fig.colorbar(pc, cax=cax)

#---------------------------------------------------------------------------

    # selecting a region of interest
    annulus_stars = (angsep_stars > r0) & (angsep_stars < 0.5)
    nbhd_stars = (angsep_stars < r0)
#---------------------------------------------------------------------------

    # Star plot  

    ax = axs[0,2]    
    plt.sca(ax)
    
    # Take region
    proj = region.proj

    # Checking to see how not isochrone filtered stars are distributed:
    # This is *not* the original setting
    x_stars_no_iso, y_stars_no_iso = proj.sphereToImage(stars[survey.catalog['basis_1']], stars[survey.catalog['basis_2']])
    ax.scatter(x_stars_no_iso, y_stars_no_iso, edgecolor='none', s=1, alpha=0.5, c='black')

    # Below plots stars that are isochrone filtered *already*
    # This is the original plot setting
    x_stars, y_stars = proj.sphereToImage(stars[iso_filter_stars][survey.catalog['basis_1']], stars[iso_filter_stars][survey.catalog['basis_2']])

    #ax.scatter(x_stars, y_stars, edgecolor='none', s=3, c='black')
    # ax.scatter(x_stars, y_stars, edgecolor='none', s=7, c='blue', label='$\Delta$CM < 0.1' )

    # Overplot stars in nbhd
    x_stars_nbhd, y_stars_nbhd = proj.sphereToImage(stars[nbhd_stars][survey.catalog['basis_1']], stars[nbhd_stars][survey.catalog['basis_2']])
    ax.scatter(x_stars_nbhd, y_stars_nbhd, edgecolor='none', s=2, c='g', label='r < {:.3f}$^\circ$'.format(r0))

    # Overplot isochrone filtered stars in nbhd
    x_stars_nbhd_iso, y_stars_nbhd_iso = proj.sphereToImage(stars[iso_filter_stars & nbhd_stars][survey.catalog['basis_1']], stars[iso_filter_stars & nbhd_stars][survey.catalog['basis_2']])
    ax.scatter(x_stars_nbhd_iso, y_stars_nbhd_iso, edgecolor='none', s=5, c='r', label='$\Delta$CM < 0.1 & r < {:.3f}$^\circ$'.format(r0))

    ax.set_xlim(0.25, -0.25)
    ax.set_ylim(-0.25, 0.25)
    ax.set_xlabel(r'$\Delta$ RA (deg)')
    ax.set_ylabel(r'$\Delta$ Dec (deg)')
    ax.legend(loc='upper right')

    #ax.set_title('Stars')
   
    #---------------------------------------------------------------------------

    # CMD stars

    ax = axs[1,0]
    plt.sca(ax)
 

    ## Plot background objects
    # ax.scatter(stars[survey.mag_dered_1][annulus_stars & iso_filter_stars] - stars[survey.mag_dered_2][annulus_stars & iso_filter_stars], stars[survey.mag_dered_1][annulus_stars & iso_filter_stars], c='k', alpha=0.1, edgecolor='none', s=1)

    #  Plot background objects
    ax.scatter(stars[annulus_stars][survey.mag_dered_1] - stars[annulus_stars][survey.mag_dered_2], stars[annulus_stars][survey.mag_dered_1], c='k', alpha=0.2, edgecolor='none', s=1)

    ## Plot isochrone
    ax.plot(iso.color, iso.mag + iso.distance_modulus, c='k', lw=1)

    ## Plot objects in nbhd
    ax.scatter(stars[nbhd_stars][survey.mag_dered_1] - stars[nbhd_stars][survey.mag_dered_2], stars[nbhd_stars][survey.mag_dered_1], c='g', s=5, label='r < {:.3f}$^\circ$'.format(r0))

    ## Plot objects in nbhd and near isochrone
    ax.scatter(stars[iso_filter_stars & nbhd_stars][survey.mag_dered_1] - stars[survey.mag_dered_2][iso_filter_stars & nbhd_stars], stars[nbhd_stars & iso_filter_stars][survey.mag_dered_1], c='r', s=5, label='$\Delta$CM < 0.1 & r < {:.3f}$^\circ$'.format(r0))

    ax.legend(loc='upper right')

    ax.set_xlim(-0.5, 1)
    ax.set_ylim(survey.catalog['mag_max'], 16)
    ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top', bbox=props)

    ax.set_xlabel('${} - {}$'.format(survey.band_1.lower(), survey.band_2.lower()))
    ax.set_ylabel('${}$'.format(survey.band_1.lower()))
    
    #---------------------------------------------------------------------------

    # CMD galaxies    
    
    ax = axs[1,1]    
    plt.sca(ax)


    angsep_galaxies = simple_adl_fig.coordinate_tools.angsep(args.ra, args.dec, galaxies[survey.catalog['basis_1']], galaxies[survey.catalog['basis_2']])

    annulus_galaxies = (angsep_galaxies > r0) & (angsep_galaxies < 0.5)
    nbhd_galaxies = (angsep_galaxies < r0)

    # isochrone for galaxies 

    iso_filter_galaxies = cut_isochrone_path(
        galaxies[survey.mag_dered_1], 
        galaxies[survey.mag_dered_2], 
        galaxies[survey.mag_err_1], 
        galaxies[survey.mag_err_2], 
        iso, 
        mag_max=survey.catalog['mag_max'], 
        radius=0.1, 
        return_all=False)

    # Plot background objects
    ax.scatter(galaxies[survey.mag_dered_1][annulus_galaxies & iso_filter_galaxies] - galaxies[survey.mag_dered_2][annulus_galaxies & iso_filter_galaxies], galaxies[survey.mag_dered_1][annulus_galaxies & iso_filter_galaxies], c='k', alpha=0.1, edgecolor='none', s=1)

    # Galaxies
    ax.scatter(galaxies[survey.mag_dered_1][annulus_galaxies] - galaxies[survey.mag_dered_2][annulus_galaxies], galaxies[survey.mag_dered_1][annulus_galaxies], c='k', alpha=0.1, edgecolor='none', s=1)

    # Plot isochrone
    ax.plot(iso.color, iso.mag + iso.distance_modulus, c='k', lw=1)

    # Plot objects in nbhd
    ax.scatter(galaxies[survey.mag_dered_1][nbhd_galaxies] - galaxies[survey.mag_dered_2][nbhd_galaxies], galaxies[survey.mag_dered_1][nbhd_galaxies], c='g', s=5, label='r < {:.3f}$^\circ$'.format(r0))

    # Plot objects in nbhd and near isochrone
    ax.scatter(galaxies[survey.mag_dered_1][nbhd_galaxies & iso_filter_galaxies] - galaxies[survey.mag_dered_2][nbhd_galaxies & iso_filter_galaxies], galaxies[survey.mag_dered_1][nbhd_galaxies & iso_filter_galaxies], c='r', s=5, label='$\Delta$CM < 0.1  & r < {:.3f}$^\circ$'.format(r0))

    ax.legend(loc='upper right')

    ax.set_xlim(-0.5, 1)
    ax.set_ylim(survey.catalog['mag_max'], 16)
    ax.set_xlabel('${} - {}$'.format(survey.band_1.lower(), survey.band_2.lower()))
    ax.set_ylabel('${}$'.format(survey.band_1.lower()))
    ax.text(0.05, 0.95, 'Galaxies', transform=ax.transAxes, verticalalignment='top', bbox=props)

 #---------------------------------------------------------------------------
    
    # Hess
    #ax = axs[1]
    ax = axs[1,2]
    plt.sca(ax)
    
    xbins = np.arange(-0.3, 1.1, 0.1)
    ybins = np.arange(16., 24.0 + 0.5, 0.5)
    foreground = np.histogram2d(color[inner], mag[inner], bins=[xbins, ybins])
    background = np.histogram2d(color[outer], mag[outer], bins=[xbins, ybins])
    fg = foreground[0].T
    bg = background[0].T
    fg_abs = np.absolute(fg)
    bg_abs = np.absolute(bg)
    mask_abs = fg_abs + bg_abs
    mask_abs[mask_abs == 0.] = np.nan # mask common zeroes
    signal = fg - bg
    signal = np.ma.array(signal, mask=np.isnan(mask_abs)) # mask nan
    pc = ax.pcolormesh(xbins, ybins, signal, cmap='viridis', rasterized=True)
    if args.mod:
        plt.plot(iso.color, iso.mag+iso.distance_modulus, lw=1, c='k')
    
    ax.set_xlim(-0.3, 1.0)
    ax.set_ylim(24.0, 16.0)
    ax.set_xlabel('${} - {}$'.format(survey.band_1.lower(), survey.band_2.lower()))
    ax.set_ylabel('${}$'.format(survey.band_1.lower()))
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0)
    cb = fig.colorbar(pc, cax=cax)
    #cb.set_label('Counts')

    #---------------------------------------------------------------------------
   
    file_name = 'candidate_{:0.2f}_{:0.2f}_{:2.2f},r={:0.2f}'.format(args.ra, args.dec, args.mod, r0)
    fig.savefig('plot_dir/{}.png'.format(file_name), bbox_inches='tight')
    # fig.savefig(os.path.join(survey.output['save_dir'],args.outfile), bbox_inches='tight')
    plt.close(fig)
    
    selected_stars = stars[nbhd_stars]

    # Convert to DataFrame if needed
    df = pd.DataFrame(selected_stars)
    # df.to_csv('files/nbhd_stars.csv', index=False)
    df.to_csv('files/{}.csv'.format(file_name), index=False)
