import sys
from saraswati_radio_plots import smearing,write_catalog,write_catalog_20sigma,spectral_index,plot_image,spatial_spectral_index,rms_plot,spectral_index_distance,spectral_distributions,radio_lumo_distance,radio_lumo_distance,radio_power_distributions,astrometry_NVSS,astrometry_VLASS,flux_scale_NVSS,flux_scale_VLASS,flux_distribution,resolved_unresolved_sources,redshift_distributions,simulation,astrometry_FIRST,flux_scale_FIRST,cumulative_rms_map,radio_cutouts,simulation_catalog,completness,source_counts,angular_distribution
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import numpy as np
from astropy.io import fits

def zwcl_main( argv ):
    
    cluster_name='Zwcl2341'
    
    output_path= cluster_name+'/plots/'

    simulation_path=cluster_name+'/simulated/'

    fits_files=cluster_name+'/input/'
    # ---------------------------------------------------------
    #  SDSS and DeCALS combined catalogues for photo vs spec tests
    # ---------------------------------------------------------

    radio_catalogue_fits = fits_files+cluster_name+'_srl.fits'
   
    x=7
    radio_optical_catalogue_fits = fits.open(fits_files+cluster_name+'_DeCALS_radio_join.fits')
    radio_optical_catalog=radio_optical_catalogue_fits[1].data


    tessel_region=fits_files+'zwcl.tessel.reg'
    app_image=fits_files+'mypipelinerun_ZwCl2341_1_p_0000_4-MFS-image.fits'#'mypipelinerun_ZwCl2341_1_p_0000_4-MFS-image.fits'
    int_image=fits_files+'mypipelinerun_ZwCl2341_1_p_0000_4-MFS-image.fits'#'mypipelinerun_ZwCl2341_1_p_0000_4-MFS-image.fits'
    fits_image=fits_files+'mypipelinerun_ZwCl2341_1_p_0000_4-MFS-image.fits'
    
    radio_HSC_catalog='Zwcl2341_HSC_radio_join.fits'
    
    rms_image=fits_files+'Zwcl2341_rms_map.fits'
    res_image=fits_files+'Zwcl2341_res_map.fits'
    # ---------------------------------------------------------------
    #  cluster general paramaters 
    # ---------------------------------------------------------------


    cluster_centre=SkyCoord(str(355.91541666667), str(0.33083333333333), frame='icrs',unit=(u.deg,u.deg))

    search_radius = 1.5

    cluster_z_lit=0.2694

    plotting_scripts(output_path,fits_files,simulation_path,cluster_name,radio_catalogue_fits,radio_optical_catalog,cluster_centre,app_image,int_image,res_image,fits_image,tessel_region,rms_image,radio_HSC_catalog)


     
def A2631_main( argv ):

    cluster_name='A2631'
    output_path= cluster_name+'/plots/'

    simulation_path=cluster_name+'/simulated/'

    fits_files=cluster_name+'/input/'

    radio_catalogue_fits = fits_files+cluster_name+'_srl.fits'
 
    COSMOS_catalogue_fits= fits_files+'vla3_cosmos_sources_160321_public5sig.fits'
    radio_optical_catalogue_fits = fits.open(fits_files+cluster_name+'_DeCALS_radio_join.fits')
    radio_optical_catalog=radio_optical_catalogue_fits[1].data

    tessel_region='A2631.tessel.reg'

    app_image=fits_files+'mypipelinerun_ABELL2631_4-MFS-image.fits'#'mypipelinerun_ABELL2631_4-MFS-image.fits'# 'image_DI_Clustered.DeeperDeconv.AP4.app.restored.fits'
    int_image=fits_files+'mypipelinerun_ABELL2631_4-MFS-image.fits'#'mypipelinerun_ABELL2631_4-MFS-image.fits'#'image_DI_Clustered.DeeperDeconv.AP4.int.restored.fits'

    rms_image=fits_files+'A2631_rms_map.fits'
    res_image=fits_files+'A2631_res_map.fits'

    fits_image=fits_files+'mypipelinerun_ABELL2631_4-MFS-image.fits'
    
    radio_HSC_catalog='A2631_HSC_radio_join.fits'

   
    cluster_centre=SkyCoord(str(354.41916666667 ), str(0.27666666666667), frame='icrs',unit=(u.deg,u.deg))

    cluster_z_lit=0.2772

    search_radius = 1.5

    R_200_Mpc=  1.93 


    plotting_scripts(output_path,fits_files,simulation_path,cluster_name,radio_catalogue_fits,radio_optical_catalog,cluster_centre,app_image,int_image,res_image,fits_image,tessel_region,rms_image,radio_HSC_catalog,COSMOS_catalogue_fits)



def combine_main( argv ):

    output_path='combine/plots/'

    fits_files='combine/input/'

    fits_image_A2631= 'combine/input/image_DD_beam_2poly_beam_brighptsrc3_mask.app.restored.fits'

    fits_image_Zwcl2341= 'combine/input/image_DD_beam_2poly_beam_brighptsrc5_mask.app.restored.fits'

    fits_image=fits_files+'image_DI_Clustered.DeeperDeconv.AP4.int.restored.fits'


    combined_NVSS_cat = fits_files+'NVSS_combine.fits'

    combined_FIRST_cat =fits_files+'FIRST_table_combined.fits'
    combined_VLASS_cat =  fits_files+ 'VLASS_combined.fits'
    
    combined_MeerKAT_cat=fits_files+'A2631_srl.fits'

    combined_MeerKAT_FIRST_cat=fits_files+'MeerKAT_FIRST_combined_20sigma.fits'#'FIRST_MeerKAT_combine.fits'
    combined_MeerKAT_VLASS_cat = fits_files+'MeerKAT_VLASS_combined_20sigma.fits'#'VLASS_MeerKAT_combine.fits'
    combined_MeerKAT_NVSS_cat=fits_files+'NVSS_MeerKAT_combine.fits'

    plotting_combine_scripts(output_path,combined_MeerKAT_cat,combined_MeerKAT_NVSS_cat,combined_MeerKAT_FIRST_cat,combined_MeerKAT_VLASS_cat,combined_VLASS_cat,combined_NVSS_cat,fits_image_A2631,fits_image_Zwcl2341)


def  plotting_scripts(output_path,fits_files,simulation_path,cluster_name,radio_catalogue_fits,radio_optical_catalog,cluster_centre,app_image,int_image,res_image,fits_image,tessel_region,rms_image,radio_HSC_catalog,COSMOS_catalogue_fits):
    
   
    #write_catalog(fits_files,cluster_name,app_image,int_image)
    #write_catalog_20sigma(fits_files,cluster_name,app_image,int_image)
    #rms_plot(output_path,cluster_name,fits_image,fits_image)
    #smearing(output_path,cluster_name,radio_catalogue_fits,cluster_centre)
    #spectral_index(output_path,cluster_name,radio_catalogue_fits
    #plot_image(output_path,cluster_name,res_image,fits_image,tessel_region)
    
    #redshift_distributions(output_path,cluster_name,radio_optical_catalog)
    
    #spectral_index_distance(output_path,cluster_name,radio_catalog,cluster_centre)
    #spectral_distributions(output_path,cluster_name,radio_catalogue_fits)
    #radio_lumo_distance(output_path,cluster_name,radio_optical_catalog,cluster_centre)
    #cumulative_rms_map(output_path,cluster_name,rms_image)
    #simulation(cluster_name,simulation_path,radio_catalogue_fits,fits_image,res_image,rms_image)
    #simulation_catalog(simulation_path)
    #completness(simulation_path,radio_catalogue_fits)
    angular_distribution(simulation_path,radio_catalogue_fits,fits_image)
    #source_counts(radio_catalogue_fits,COSMOS_catalogue_fits,output_path,cluster_name)

def  plotting_combine_scripts(output_path,combined_MeerKAT_cat,combined_MeerKAT_NVSS_cat,combined_MeerKAT_FIRST_cat,combined_MeerKAT_VLASS_cat,combined_VLASS_cat,combined_NVSS_cat,fits_image_A2631,fits_image_Zwcl2341):
    
    #astrometry_NVSS(output_path,combined_MeerKAT_NVSS_cat)
    #astrometry_VLASS(output_path,combined_MeerKAT_VLASS_cat)
    #astrometry_FIRST(output_path,combined_MeerKAT_FIRST_cat)
    #radio_cutouts(output_path,fits_image_A2631,fits_image_Zwcl2341)
    #flux_scale_VLASS(output_path,combined_MeerKAT_VLASS_cat,flux_colname='Flux',Pflux_colname='Fpk',e_flux_colname='e_Flux',survey='VLASS',freq=3e9)
    #flux_scale_FIRST(output_path,combined_MeerKAT_FIRST_cat,flux_colname='Fint',Pflux_colname='Peak_flux',e_flux_colname='E_Total_flux',survey='FIRST',freq=1.4e9)
    #flux_scale_NVSS(output_path,combined_MeerKAT_NVSS_cat,flux_colname='S1.4',Pflux_colname='Peak_flux',e_flux_colname='e_S1.4',survey='NVSS',freq=1.4e9)
    #resolved_unresolved_sources(output_path,combined_MeerKAT_cat)
    #flux_distribution(output_path,combined_MeerKAT_cat)
    #spectral_distributions(output_path,combined_MeerKAT_cat)
    #spatial_spectral_index(output_path,combined_MeerKAT_cat)

if __name__ == "__main__":
    #zwcl_main(sys.argv) 
    #A2631_main(sys.argv)
    combine_main(sys.argv)

