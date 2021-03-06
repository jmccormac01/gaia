"""
Script to ingest the Split csv files from DR2

We try to use LOAD FILE so as not to have billions of
single commits to the database

An example query for loading from a file can be found here

http://www.mysqltutorial.org/import-csv-file-mysql-table/

Two columns in the table are sepecified as true/false
and need to be converted to 1/0 when ingesting. This is
done using the @dummy1 and @dummy2 entries below
"""
import os
import glob as g
import pymysql

top_dir = "/wasp/scratch/gaia_dr2"
os.chdir(top_dir)
files_to_ingest = sorted(g.glob("*.csv"))
nfiles = len(files_to_ingest)
for i, file_to_ingest in enumerate(files_to_ingest[61105:]):
    qry2 = """
            LOAD DATA LOCAL INFILE '{}'
            INTO TABLE gaia_dr2
            FIELDS TERMINATED BY ','
            LINES TERMINATED BY '\n'
            IGNORE 1 ROWS
            (solution_id,designation,source_id,random_index,ref_epoch,ra_deg,
             ra_deg_error,dec_deg,dec_deg_error,parallax,parallax_error,
             parallax_over_error,pmra,pmra_error,pmdec,pmdec_error,ra_dec_corr,
             ra_parallax_corr,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,
             dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr,
             pmra_pmdec_corr,astrometric_n_obs_al,astrometric_n_obs_ac,
             astrometric_n_good_obs_al,astrometric_n_bad_obs_al,astrometric_gof_al,
             astrometric_chi2_al,astrometric_excess_noise,astrometric_excess_noise_sig,
             astrometric_params_solved,@dummy1,astrometric_weight_al,
             astrometric_pseudo_colour,astrometric_pseudo_colour_error,
             mean_varpi_factor_al,astrometric_matched_observations,visibility_periods_used,
             astrometric_sigma5d_max,frame_rotator_object_type,matched_observations,
             @dummy2,phot_g_n_obs,phot_g_mean_flux,phot_g_mean_flux_error,
             phot_g_mean_flux_over_error,phot_g_mean_mag,phot_bp_n_obs,phot_bp_mean_flux,
             phot_bp_mean_flux_error,phot_bp_mean_flux_over_error,phot_bp_mean_mag,
             phot_rp_n_obs,phot_rp_mean_flux,phot_rp_mean_flux_error,
             phot_rp_mean_flux_over_error,phot_rp_mean_mag,phot_bp_rp_excess_factor,
             phot_proc_mode,bp_rp,bp_g,g_rp,radial_velocity,radial_velocity_error,
             rv_nb_transits,rv_template_teff,rv_template_logg,rv_template_fe_h,
             phot_variable_flag,l,b,ecl_lon,ecl_lat,priam_flags,teff_val,
             teff_percentile_lower,teff_percentile_upper,a_g_val,a_g_percentile_lower,
             a_g_percentile_upper,e_bp_min_rp_val,e_bp_min_rp_percentile_lower,
             e_bp_min_rp_percentile_upper,flame_flags,radius_val,radius_percentile_lower,
             radius_percentile_upper,lum_val,lum_percentile_lower,lum_percentile_upper)
            SET astrometric_primary_flag := @dummy1 = 'true', duplicated_source := @dummy2 = 'true';
        """.format(file_to_ingest)
    with pymysql.connect(host='ngtsdb',
                         db='catalogues',
                         local_infile=True) as cur:
        cur.execute(qry2)
    print("[{}:{}]".format(i+1, nfiles))
