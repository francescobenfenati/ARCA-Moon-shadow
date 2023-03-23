Cosmic rays are blocked by nearby celestial bodies such as the Moon and the Sun. This induces a deficit in the atmospheric muon flux coming from the direction of these objects. Its observation can be used to verify the pointing accuracy and angular resolution of detectors which are able to measure the secondary muon flux from cosmic ray interactions.
This analysis is based on a comparison between the observed atmospheric muon flux coming from the Moon/Sun directions to the expected flux in case of no Moon/Sun. The search for the shadow of the Moon/Sun is done in the phase space of the angular differences between the reconstructed track coordinates and the celestial object. Both a 1-dimensional and a 2-dimensional analysis have been performed. The 1-dimensional analysis uses the space angle between the Moon/Sun and the track. The 2-dimensional analysis uses two Cartesian coordinates [x, y], starting from the zenith and azimuth angles of the sky object and the tracks.





#####################################################
#####EXTRACT AND CUT DATA ROOT FILES#####
######################################################

./extractor.sh <list_of_root_files> <name_of_output_file> <name_of_logfile>

this calls the script slurm_extractor.sh which in turn calls arca_data_extractor.py
Note: in the python script, change the directory where csv files are written according to the detid.
Also: if using the NN for classification, or if you use dst files, change the script accordingly


(Note:
If you need to separate runs from a merged dst file (too big to be analyzed as one):

python dst_merged_separator.py

To analyze then the single h5 dataframes:
arca_h5_extractor.py <file> <output>
)

#####################################################
#####PREPARATION OF MC SAMPLE FOR NN TRAINING#####
######################################################
I use "MC_extract_and_cut.py", but note that some df keys are removed e.g. run_id, faket...
usage:

./MC_extractor.sh <list_of_root_files> <name_of_output_file> <name_of_logfile>

e.g.:
./MC_extractor.sh ARCA8/MC/arca8rbr_mu_v8_list_subsample.txt /sps/km3net/users/fbenfe/Moon_shadow/CNN/arca8rbr_mu_v8.txt /sps/km3net/users/fbenfe/Moon_shadow/logs/arca8/MC/

######################################################
#####1D ANALYSIS ARCA#####
######################################################
Plot 1D hist. for ARCA Sun and Moon:
1) Prepare list of csv and plot names:
e.g. ARCA8/data/v8.0/1d_list.txt :
/sps/km3net/users/fbenfe/Moon_shadow/csv/arca21/data/arca21_grbv8_pre+anoise_sun.txt arca21_grbv8_0_sun
/sps/km3net/users/fbenfe/Moon_shadow/csv/arca21/data/arca21_grbv8_pre+anoise_moon.txt arca21_grbv8_0_moon

3) Launch:
./launcher_1D.sh <list of csv> <beta_cut> <logfile_path>
e.g.
./launcher_1D.sh ARCA8/data/v8.0/1d_list.txt 0.5 logs/arca8


######################################################
#####1D ANALYSIS ORCA#####
######################################################
Plot 1D hist. for ORCA Sun and Moon:
1) Prepare list of csv and plot names:
e.g. ORCA6/list_orca.txt :
/sps/km3net/users/fbenfe/Moon_shadow/csv/ORCA/orca_data_pre+anoise+sun.txt orca_v7_0_sun
/sps/km3net/users/fbenfe/Moon_shadow/csv/ORCA/orca_data_pre+anoise+moon.txt orca_v7_0_moon

3) Launch:
./orca_launcher_1D.sh <list of csv> <beta_cut>
e.g.
./orca_launcher_1D.sh ORCA6/list_orca.txt 0.27

######################################################
#####2D ANALYSIS ARCA#####
######################################################

3) Launch:
./launcher_2D.sh <list of csv> <sigma_fit_start> <beta_cut> <maxdeg> <bin_size> <logfile_path>
e.g.
./launcher_2D.sh ARCA8/data/v8.0/2d_list.txt 0.5 0.5 6 0.1 logs/arca8

######################################################
#####Other useful scripts#####
######################################################

1 --> scan_beta.cc

      create csv with lambda, sigma, amp scanning for betas
      otuput ex: "beta_scan_12s_shad_loglik.txt"

2 --> scan_beta_plot.cc

      plotting variables as a function of beta cuts

######################################################
#####COMPILE ROOT SCRIPT#####
######################################################
g++ data_density_histogram.cc -o data_density_histogram.exe `root-config --cflags --glibs
g++ arca_2d_data_analysis.cc -o contour_data.exe `root-config --cflags --glibs`

