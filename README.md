
######################################################
#####STEPS FOR PREPARATION OF SUN/MOON DATAFRAMES#####
######################################################

Note: divide file list in batches to speed up and avoid resources consumption 

1 --> python3 data_extractor.py 
      
      obtain "arca6rbr_mu_RUNi_RUNf.h5" : no cuts at all

2 --> on Jupyter (faster): pre+anoise cuts functions:

      obtain "arca6rbr_mu_pre+anoise_RUNi_RUNf.h5" : pre+anoise cuts

3 --> python3 extract_and_cut_correct.py 1 (1==Sun,0==Moon) (select # of time shifts)

      obtain "arca6rbr_mu_pre+anoise_RUNi_RUNf_nshifts.h5" : make time shifts and select events with reco_dist < 12 deg to Sun/Moon

4 --> python3 h1_prepare_sample.py (select Moon/Sun)

      updates "arca6rbr_mu_pre+anoise_RUNi_RUNf_nshifts.h5" : add dist. reco-MC and MC dist to sun/moon

5 --> (manually add x,y for data files) + cut only ms alt > 6 deg


######################################################
#####STEPS FOR ANALYSIS#####
######################################################


1 --> df_to_csv.py 
      
      create csv with "x","y","fitinf0","sun_dist", ( "MC_sun_dist" if MC )
      output ex: "arca6rbr_sun_allbeta_12shifts.txt"

2 --> scan_beta.cc
  
      create csv with lambda, sigma, amp scanning for betas
      otuput ex: "beta_scan_12s_shad_loglik.txt"

3 --> scan_beta_plot.cc

      plotting variables as a function of beta cuts


4 --> ./test_data_analysis.exe list_csv.txt


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
#####COMPILE ROOT SCRIPT#####
######################################################
g++ arca_2d_data_analysis.cc -o contour_data.exe `root-config --cflags --glibs`




#####################################################
#####PREPARATION OF MC SAMPLE FOR NN TRAINING#####
######################################################
I use "MC_extract_and_cut.py", but note that some df keys are removed e.g. run_id, faket...

e.g.:
./MC_extractor.sh ARCA8/MC/arca8rbr_mu_v8_list_subsample.txt /sps/km3net/users/fbenfe/Moon_shadow/CNN/arca8rbr_mu_v8.txt /sps/km3net/users/fbenfe/Moon_shadow/logs/arca8/MC/


#####################################################
#####EXTRACT AND CUT DATA ROOT FILES#####
######################################################

./extractor.sh <list_of_root_files> <name_of_output_file> <name_of_logfile>

this calls the script slurm_extractor.sh which in turn calls arca_data_extractor.py
Note: in the python script, change the directory where csv files are written according to the detid.
Also: if using the NN for classification, or if you use dst files, change the script accordingly


If you need to separate runs from a merged dst file (too big to be analyzed alone):
python dst_merged_separator.py

To analyze then the single h5 dataframes:
arca_h5_extractor.py <file> <output>

