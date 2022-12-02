
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
