import pandas as pd
import numpy as np
import sys
shift = sys.argv[1]
shadow = sys.argv[2]
shad = {"0":"moon","1":"sun"}
#df=pd.read_hdf("/sps/km3net/users/fbenfe/Moon_shadow/dataframes/MC/arca8/sun/arca8rbr_mu_pre+anoise+sun.h5")
df=pd.read_hdf("/sps/km3net/users/fbenfe/Moon_shadow/dataframes/data/arca19/arca19_data_pre+anoise+"+shad[shadow]+"_fake_"+shift+"_nosun.h5")
#df = df[(df["likelihood"]>50)] #& ((df["fitinf0"])*180/np.pi<0.3)]

#df2 = df[["x","y","fitinf0","sun_dist","MC_sun_dist"]]
#df2 = df[["x","y","fitinf0","ms_dist","likelihood","fitinf10"]]
#df2 = df[["x","y","fitinf0","ms_dist","MC_ms_dist","likelihood","fitinf10"]]
df2 = df[["x","y","fitinf0",shad[shadow]+"_dist_"+shift,"likelihood","fitinf10"]]

df2.to_csv('/sps/km3net/users/fbenfe/Moon_shadow/csv/arca19/data/arca_data_pre+anoise+'+shad[shadow]+'_fake_'+shift+'_nosun.txt',sep=' ',index=False,header=0)
