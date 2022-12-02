import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
import numpy as np
import pandas as pd
import astropy
import km3io
import matplotlib.pyplot as plt
import km3astro as km3
from km3astro.coord import local_event, sun_local, moon_local
from astropy.coordinates import (
    EarthLocation,
    SkyCoord,
    AltAz,
    Longitude,
    Latitude,
    get_sun,
    get_moon,
)

def get_relative_coordinates(df,shadow):
    shadow = int(shadow)
    """Return x=Delta_azimuth and y=Delta_altitude wrt moon/sun"""
    if shadow == 0:
        df["x"]=(df["local_azimuth"]*180/np.pi-df["az_moon"])*np.cos(df["local_altitude"])
        df["y"]=(df["local_altitude"]*180/np.pi-df["alt_moon"])
    elif shadow == 1:
        df["x"]=(df["local_azimuth"]*180/np.pi-df["az_sun"])*np.cos(df["local_altitude"])
        df["y"]=(df["local_altitude"]*180/np.pi-df["alt_sun"])
    return df

def get_relative_coordinates_fake_source(df,shadow,shift):
    """Return x=Delta_azimuth and y=Delta_altitude wrt moon/sun"""
    shadow = int(shadow)
    if shadow == 0:
        df["x"]=(df["local_azimuth"]*180/np.pi-df["az_moon_"+shift])*np.cos(df["local_altitude"])
        df["y"]=(df["local_altitude"]*180/np.pi-df["alt_moon_"+shift])
    elif shadow == 1:
        df["x"]=(df["local_azimuth"]*180/np.pi-df["az_sun_"+shift])*np.cos(df["local_altitude"])
        df["y"]=(df["local_altitude"]*180/np.pi-df["alt_sun_"+shift])
    return df

filename = sys.argv[1]
shadow = sys.argv[2]
#for fake sources files:
shift = filename.split("/")[-1].split(".")[0].split("_")[-1]

df=pd.read_hdf('/sps/km3net/users/fbenfe/Moon_shadow/dataframes/data/arca19/'+filename)
print("shift = ",shift)
print(len(df.keys()))
### SUN == 1 , MOON == 0 

#df2 = get_relative_coordinates(df,1)
#df3 = get_relative_coordinates(df,shadow)
df3 = get_relative_coordinates_fake_source(df,shadow,shift)
print(len(df3.keys()))

#df2.to_hdf('/sps/km3net/users/fbenfe/Moon_shadow/dataframes/data/arca19/arca19_data_pre+anoise+sun.h5', key='df', mode='w')
df3.to_hdf('/sps/km3net/users/fbenfe/Moon_shadow/dataframes/data/arca19/'+filename, key='df', mode='w')
#print("Total events: ",len(df2))
print("Total events: ",len(df3))
#print(len(df.keys()))
#print(df2.keys())
#df.to_hdf('/sps/km3net/users/fbenfe/Moon_shadow/arca6rbr_mu_9635-9918_h1_sun_sample.h5', key='df', mode='w')
