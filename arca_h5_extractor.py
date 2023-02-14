import datetime
import tracemalloc
import random
import uproot
import math
import time
import numpy as np
import km3io
import km3pipe as kp
import sys
import pandas as pd
from matplotlib import pyplot as plt
import os
from tqdm import tqdm
import re
import awkward as ak
import utm
import scipy.stats as st
#import aa
#from scipy.stats import poisson
#import random
import km3astro as km3
from km3astro.coord import local_event, sun_local, moon_local
from astropy.coordinates import (
    concatenate,
    FK5,
    EarthLocation,
    SkyCoord,
    AltAz,
    Longitude,
    Latitude,
    get_sun,
    get_moon,
)
import astropy.units as u
from astropy.time import Time
import km3db
#db = km3db.DBManager()
#Preferable choice - already in Pandas DF
#sds = km3db.StreamDS(container="pd")

start_time = time.time()

def data_extractor(filetype,filename):

    dfb=pd.read_hdf(filename)

    print(dfb)

    dfb=moon_sun_dist(dfb,dfb["phi"],dfb["time"],dfb["theta"],loc="arca")  #git issue, local event requires phi insted of azi
    print("moon/sun dist calculated")
    print(dfb["moon_dist"].min())
    print(dfb["sun_dist"].min())
    

    dfb_sun,dfb_moon=cut_ms_dist(dfb)
    print("moon/sun df created")
    dfb_sun=get_relative_coordinate(dfb_sun,1)
    dfb_moon=get_relative_coordinate(dfb_moon,0)
    print("x y coords. done")
    
    dfb = dfb.reset_index(drop=True)    
    dfb_sun = dfb_sun.reset_index(drop=True)
    dfb_moon = dfb_sun.reset_index(drop=True)

    print(dfb_sun.keys())
    print(len(dfb_sun.keys()))

    return dfb,dfb_sun,dfb_moon


def get_angles(df):
    #RECO:
    #THETA muon and THETA source
    cos_theta = list(df["dir_z"])
    theta=[]
    theta_inv =[] #this is the "real" theta from zenith to source direction (180 - theta_neutrino)
    for i in cos_theta:
        theta_inv.append((math.pi - (math.acos(i))))#*360/(2*math.pi))
        theta.append(math.acos(i))#*360/(2*math.pi))
        #cos_theta.append(-1*i)
    theta=np.array(theta)
    theta_inv=np.array(theta_inv)
    df["theta"]=theta
    df["theta_inv"]=theta_inv

    
    #RECO:
    #ALTITUDE
    alt=[]
    for i in cos_theta:
        alt.append(math.acos(i)-math.pi/2)#*360/(2*math.pi)-90)
    df["alt"] = np.array(alt)
            

    #RECO:
    #PHI muon track
    #y = list(df["dir_y"])
    #x = list(df["dir_x"])
    #phi = []
    #for i in range(len(df["dir_y"])):
    #    phi.append(math.atan2(y[i],x[i]))#*360/(2*math.pi))
    #phi=np.array(phi)
    phi = kp.math.phi(np.array([df["dir_x"], df["dir_y"], df["dir_z"]]).T)    
    df["phi"] = phi

    return df
                
def pre_cuts(df):
    df1=df[(df['rec_type']==4000) & (df['dir_z']<0) & (df['rec_stages']>4) & (df['likelihood']>0) & (df['fitinf10']>0) & (df['fitinf0']>0)]
    return df1

def anti_noise_cuts(df):
    df2=df[(df['fitinf3']>15) & ((df['likelihood']/df["fitinf3"])>0.5)] #Luc cut
    return df2

def logbeta0_cut(df):
    df3=df[np.log10(np.array(df['fitinf0'])*180/np.pi)<0]
    return df3

def Energy_cut(df,Energy_value):
    df4=df[np.log10(np.array(df['Energy']))>np.log10(Energy_value)]
    return df4

def Direction_cut(df,deg_dir_z_cut):
    df5=df[(180/np.pi)*np.arccos(df['dir_z'])<deg_dir_z_cut]
    return df5

#local_event as defined by Tamas requires zenith but the definition is wrong as it takes the zenith but then treats it like if it was theta (180-zenith)
#check it here:  https://git.km3net.de/km3py/km3astro/-/blob/master/km3astro/coord.py

def moon_sun_dist(df, azimuth, time, theta, loc):   #git issue, local event requires phi insted of azi
    """Return distance of event to moon, in detector coordinates."""
    evt = local_event(azimuth, time, theta, loc)
    df["local_altitude"] = evt[:].alt.rad
    df["local_azimuth"] = evt.az.rad
    moon = moon_local(time, loc)
    alt_moon = moon.alt.deg #+ 90
    az_moon = moon.az.deg
    df["alt_moon"]=alt_moon
    df["az_moon"]=az_moon
    sun = sun_local(time, loc)
    alt_sun = sun.alt.deg #+ 90
    az_sun = sun.az.deg
    df["alt_sun"]=alt_sun
    df["az_sun"]=az_sun
    dist_sun = evt.separation(sun)
    dist_moon = evt.separation(moon)
    df["moon_dist"]=dist_moon.deg
    df["sun_dist"]=dist_sun.deg
    return df

def get_relative_coordinate(df,shadow):
    shadow = int(shadow)
    """Return x=Delta_azimuth and y=Delta_altitude wrt moon/sun"""
    if shadow == 0:
        df["x"]=(df["local_azimuth"]*180/np.pi-df["az_moon"])*np.cos(df["local_altitude"])
        df["y"]=(df["local_altitude"]*180/np.pi-df["alt_moon"])
    elif shadow == 1:
        df["x"]=(df["local_azimuth"]*180/np.pi-df["az_sun"])*np.cos(df["local_altitude"])
        df["y"]=(df["local_altitude"]*180/np.pi-df["alt_sun"])
    return df

def sun_dist(df, azimuth, time, theta, loc):
    """Return distance of event to sun, in detector coordinates."""
    evt = local_event(azimuth, time, theta, loc)
    df["local_altitude"] = evt[:].alt.rad
    df["local_azimuth"] = evt.az.rad
    sun = sun_local(time, loc)
    alt_sun = sun.alt.deg #+ 90
    az_sun = sun.az.deg
    df["alt_sun"]=alt_sun
    df["az_sun"]=az_sun    
    dist = evt.separation(sun)
    return dist

def cut_ms_dist(df):
    df_sun=df[df["sun_dist"]<12]
    df_sun=df_sun[df_sun["alt_sun"]>0]
    df_moon=df[df["moon_dist"]<12]
    df_moon=df_moon[df_moon["alt_moon"]>0]
    return df_sun,df_moon

def get_source_coords(df):
    source_coord = km3.coord.neutrino_to_source_direction(df["phi"],df["theta"])
    azimuth = source_coord[0]
    zenith = source_coord[1]
    df["source_azimuth"] = azimuth
    df["source_zenith"] = zenith
    return df

'''
def moon_distance_cut(df,angle):
    dist_to = moon_dist(df, df["phi"], df["faket"], df["theta"], loc="arca")  #git issue, local event requires phi insted of azi
    df["moon_dist"]=dist_to.deg
    df6=df[df["alt_moon"]>0]
    df6=df6[df6["moon_dist"]<angle]
    return df6

def sun_distance_cut(df,angle):
    dist_to = sun_dist(df, df["phi"], df["faket"], df["theta"], loc="arca")  #git issue, local event requires phi insted of azi
    df["sun_dist"]=dist_to.deg
    df6=df[df["alt_sun"]>0]
    df6=df6[df6["sun_dist"]<angle]
    return df6
'''

if __name__ == "__main__":

    filety='data'
    filename = sys.argv[1]
    outfile = sys.argv[2]

    data,data_sun,data_moon = data_extractor(filety,filename)
    #data = data_extractor(filety,filename)
    
    data_sun.to_csv('dataframes/data/arca21/'+outfile+'_sun.txt',mode='a',sep=' ',index=False,header=0)
    data_moon.to_csv('dataframes/data/arca21/'+outfile+'_moon.txt',mode='a',sep=' ',index=False,header=0)

    data_sun = data_sun[["x","y","fitinf0","sun_dist","likelihood","fitinf10","t_sec"]]
    data_moon = data_moon[["x","y","fitinf0","moon_dist","likelihood","fitinf10","t_sec"]]

    data_sun.to_csv('csv/arca21/data/'+outfile+'_sun.txt',mode='a',sep=' ',index=False,header=0)
    data_moon.to_csv('csv/arca21/data/'+outfile+'_moon.txt',mode='a',sep=' ',index=False,header=0)
    
    
    #save all events, no cut near moon/sun
    data.to_csv('csv/arca21/data/'+outfile+'.txt', sep=',', header=False, index=False, mode='a')   

    print("---------",time.time() - start_time," seconds ----------")
