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
    angular_separation,
    Angle
)
import astropy.units as u
from astropy.time import Time
import km3db

#remove any previous cookie
#os.remove("/pbs/home/f/fbenfe/.km3netdb_cookie")

#db = km3db.DBManager()
#Preferable choice - already in Pandas DF
#sds = km3db.StreamDS(container="pd")

start_time = time.time()

def data_extractor(filetype,filename,outputfile):

    interaction_type=[]
    #dstfile="/sps/km3net/repo/data_processing/tag/v1_ARCA19_run_by_run/data_processing/prod/data/KM3NeT_00000116/v1.1_test/reco/"+filename
    dstfile = filename
    print(dstfile)
    # at the end of the with block, the file is automatically closed
    '''
    with uproot.open(dstfile)["T"] as ifile:
    
        coords=ifile['coords']
        listdf.append(coords.arrays(library="pd"))
    
        #sumj=ifile['sum_jgandalf']  #missing in retr
        #sumj=ifile['sum_jpptrack']
        #listdf.append(sumj.arrays(library="pd"))
    
        if filetype=='MC':
            #print('This is a MC file.')
            sum_mc_evts=ifile['sum_mc_evt']
            listdf.append(sum_mc_evts.arrays(library="pd"))
    
            sumhits=ifile['sum_hits']
            sumHits=sumhits.arrays(library="pd")
            sumHits = sumHits.add_prefix('sH_')
            listdf.append(sumHits)
    
            sum_tr_hits=ifile['sum_trig_hits']
            sumTrHits=sum_tr_hits.arrays(library="pd")
            sumTrHits = sumTrHits.add_prefix('Tr_')
            listdf.append(sumTrHits)
    
            #print('T tree done!')
            dfb=pd.concat(listdf,axis=1)
            #print(dfb.keys())
            #print(len(listdf))
            #print(dfb)
    '''

    dfb = pd.DataFrame()
        
    offline=km3io.OfflineReader(dstfile)
    has_trks = offline.n_tracks[:] > 0
    #print(len(has_trks))
    
    run_id=offline.events.run_id[has_trks]
    dfb['run_id']=np.array(run_id)
        
    likelihood=offline.events.tracks.lik[has_trks][:,0]
    dfb['likelihood']=np.array(likelihood)
    
    dir_x=offline.events.tracks.dir_x[has_trks][:,0]
    dfb['dir_x']=np.array(dir_x)
    
    dir_y=offline.events.tracks.dir_y[has_trks][:,0]
    dfb['dir_y']=np.array(dir_y)
    
    dir_z=offline.events.tracks.dir_z[has_trks][:,0]
    dfb['dir_z']=np.array(dir_z)
    #shower_dir_z=offline.events.tracks.dir_z[has_trks][:,1]
    #dfb['aa_dir_z']=np.array(shower_dir_z)
    
    pos_x=offline.events.tracks.pos_x[has_trks][:,0]
    dfb['pos_x']=np.array(pos_x)
    
    pos_y=offline.events.tracks.pos_y[has_trks][:,0]
    dfb['pos_y']=np.array(pos_y)
    
    pos_z=offline.events.tracks.pos_z[has_trks][:,0]
    dfb['pos_z']=np.array(pos_z)
    
    t_sec = offline.events.t_sec[has_trks]
    dfb['t_sec'] = np.array(t_sec)
    dfb['time'] = pd.to_datetime(dfb["t_sec"], unit="s")
    
    
    E=offline.events.tracks.E[has_trks][:,0]
    dfb['Energy']=np.array(E)
    
    #print('1D array extraction done!')
    rec_stage=offline.events.tracks.rec_stages[has_trks][:,0]
    dfb['rec_stages']=ak.count(rec_stage, axis=-1)
    
    #print('len rec_stages done')
    rec_type=offline.events.tracks.rec_type[has_trks][:,0]
    dfb['rec_type']=np.array(rec_type)
    
    
    #mask=(dfb.trackfit_ra != -999.000000) & (dfb.trackfit_dec != -999.000000)
    mask=dfb.rec_type == 4000
    dfb = dfb[mask]
    fitinfs=offline.events.tracks.fitinf[has_trks][:,0]
    dfb['fitinf0']=np.array(fitinfs[mask][:,0])
    dfb['fitinf1']=np.array(fitinfs[mask][:,1])
    dfb['fitinf2']=np.array(fitinfs[mask][:,2])
    dfb['fitinf3']=np.array(fitinfs[mask][:,3])
    dfb['fitinf4']=np.array(fitinfs[mask][:,4])
    dfb['fitinf10']=np.array(fitinfs[mask][:,10])
    
    dfb=pre_cuts(dfb)
    #print("pre_cuts done")
    #dfb=logbeta0_cut(dfb)
    #print("logbeta cut done")
    dfb=anti_noise_cuts(dfb)
    #print("anoise cut done")
    dfb=get_angles(dfb)
    #print("get angles done")
    dfb=get_source_coords(dfb)
    #print("get source coord done")
    #print(dfb)
    
    #dfb=dfb.drop(["showerfit_ra","showerfit_dec","sH_nhits","sH_atot","sH_tmin","sH_tmax","sH_ndoms","sH_nlines"],axis=1)
    
    dfb=moon_sun_dist(dfb,dfb["phi"],dfb["time"],dfb["theta"],loc="arca")  #git issue, local event requires phi insted of azi
            
    print("total events = ",len(dfb))
    dfb = dfb.reset_index(drop=True)
    
    shadow = {0:"moon",1:"sun"} 
    for shad in [0,1]:
        print(shadow[shad])
        for shift in [-8,-4,0,4,8]:
            print("Cutting on shift: ",shift)
            df = distance_cut(dfb,12,shad,shift)
            df = get_relative_coordinates(df,shad,shift)
            print("events: ",len(df))
            if shift == 0:
                if len(df) != 0:
                    df.to_csv(outputfile+'_'+shadow[shad]+'.txt',columns=["x","y","fitinf0",shadow[shad]+"_dist","likelihood","fitinf10","t_sec","run_id"],sep=' ',index=False,header=0,mode='a')
            else:
                if len(df) != 0:
                    df.to_csv(outputfile+'_'+shadow[shad]+'_fake_'+str(shift)+'.txt',columns=["x","y","fitinf0",shadow[shad]+"_dist_"+str(shift),"likelihood","fitinf10","t_sec","run_id"],sep=' ',index=False,header=0,mode='a')
    return None

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

def get_relative_coordinates(df,shadow,tshift):
    """Return x=Delta_azimuth and y=Delta_altitude wrt moon/sun"""
    if tshift == 0:
        if shadow == 0:
            df["x"]=(df["local_azimuth"]*180/np.pi-df["az_moon"])*np.cos(df["local_altitude"])
            df["y"]=(df["local_altitude"]*180/np.pi-df["alt_moon"])
        elif shadow == 1:
            df["x"]=(df["local_azimuth"]*180/np.pi-df["az_sun"])*np.cos(df["local_altitude"])
            df["y"]=(df["local_altitude"]*180/np.pi-df["alt_sun"])
    else:
        if shadow == 0:
            df["x"]=(df["local_azimuth"]*180/np.pi-df["az_moon_"+str(tshift)])*np.cos(df["local_altitude"])
            df["y"]=(df["local_altitude"]*180/np.pi-df["alt_moon_"+str(tshift)])
        elif shadow == 1:
            df["x"]=(df["local_azimuth"]*180/np.pi-df["az_sun_"+str(tshift)])*np.cos(df["local_altitude"])
            df["y"]=(df["local_altitude"]*180/np.pi-df["alt_sun_"+str(tshift)])
    return df


#local_event as defined by Tamas requires zenith but the definition is wrong as it takes the zenith but then treats it like if it was theta (180-zenith)
#check it here:  https://git.km3net.de/km3py/km3astro/-/blob/master/km3astro/coord.py

def moon_sun_dist(df, azimuth, time, theta, loc):   #git issue, local event requires phi insted of azi (check phi is given when function is called!)
    """Return distance of event to moon, in detector coordinates."""
    """Also computes 4 fake sources, 4h aparts from each other"""
    evt = local_event(azimuth, time, theta, loc)
    df["local_altitude"] = evt[:].alt.rad
    df["local_azimuth"] = evt.az.rad
    #shifts = [0]
    shifts = [-8,-4,0,4,8]
    for tshift in shifts:
        print("Processing tshift = ",tshift,"...")
        moon = moon_local(time + datetime.timedelta(hours = tshift), loc)
        alt_moon = moon.alt.deg #+ 90
        az_moon = moon.az.deg
        if tshift == 0:
            df["alt_moon"]=alt_moon
            df["az_moon"]=az_moon
        else: 
            df["alt_moon_"+str(tshift)]=alt_moon
            df["az_moon_"+str(tshift)]=az_moon
        sun = sun_local(time + datetime.timedelta(hours = tshift), loc)
        alt_sun = sun.alt.deg #+ 90
        az_sun = sun.az.deg
        if tshift == 0:
            df["alt_sun"]=alt_sun
            df["az_sun"]=az_sun
        else:
            df["alt_sun_"+str(tshift)]=alt_sun
            df["az_sun_"+str(tshift)]=az_sun
        #dist_sun = evt.separation(sun)
        #dist_moon = evt.separation(moon)
        dist_sun = sep(evt,sun)
        dist_moon = sep(evt,moon)
        if tshift == 0:
            df["moon_dist"]=dist_moon.deg
            df["sun_dist"]=dist_sun.deg
            #print(df["moon_dist"])
            #print(df)
        else:
            df["moon_dist_"+str(tshift)]=dist_moon.deg
            df["sun_dist_"+str(tshift)]=dist_sun.deg
        print("Done!")
    
    return df


def sep(one,other):
    """This is the function used inside SkyCoord.separation() but with that function there is the problem with the frame transformation --> no differences for normal time but don't know why it computes different values when I calculate separations for fake sources i.e. shifted times..."""
    lon1 = one.spherical.lon
    lat1 = one.spherical.lat
    lon2 = other.spherical.lon
    lat2 = other.spherical.lat

    # Get the separation as a Quantity, convert to Angle in degrees
    sep = angular_separation(lon1, lat1, lon2, lat2)
    dist = Angle(sep, unit=u.degree)
    return dist

def get_source_coords(df):
    source_coord = km3.coord.neutrino_to_source_direction(df["phi"],df["theta"])
    azimuth = source_coord[0]
    zenith = source_coord[1]
    df["source_azimuth"] = azimuth
    df["source_zenith"] = zenith
    return df

def distance_cut(df,angle,shadow,tshift):
    if tshift == 0:
        if shadow == 0:
            dfc=df[df["alt_moon"]>0]
            dfc=dfc[dfc["moon_dist"]<angle]
        elif shadow == 1:
            dfc=df[df["alt_sun"]>0]
            dfc=dfc[dfc["sun_dist"]<angle]
    else:
        if shadow == 0:
            dfc=df[df["moon_dist_"+str(tshift)]<angle]
        elif shadow == 1:
            dfc=df[df["sun_dist_"+str(tshift)]<angle]
    return dfc

if __name__ == "__main__":
    
    file_name = sys.argv[1]
    filety='data'
    outfile = sys.argv[2]
    data_extractor(filety,file_name,outfile)
            
print("---------",time.time() - start_time," seconds ----------")
