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
db = km3db.DBManager()
#Preferable choice - already in Pandas DF
sds = km3db.StreamDS(container="pd")



def data_extractor(filetype,filename):

    dstfile=filename
    # at the end of the with block, the file is automatically closed
    #remove for ORCA
    
    listdf=[]
    with uproot.open(dstfile)["T"] as ifile:
        #print(ifile.keys())
        coords=ifile['coords']
        listdf.append(coords.arrays(library="pd"))

        #sumj=ifile['sum_jgandalf'] #<------- this is not present in MC of bad period runs!
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
    
    #dfb = pd.DataFrame()
    offline=km3io.OfflineReader(dstfile)
    has_trks = offline.n_tracks[:] > 0
    print(len(has_trks))
    print(offline.events.tracks.fitinf)

    run_id=offline.events.run_id[has_trks]
    dfb['run_id']=np.array(run_id)
    #print(run_id)
    #print(dfb["run_id"])

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
    
    E=offline.events.tracks.E[has_trks][:,0]
    dfb['Energy']=np.array(E)
    
    #print('1D array extraction done!')
    rec_stage=offline.events.tracks.rec_stages[has_trks][:,0]
    dfb['rec_stages']=ak.count(rec_stage, axis=-1)
    
    #print('len rec_stages done')
    rec_type=offline.events.tracks.rec_type[has_trks][:,0]
    dfb['rec_type']=np.array(rec_type)
    
    if filetype=='MC':
        #print('This is a MC file.')

        mc_trk_type=offline.events.mc_trks.pdgid[has_trks][:,0]
        dfb['MC_trk_type']=np.array(mc_trk_type)
        
        mce=offline.events.mc_trks.E[has_trks][:,0]
        #print(offline.events.mc_trks.keys())
        dfb['MC_E']=np.array(mce)
        
        
        MC_dir_x=offline.events.mc_tracks.dir_x[has_trks][:,0]
        dfb['MC_dir_x']=np.array(MC_dir_x)
        
        MC_dir_y=offline.events.mc_tracks.dir_y[has_trks][:,0]
        dfb['MC_dir_y']=np.array(MC_dir_y)
        
        MC_dir_z=offline.events.mc_tracks.dir_z[has_trks][:,0]
        dfb['MC_dir_z']=np.array(MC_dir_z)
        
        angular_resolution = np.arccos(dir_x*MC_dir_x+dir_y*MC_dir_y+dir_z*MC_dir_z)*180/np.pi
        dfb['angular_resolution'] = angular_resolution
        
    
    #mask=(dfb.trackfit_ra != -999.000000) & (dfb.trackfit_dec != -999.000000)
    mask=dfb.rec_type == 4000
    dfb = dfb[mask]
    fitinfs=offline.events.tracks.fitinf[has_trks][:,0]
    print(type(fitinfs))
    print(len(fitinfs[mask]))
    #print("fitinfs = ",fitinfs[mask][0])
    #print("fitinfs = ",fitinfs[mask][:,21])
    #print(len(fitinfs[mask][:,21]))
    for i in range(23):
        dfb['fitinf'+str(i)]=np.array(fitinfs[mask][:,i])
        print("fitinf",str(i)," = ",fitinfs[mask][:,i])
    '''
    dfb['fitinf0']=np.array(fitinfs[mask][:,0])
    dfb['fitinf1']=np.array(fitinfs[mask][:,1])
    dfb['fitinf2']=np.array(fitinfs[mask][:,2])
    dfb['fitinf3']=np.array(fitinfs[mask][:,3])
    dfb['fitinf4']=np.array(fitinfs[mask][:,4])
    dfb['fitinf10']=np.array(fitinfs[mask][:,10])
    '''

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
    
    '''
    detector_name = "D0ARCA021"
    #Generate uniform times
    run_start=int(str(sds.runs(detid=detector_name,run=run_id[0])['UNIXSTARTTIME'][0])[:-3])
    start_time=int(str(sds.runs(detid=detector_name,run=run_id[0])['UNIXJOBSTART'][0])[:-3])
    end_time =int(str(sds.runs(detid=detector_name,run=run_id[0])['UNIXJOBEND'][0])[:-3])
    
    #generate a fake time for each muon event, randomly between start date and end date
    start_date = pd.to_datetime(start_time, unit="s")
    end_date = pd.to_datetime(end_time, unit="s")
    
    seconds_between_dates = end_time-start_time
    times = []
    
    for evt in range(len(dfb)):
        random_number_of_time = random.randrange(seconds_between_dates)
        random_date = pd.to_datetime(start_time+random_number_of_time, unit="s")
        times.append(random_date)

    dfb['faket']=np.array(times)
    '''
    #not for ORCA
    #dfb=dfb.drop(["mjd","nu_ra","nu_dec","MC_run","weight","weight_noOsc"],axis=1)
    dfb = dfb.drop(["mjd","nu_ra","nu_dec","MC_run","weight","weight_noOsc","t_sec","MC_trk_type",'t_sec','run_id','livetime_DAQ','livetime_sim','MC_E','MC_dir_x', 'MC_dir_y', 'MC_dir_z','MC_theta', 'MC_alt','MC_phi','fitinf8','fitinf11','fitinf12','fitinf17'],axis=1)
    dfb = dfb.drop(['fitinf5','fitinf6','fitinf7','fitinf9','fitinf13','fitinf14','fitinf15','fitinf16','fitinf18','fitinf19','fitinf20','fitinf21','fitinf22'],axis=1)
    '''
    #making shifts
    #shifts=np.linspace(0,24*3600,12,endpoint=False)
    if shadow == "0":
        listofshifts=[]
        print("Moon selected")
        for i in range(52): #<---------- n of shifts !
            cp = dfb.copy()
            cp["faket"] = cp["faket"] + datetime.timedelta(hours = 1*i)
            #cp["faket"] = cp["faket"] + datetime.timedelta(hours = i)
            cp=moon_distance_cut(cp,12)
            cp["shift"]=np.full(len(cp),i)
            print("Events in shift ",i,": ",len(cp))
            listofshifts.append(cp)
        dfb=pd.concat(listofshifts)
        dfb = dfb.reset_index(drop=True)
    elif shadow == "1":
        listofshifts=[]
        print("Sun selected")
        for i in range(52):
            cp = dfb.copy()
            cp["faket"] = cp["faket"] + datetime.timedelta(hours = 1*i)
            #    cp["faket"] = cp["faket"] + datetime.timedelta(hours = i)
            cp=sun_distance_cut(cp,12)
            cp["shift"]=np.full(len(cp),i)
            print("Events in shift ",i,": ",len(cp))
            listofshifts.append(cp)
        dfb=pd.concat(listofshifts)
        dfb = dfb.reset_index(drop=True)
        
        #if shadow == "0":
        #    dfb=moon_distance_cut(dfb,12)
        #elif shadow == "1":
        #    dfb=sun_distance_cut(dfb,12)
        
        '''

    print("Events after cuts: ",len(dfb))
    dfb = dfb.reset_index(drop=True)
    #print(result)

    return dfb

def get_angles(df):
    #RECO:
    #THETA muon and THETA source
    cos_theta = list(df["dir_z"])
    theta=[]
    theta_inv =[] #this is the "real" theta from zenith to source direction (180 - theta_neutrino)
    for i in cos_theta:
        theta_inv.append((math.pi - (math.acos(i))))
        theta.append(math.acos(i))
        #cos_theta.append(-1*i)
    theta=np.array(theta)
    theta_inv=np.array(theta_inv)
    df["theta"]=theta
    #df["theta_inv"]=theta_inv

    #Monte Carlo:
    #THETA muon and THETA source
    cos_theta_mc = list(df["MC_dir_z"])
    theta_mc=[]
    theta_inv_mc =[] #this is the "real" theta from zenith to source direction (180 - theta_neutrino)
    for i in cos_theta_mc:
        theta_inv_mc.append((math.pi - (math.acos(i))))#*360/(2*math.pi))
        theta_mc.append(math.acos(i))#*360/(2*math.pi))
        #cos_theta_mc.append(-1*i)
    theta_mc=np.array(theta_mc)
    theta_inv_mc=np.array(theta_inv_mc)
    df["MC_theta"]=theta_mc
    #df["MC_theta_inv"]=theta_inv_mc

    
    #RECO:
    #ALTITUDE
    alt=[]
    for i in cos_theta:
        alt.append(math.acos(i)-math.pi/2)#*360/(2*math.pi)-90)
    df["alt"] = np.array(alt)
            

    #Monte Carlo:
    #ALTITUDE
    alt_mc=[]
    for i in cos_theta_mc:
        alt_mc.append(math.acos(i)-math.pi/2)#*360/(2*math.pi)-90)
    df["MC_alt"] = np.array(alt_mc)
    

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


    #Monte Carlo:
    #PHI muon track
    #y_mc = list(df["MC_dir_y"])
    #x_mc = list(df["MC_dir_x"])
    #phi_mc = []
    #for i in range(len(df["MC_dir_y"])):
    #    phi_mc.append(math.atan2(y_mc[i],x_mc[i]))#*360/(2*math.pi))
    #phi_mc=np.array(phi_mc)
    
    phi_mc = kp.math.phi(np.array([df["MC_dir_x"], df["MC_dir_y"], df["MC_dir_z"]]).T)
    df["MC_phi"] = phi_mc


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

def moon_dist(df, azimuth, time, theta, loc): 
    """Return distance of event to moon, in detector coordinates."""
    evt = local_event(azimuth, time, theta, loc)
    df["local_altitude"] = evt[:].alt.rad
    df["local_azimuth"] = evt.az.rad
    moon = moon_local(time, loc)
    alt_moon = moon.alt.deg #+ 90
    az_moon = moon.az.deg
    df["alt_moon"]=alt_moon
    df["az_moon"]=az_moon
    dist = evt.separation(moon)
    return dist


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

def get_source_coords(df):
    source_coord = km3.coord.neutrino_to_source_direction(df["phi"],df["theta"])
    azimuth = source_coord[0]
    zenith = source_coord[1]
    df["source_azimuth"] = azimuth
    df["source_zenith"] = zenith
    return df

def moon_distance_cut(df,angle):
    dist_to = moon_dist(df, df["phi"], df["faket"], df["theta"], loc="arca") #git issue, local event requires phi insted of azi
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


if __name__ == "__main__":

    filety='MC'
    
    #Shadow ==> Moon = 0, Sun = 1 
    file_name = sys.argv[1]
    outfile = sys.argv[2]
    print(sys.argv[1])
    

    data=data_extractor(filety,file_name)
    print(data.keys())
    print(len(data.keys()))
    data.to_csv(outfile, sep=',', header=False, index=False, mode='a')
    #data.to_hdf("/sps/km3net/users/fbenfe/Moon_shadow/dataframes/MC/arca8/sun/arca8rbr_mu_pre+anoise+sun_52shifts.h5", key='df', mode='w')
    
