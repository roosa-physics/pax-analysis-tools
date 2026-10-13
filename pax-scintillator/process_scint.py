#!/usr/bin/env python
# coding: utf-8
import lecroyparser
from lecroyparser import ScopeData

import os
import matplotlib.pyplot as plt
import pandas as pd

import calendar
import time
from datetime import timedelta

import numpy as np

def fetch_raw(trace_no, run_dir, path):
    trace_no = str(trace_no).zfill(5)
    fname = '/C1--Trace--{}.trc'.format(trace_no)
    fpath = path+run_dir+fname
    return fpath

def get_sec_date(dd):
    return float (calendar.timegm(time.strptime(dd,"%Y-%m-%d")))

def get_sec_time(tt):
    hh = int(tt.split(":")[0])
    mm = int(tt.split(":")[1])
    ss = float(tt.split(":")[2])
    return 3600*hh+60*mm+ss


def get_ts(lecroyTS):
    dd = lecroyTS.split()[0]
    tt = lecroyTS.split()[1]
    return get_sec_date(dd)+get_sec_time(tt)

def get_timedelta(stop, start=0): #define from the beginning of the run                                                   
    dt  = timedelta(seconds = get_ts(stop))
    if start==0: return dt
    else:        return dt - start


def calc_baseline(v_arr, lo = -2e-2, hi=2e-2, step=4e-5): 
    #take array of voltages and calculate the baseline with the mode                                                                                                                         
    bins = int((hi-lo)/step)
    h = TH1F("hV","hV",bins,lo,hi)
    for v in v_arr: h.Fill(v)
    return h.GetBinCenter(h.GetMaximumBin())

def calc_integral(V, i_strt, i_stop, dt=1e-10, bl=0):
    sum_bl = bl*(i_stop-i_strt)
    return (sum(V[i_strt:i_stop])-sum_bl)*dt

#lets you call as                                                                                                         
#plt.figure                                                                                                               
#plot_wvf()                                                                                                               
#more stuff                                                                                                               
#plt.show()                                                                                                               
def plot_wvf(T, V, plot_baseline=True):
    plt.plot(T, V, label='scope data')
    if plot_baseline:
        bl = calc_baseline(V)
        print('baseline = ', bl)
        plt.axhline(bl, color='r', linestyle='--', label="Baseline = "+str(bl))


#Calculate the centroid of the summed waveform
def calc_wvf_centroid(wvf_tim, wvf_sum):
    i_mid = int(len(wvf_tim)/2)
    t_mid = wvf_tim[i_mid]
    #print(i_mid, t_mid)

    for i, t in enumerate(wvf_tim): 
        if i<7000 or i>11000: continue 
        int_ft = calc_integral(wvf_sum,0,i)
        int_bk = calc_integral(wvf_sum,i,-1)
        if int_ft>int_bk:
            i_mid = i
            t_mid = t
            print(i, t)
            break
    return i_mid, t_mid


#assign the spill and MB indices for the external trigs and CR trigs

#make a function that takes an array of times for the triggers and "output arrays" for the indices
#Forces the output arrays to be 0s with len = input arr 
#allows to set a threshold, set to 1 sec by default
def calc_spill_mb_no(aTimes, fThreshold=10): 
    aSpillNo = np.zeros_like(aTimes)
    aMBNo = np.zeros_like(aTimes)
    spill_ct = 0
    for i in range(len(aTimes)): #loop thru time array 
        if aTimes[i]<0: 
            aMBNo[i] = -1 # flag for removal
            continue 
        if i==0: continue

        dt = aTimes[i] - aTimes[i-1] #calculate dt 
        if dt > fThreshold and aTimes[i-1]>0: spill_ct+=1     
        else: aMBNo[i] = aMBNo[i-1]+1
        aSpillNo[i] = spill_ct
    return aSpillNo, aMBNo


def plot_spill_profile(df, run):

    f, ax = plt.subplots(1,1,figsize=(8,6))

    ax.spines[['right', 'top']].set_visible(False)
    dts = df.get('dTime')
    spill_no, mb_no = calc_spill_mb_no(dts)

    integral = df.get('Integral')
    #good = integral > 7.5E-10 #cut 0 signals
    #plt.hist2d(mb_no[good], integral[good], 197)
    cmap = plt.cm.viridis.copy()
    cmap.set_under('white')  # Color everything below vmin white
    plt.hist2d(mb_no, integral*1E9, 197, range = [[0, 197], [0, 400]], cmap=cmap, vmin=1)
    #plt.title("run "+run+" spill profile")
    plt.xlabel(('Minibunch Index'))
    plt.ylabel(('Pulse Integral (nVs)'))

    plt.savefig("run_"+run+"_spillprofile.pdf")
    plt.show()    


#plot the summed waveform with the centroid and guesses for the integration window halfwidth
def plot_summed_wvf(wvf_tim, wvf_sum, i_mid, t_mid, run):
    hwidth = 4000

    plt.figure(figsize=[18,4])
    plot_wvf(wvf_tim, wvf_sum)

    plt.title('Accumulated wvf: run '+run+' -- centroid and bounds are given as indices')

    plt.axvline(t_mid, color='g', linestyle='--', label="Centroid = "+str(i_mid))
    plt.axvline(wvf_tim[i_mid-hwidth], color='b', linestyle='--', label="Integral bound low = "+str(i_mid-hwidth))
    plt.axvline(wvf_tim[i_mid+hwidth], color='b', linestyle='--', label="Integral bound high = "+str(i_mid+hwidth))

    #plt.yscale('log')
    plt.xlabel(('Time (s)'))
    plt.ylabel(('Voltage (V)'))
    plt.legend()
    plt.savefig("run"+run+"sumpulse.png")
    plt.show()


runlist = [rundir1, rundir2]

wvf_sum = []
wvf_tim = []
nwvfs = 0


for run in runlist:
    fdir = f'/PATH/TO/MY/ACQ/DIR/'+run
    #Start the CSV and write only the header
    #We will fill each row individually to avoid storing the full frame in memory
    outfname = 'run_'+run+'.csv'
    colnames = ['RunNo', 'FNo', 'TraceNo', 'Date', 'Time', 'dTime', 'Integral','Baseline']
    df = pd.DataFrame([], columns=colnames)
    df.to_csv(outfname, mode='w', index=False, header=True)
    td_start = 0
    for n,fname in enumerate(os.listdir(fdir)):
        if fname.endswith(".trc"): #skip non-data files

            trace_no = fname.split('--')[2][:5] #take the file number from the file name
            fpath = os.path.join(fdir, fname)
            try:
                #load the waveform
                waveform = ScopeData(fpath)

                #print(waveform)

                #parse the timestamp, get date, time and time since start
                timestamp=waveform.triggerTime

                dd = timestamp.split()[0]
                tt = timestamp.split()[1]

                if td_start == 0: td_start = get_timedelta(timestamp)
                ts = get_timedelta(timestamp, td_start)
                #print(dd, tt, dt)

                #parse the time array
                t_arr = waveform.x         #s 
                dt = t_arr[1]-t_arr[0]

                #parse the voltage array
                v_arr = -waveform.y         #V

                #Calculate the baseling and integral
                baseline  = calc_baseline(v_arr)
                integral  = calc_integral(v_arr, 0, -1, dt=dt, bl = baseline)

                #make an entry list
                entry = [run, n, trace_no, dd, tt, ts, integral, baseline]

                #write entry to the csv
                temp_df = pd.DataFrame([entry], columns=colnames)
                temp_df['dTime'] = temp_df['dTime'].dt.total_seconds()
                temp_df.to_csv(outfname, mode='a', index=False, header=False)

                #Make a summed waveform to benchmark processing params
                nwvfs+=1
                if len(wvf_sum)==0: 
                    wvf_sum = v_arr-baseline
                    wvf_tim = t_arr
                else: wvf_sum += v_arr-baseline

            except Exception as e:
                print(f"Error reading {fname}: {e}")

    df= pd.read_csv(outfname)
    df['dTime'] = pd.to_numeric(df['dTime'])
    df = df.sort_values(by='dTime', ascending=True, ignore_index=True)
    df.to_csv("run_"+run+"_newInt_sorted.csv")


    i_mid, t_mid = calc_wvf_centroid(wvf_tim, wvf_sum)
    plot_summed_wvf(wvf_tim, wvf_sum, i_mid, t_mid, run)
    plot_spill_profile(df, run)




