import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from matplotlib.animation import FuncAnimation
import sys
import math

#constants ---------------------
minute = 60.
hour = minute*60.
#--------------------------------

max_time_window = 10*minute # <--------------------- window duration
max_rows = int(max_time_window/2.) #1 row / 2s 
#Load last X mins/hours of data set above at the arrow

#choose file
root_dir = "/data/ADRLogs"
fname = root_dir+'/ADRLog_20260929_t133604.txt'
if len(sys.argv)>1: fname = sys.argv[1] #overwrite to input

#build fig
fig,ax_t = plt.subplots(1,1,figsize=(8,5))
ax_h = ax_t.twinx() # twin the temp axis
ax_h.set_ylabel('Heater Out [%]', color='red')        
  

fig.suptitle("ADR Temperature Monitoring")


def get_data(fname, max_rows): 
        n_rows = sum(1 for row in open(fname, 'r'))
        if(n_rows>max_rows):
                #delimiter='\t',
                df = pd.read_csv(fname,  header=None, skiprows=range(0, n_rows-max_rows))
        else:
                df = pd.read_csv(fname, header=None)
                #, delimiter='\t'
        return df



df = get_data(fname, max_rows)

def get_time(df): return (df[1]-df[1][np.argmax(df[1])])
def get_temp(df): return df[2]
def get_htro(df): return df[3]


def plot_temp(ax, df):
        ax.clear()
        ax.plot(get_time(df)/minute, get_temp(df), color='blue', label="Temp. K")
        ax.relim()
        ax.autoscale_view()
        ax.set_ylabel('Temperature [K]')
        ax.set_xlabel('Time [min]')

        
def plot_htro(ax, df):
        ax.clear()
        ax.plot(get_time(df)/minute, get_htro(df), color='red', label="H.O. %")
        ax.relim()
        ax.autoscale_view()
        ax.set_ylabel('Heater Out [%]', color='red')        
        ax.yaxis.set_label_position('right')        

        #ax.set_xlabel('Time [min]')
        #ax_h.set_ylim(0,80)

        
def calc_drms_mean(arr):
        cntr = len(arr)
        sum1 = 0.
        sum2 = 0.
        for entry in arr:
                sum1 += entry
        mean = sum1/float(cntr)
        for entry in arr:
                dif  = entry-mean
                sum2 += dif*dif
        drms = math.sqrt(sum2/float(cntr))
        return drms, mean

def print_status(df, status_prd = 60): # default is 2 mins at 30 sample / min 
        status_arr = get_temp(df)[-status_prd:]     
        drms, mean = calc_drms_mean(status_arr)
        
        top = max(get_temp(df))
        btm = min(get_temp(df))
        scale =  top - btm

        print_mean = "{:.2f}".format(mean*1e3)
        print_drms = "{:.2f}".format(drms*1e6)

        text_mean = ax_t.text(-11,top+.16*scale,"Temp. Mean: "+print_mean+" mK")
        text_drms = ax_t.text(-11,top+.10*scale,"Temp. dRMS: "+print_drms+" uK")
 
        
#Plot the temp and htro

plot_temp(ax_t, df)
plot_htro(ax_h, df)

#add a nice indicator line for setpoint
#ax_t.axhline(0.080, color='k', alpha=0.7,  label='Setpoint')


print_status(df)
fig.legend()

def update(frame):
        #update df to have only the last X mins 
        ax_t.clear()
      
        df = get_data(fname, max_rows)
        
        plot_temp(ax_t, df)
        plot_htro(ax_h, df)
 
      
        print_status(df, int(2*minute/2)) # 2 mins at 30/min
        
        fig.canvas.draw_idle()

#interval in ms
ani = FuncAnimation(fig,update, interval=2*1000, cache_frame_data=False)

plt.show()

