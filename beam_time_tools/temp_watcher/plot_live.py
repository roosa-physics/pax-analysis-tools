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

#choose file
root_dir = "/data/ADRLogs"
fname = root_dir+'/ADRLog_20260929_t133604.txt'
if len(sys.argv)>1: fname = sys.argv[1] #overwrite to input

#build fig
fig,ax = plt.subplots(1,1,figsize=(10,8))
fig.suptitle("ADR Temperature Monitoring")

#add a nice indicator line for setpoint
ax.axhline(0.080, color='k', alpha=0.7,  label='Setpoint')


#Load last X mins/hours of data set above at the arrow
def get_data(fname, max_rows): 
        n_rows = sum(1 for row in open(fname, 'r'))
        if(n_rows>max_rows):
                df = pd.read_csv(fname, delimiter='\t', header=None, skiprows=range(0, n_rows-max_rows))
        else:
                df = pd.read_csv(fname, delimiter='\t', header=None)
        return df


df = get_data(fname, max_rows)

def get_time(df): return (df[1]-df[1][np.argmax(df[3])])
def get_temp(df): return df[2]
def get_htro(df): return df[3]


#initialize the plot to data
line, = ax.plot(get_time(df)/minute, get_temp(df), color='blue', label="Temp. K")

# Always lable your axes
ax.set_ylabel('Temperature [K]')
ax.set_xlabel('Time [min]')

#prep the HO line obj
axHO = ax.twinx()
lineHO, = axHO.plot(get_time(df)/minute, get_htro(df), color='red', label="H.O. %")

#Always label your axes
axHO.set_ylabel('Heater Out [%]', color='red')        
axHO.set_ylim(0,80)



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

def print_status(df, status_prd = 4*60): 
        status_arr = get_temp(df)[-status_prd:-1]     
        drms, mean = calc_drms_mean(status_arr)
        
        top = max(get_temp(df))
        btm = min(get_temp(df))
        scale =  top - btm
        print_drms = "{:.2f}".format(drms*1e6)
        print_mean = "{:.2f}".format(mean*1e3)
        ax.text(-1, top+.10*scale , "Temp. dRMS: " + print_drms + " uK")
        ax.text(-1, top+.06*scale , "Temp. Mean: " + print_mean + " mK")
        
print_status(df)
fig.legend()



def update(frame):
        #update df to have only the last X mins 
        df = get_data(fname, max_rows)
                
        line.set_data((df[1]-df[1][np.argmax(df[3])])/minute,df[2])
        ax.relim()
        ax.autoscale_view()

        lineHO.set_data((df[1]-df[1][np.argmax(df[3])])/minute,df[3])
        axHO.relim()
        axHO.autoscale_view()

        status_prd  = int(2*minute*2) # X minutes * 2 s bin
        print_status(df, status_prd)
        
        fig.canvas.draw_idle()

#interval in ms
ani = FuncAnimation(fig,update, interval=2*1000, cache_frame_data=False)

plt.show()

