#some utils
import numpy as np
import matplotlib.pyplot as plt

def give_escapes(gamma, x_list):
    return (-1*np.array(x_list)+gamma).tolist()


#Make a channel-by-channel rise-time normalization using the most common rise-time
'''
for ch in data.good_channels: #loop channels
    ds = data.channel[ch] 
    ds.p_rise_time[:] = ds.p_rise_time[:]/np.median(ds.p_rise_time[:])
'''

def mode(mass_variable, mass_good_mask, no_bins, search_range): 
    c,b   = np.histogram(mass_variable[:], bins=no_bins, range=search_range)
    cg,bg = np.histogram(mass_variable[mass_good_mask], bins=no_bins, range=search_range)
    return bg[np.argmax(cg)]


def mode_norm(mass_variable, mass_good_mask, no_bins, search_range, norm = 1): 
    pk_val = mode(mass_variable, mass_good_mask, no_bins, search_range)
    if(pk_val==0): pk_val==-1
    return mass_variable[:]*norm/pk_val
    
def mode_sub(mass_variable, mass_good_mask, no_bins, search_range, scale = 1000, show_plots=False): 
    pk_val = mode(mass_variable, mass_good_mask, no_bins, search_range)
    return (mass_variable[:]-pk_val)*scale


def calc_spill_mb_no(aTimes, fThreshold=101): 
    aSpillNo = np.zeros_like(aTimes)
    aMBNo = np.zeros_like(aTimes)
    spill_ct = 0
    for i in range(1,len(aTimes)): #loop thru time array 
        #if aTimes[i]<0: 
        #    aMBNo[i] = -1 # flag for removal
        #    continue 
  
        dt = aTimes[i] - aTimes[i-1] #calculate dt 
        if dt > fThreshold: spill_ct+=1 #removed a positive time condition here     
        else: aMBNo[i] = aMBNo[i-1]+1
        aSpillNo[i] = spill_ct

    #Now correct the first extraction
    for i in range(1,200):
        if aMBNo[i]==0:          #find the start of first full spill
            offset = 198-i       #calc number of missed bunches; tot is 198 / extraction  
            for j in range(i):   #Push the offset 
                #print(mb_id[j],"to",mb_id[j]+offset)
                aMBNo[j] += offset 
            break #all done
    return aSpillNo, aMBNo


def check_ptm(data, ch, run_no):
    chan_data = data.channel[ch]
    g = chan_data.good()
    #mode(mass_variable, mass_good_mask, no_bins, search_range)
    ptm_range = (0,16000)
    bins = ptm_range[1]-ptm_range[0]
    mode_ptm = mode(chan_data.p_pretrig_mean, g, bins,ptm_range)

    start = chan_data.p_timestamp[0]
    

    plt.figure(figsize=(7,3.5))
    plt.suptitle("run "+str(run_no)+" channel "+str(ch)+" ptm check", fontsize=16)
    plt.plot((chan_data.p_timestamp[:] - start)/60, chan_data.p_pretrig_mean[:]-mode_ptm, '.',ms=0.5)
    plt.plot((chan_data.p_timestamp[g] - start)/60, chan_data.p_pretrig_mean[g]-mode_ptm, '.',ms=0.5)
    plt.ylabel("dPTM from mode (daq U)")
    plt.xlabel("Time from start (min)")
    #plt.ylim(ptm_range)

def check_dc(data, ch, run_no):
    ds = data.channel[ch]
    g = ds.good()
    mode_ptm    = mode(ds.p_pretrig_mean,    g, 12000, (0,12000))
    mode_fv     = mode(ds.p_filt_value,      g, 12000, (0,12000))    
    mode_fv_dc  = mode(ds.p_filt_value_dc,   g, 12000, (0,12000))    

    #dc_check_xlim = (-800,800)
    dc_check_ylim = (-30,30)  

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12, 4))
    
    ax1.plot(ds.p_pretrig_mean[:]-mode_ptm, ds.p_filt_value[:]-mode_fv, '.',ms=0.5)
    ax1.plot(ds.p_pretrig_mean[g]-mode_ptm, ds.p_filt_value[g]-mode_fv,'.',ms=0.5)
    ax1.set_xlabel('Pretrigger mean (arb)')
    ax1.set_ylabel('filt. value (arb)')
    ax1.set_ylim(dc_check_ylim)
    #ax1.set_xlim(dc_check_xlim)
    
    
    ax2.plot(ds.p_pretrig_mean[:]-mode_ptm, ds.p_filt_value_dc[:]-mode_fv_dc, '.',ms=0.5)
    ax2.plot(ds.p_pretrig_mean[g]-mode_ptm, ds.p_filt_value_dc[g]-mode_fv_dc,'.',ms=0.5)
    ax2.set_xlabel('Pretrigger mean (arb)')
    ax2.set_ylabel('filt. value drift corrected (arb)')
    ax2.set_ylim(dc_check_ylim)
    #ax2.set_xlim(dc_check_xlim)
   
    
    plt.suptitle("run "+str(run_no)+" channel "+str(ch)+" drift correct check", fontsize=16)
