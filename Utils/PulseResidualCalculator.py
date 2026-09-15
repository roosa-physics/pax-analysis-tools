def pulse_vals(pulse, nSamples, nPresamples, pretrigger_ignore, int_len=10e6): #mostly thieved from MASS
    e_nPresamples = nPresamples - pretrigger_ignore
    
    pretrig_sum    = 0
    pulse_sum      = 0 
    pulse_lead_sum = 0 
    pulse_rms_sum  = 0
 

    k_max = np.argmax(pulse)
    for k in range(nSamples):
        signal =  pulse[k]
        if k < e_nPresamples:
                pretrig_sum += signal
        if k >= nPresamples - 1 and k<nPresamples+int_len:
            pulse_sum += signal
            pulse_rms_sum += signal ** 2
            if k <= k_max: pulse_lead_sum += signal
 
    ptm = pretrig_sum / e_nPresamples
    
    pulse_avg = pulse_sum / (nSamples - nPresamples + 1) - ptm
    pulse_rms = np.sqrt(pulse_rms_sum / (nSamples - nPresamples + 1) - ptm * pulse_avg * 2 - ptm ** 2)
    pulse_lead_int = pulse_lead_sum - ptm*(k_max-nPresamples+1)

    return  pulse_avg, pulse_rms, pulse_lead_int