"""
Compute and apply an appropriate set of automatically generated cuts.

The peak time and rise time come from the measured most-common peak time.
        The pulse RMS and postpeak-derivative cuts are based on what's observed in
        the (presumably) pulse-free noise file associated with this data file.

        Args:
            nsigma_pt_rms (float):  How big an excursion is allowed in pretrig RMS
                (default 8.0).
            nsigma_max_deriv (float): How big an excursion is allowed in max
                post-peak derivative (default 8.0).
            pretrig_rms_percentile (float): Make upper limit for pretrig_rms at
                least as large as this percentile of the data. I.e., if you
                pass in 99, then the upper limit for pretrig_rms will exclude
                no more than the 1 % largest values. This number is a
                percentage, *not* a fraction. This should not be routinely used
                - it is intended to help auto_cuts work even if there is a
                problem during a data acquisition that causes large drifts in
                noise properties.
            forceNew (bool): Whether to perform auto-cuts even if cuts already exist.
            clearCuts (bool): Whether to clear any existing cuts first (default
                True).

        The two excursion limits are given in units of equivalent sigma from the
        noise file. "Equivalent" meaning that the noise file was assessed not for
        RMS but for median absolute deviation, normalized to Gaussian distributions.

        Returns:
            The cut object that was applied.
"""

import numpy as np
import mass 

def auto_cuts(dataset,nsigma_pt_rms=8.0, nsigma_max_deriv=8.0, pretrig_rms_percentile=None, forceNew=False, clearCuts=True):

    if dataset.saved_auto_cuts is None:
        forceNew = True
    if not forceNew:
        # LOG.info("channel %g skipping auto cuts because cuts exist", dataset.channum)
        return
    
    if clearCuts:
        dataset.clear_cuts()
        
    # Step 1: peak and rise times
    if dataset.peak_samplenumber is None:
        dataset._compute_peak_samplenumber()

    MARGIN = 3  # step at least this many samples forward before cutting.
    peak_time_ms = (MARGIN + dataset.peak_samplenumber - dataset.nPresamples) * dataset.timebase * 1000
            
    # Step 2: analyze *noise* so we know how to cut on pretrig rms postpeak_deriv
    pretrigger_rms = np.zeros(dataset.noise_records.nPulses)
    for i in range(dataset.noise_records.nPulses):
        ndata = dataset.noise_records.datafile.alldata[i]
        pretrigger_rms[i] = ndata[:dataset.nPresamples].std()

    ndata_all = dataset.noise_records.datafile.alldata
    ndata_all.reshape((dataset.nSamples,dataset.noise_records.nPulses))
    max_deriv = mass.analysis_algorithms.compute_max_deriv(ndata_all, ignore_leading=0)
                
    # Multiply MAD by 1.4826 to get into terms of sigma, if distribution were Gaussian.
    md_med = np.median(max_deriv)
    pt_med = np.median(pretrigger_rms)
    md_madn = np.median(np.abs(max_deriv - md_med)) * 1.4826
    pt_madn = np.median(np.abs(pretrigger_rms - pt_med)) * 1.4826
    md_max = md_med + md_madn * nsigma_max_deriv
    pt_max = max(0.0, pt_med + pt_madn * nsigma_pt_rms)
                
    # Step 2.5: In the case of pretrig_rms, cut no more than pretrig_rms_percentile percent
    # of the pulses on the upper end. This appears to be appropriate for
    # SLEDGEHAMMER gamma devices, but may not be appropriate in cases where
    # there are many pulses riding on tails, so by default we don't do
    # this.
    if pretrig_rms_percentile is not None:
        pt_max = max(pt_max, np.percentile(self.p_pretrig_rms, pretrig_rms_percentile))
        
    # Step 3: make the cuts
    cuts = mass.core.controller.AnalysisControl(
        peak_time_ms=(0, peak_time_ms * 1.25),
        rise_time_ms=(0, peak_time_ms * 1.10),
        pretrigger_rms=(None, pt_max),
        postpeak_deriv=(None, md_max),
    )
    cuts._pretrig_rms_median = pt_med  # store these so we can acess them when writing projectors to hdf5
    cuts._pretrig_rms_sigma = pt_madn
                    
    dataset.apply_cuts(cuts, forceNew=True, clear=False)
    """Store the results of auto-cuts internally and in HDF5."""
    dataset.saved_auto_cuts = cuts
    #g = dataset.hdf5_group["cuts"].require_group("auto_cuts")
    #for attrname in ("peak_time_ms", "rise_time_ms", "pretrigger_rms", postpeak_deriv"):
    #g.attrs[attrname] = cuts.cuts_prm[attrname][1]

    return cuts
