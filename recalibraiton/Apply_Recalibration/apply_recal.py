#!/usr/bin/env python
# coding: utf-8

import numba
import ROOT
ROOT.EnableImplicitMT(3)
import mass
from mass.calibration.energy_calibration import EnergyCalibration
from mass.calibration.energy_calibration import EnergyCalibrationMaker
import numpy as np

linenames = ["es122kb", "es122ka", "es136kb", "es136ka", "g122", "g136"]
curvenames = ["linear", "linear+0", "loglog", "gain", "invgain", "loggain"]
curve = curvenames[3]

# --- Load in the fit result TGraphe Errors
run_no = 83
fdir = "./rc_cal_graphs/merged/" 
fname = "escal_run"+str(run_no)+".root"
path = fdir+fname
file  = ROOT.TFile(path)

# store graphs and 122 peak fit results
graph_names = {}
fit_results = {}

for key in file.GetListOfKeys(): 
    name= ROOT.TString(key.GetName())
    if(name.Contains("_cal")): 
        name_cpy = name.Data()
        num = int(name_cpy.split("_")[2][2:])
        graph_names[num] = name_cpy
        #print(name_cpy, num)
    if(name.Contains("_fr")):
        name_cpy = name.Data()
        num = int(name_cpy.split("_")[2][2:])
        fit_results[num] = name_cpy
        


# --- utility to find NANs fast
def hasnan(arr): 
    return np.isnan(np.sum(arr))

# store our fit results
cal_crvs = {}
cal_fwhm = {}
cal_u130 = {}

# loop channels
for channel in graph_names: 
    print("Channel: ", channel)
    graph = file.Get(graph_names[channel])
    v_ph = np.array(graph.GetX())
    d_ph = np.array(graph.GetEX())
    v_en = np.array(graph.GetY())
    d_en = np.array(graph.GetEY())

    #catch nans for reprocessing
    if hasnan(v_ph) or hasnan(d_ph) or hasnan(v_en) or hasnan(d_en): 
        print("Problem!")
        print(v_ph)
        print(v_en)
        print(d_ph)
        print(d_en)
        continue


    cal_maker = EnergyCalibrationMaker(v_ph,v_en,d_ph,d_en, linenames)

    #if bad stats produce a >50eV uncertainty, drop the 136kbes
    #if run_no == 83: cal_maker = cal_maker.remove_cal_point_name("es136kb") # demo 

    # Calculate the spline
    cal = cal_maker.make_calibration(curvename = curve, approximate=True)

    # Calculate the interpolation uncertainty at 130 keV
    unc130 = cal.energy2uncertainty(130000)

    # Calculate the estimated FWHM of the 122 keV line
    fit_result = file.Get(fit_results[channel])
    fwhm = 2.35*122060.65*fit_result.Parameter(2)/fit_result.Parameter(1)
    if unc130>50: 
        print(unc130, "eV", fwhm, "eV")

    # Save to dicts
    cal_crvs[channel] = cal
    cal_fwhm[channel] = fwhm
    cal_u130[channel] = unc130

# --- Load old tree
fTES = ROOT.TFile("./coadd_0708/run"+str(run_no)+"_0708.root")
tTES = fTES.Get("tTES")
#tTES.Print() # check tree


# --- Generalized tool for Set Address of all branches
from array import array

#Dict allows us to avoid typing out all var names
branches = {} # Store all variables

for branch in tTES.GetListOfBranches(): # loop old tree 
    NAME = branch.GetName()                          # Pull branch name
    #TITLE = branch.GetTitle()
    TYPE = branch.GetTitle().split("/")[-1].lower()  # Pull branch data type
    branches[NAME] = array(TYPE,[0])                 # Make a variable of the correct type in memory
    tTES.SetBranchAddress(NAME, branches[NAME])      # Set branch address call
    #print(NAME, TITLE, TYPE)



# --- Make a new tree with all the old data plus 3 new branches
fTES_new = ROOT.TFile("./recald_roots/run"+str(run_no)+"_rcal.root","RECREATE")
tTES_new = tTES.CloneTree(0) 

new_en = array("f", [0.0])
pxFWHM = array("f", [0.0])
u130 = array("f", [0.0])
b_new_en = tTES_new.Branch("new_en", new_en, "new_en/F")
b_pxfwhm = tTES_new.Branch("pxFWHM", pxFWHM, "pxFWHM/F")
b_calunc = tTES_new.Branch("u130", u130, "u130/F")
#tTES_new.Print() #check tree

# --- Loop over the existing entries ---
n_entries = tTES.GetEntries()
for i in range(n_entries):
    tTES.GetEntry(i)  # loads existing branches; access them as tree.<branch_name>

    # Set your new values here. Example using hypothetical existing branches:
    new_en[0] = cal_crvs[tTES.Chan].ph2energy(tTES.fValue_tdc)
    pxFWHM[0] = cal_fwhm[tTES.Chan]
    u130[0] = cal_u130[tTES.Chan]

    #b_new_en.Fill()


    tTES_new.Fill()
# --- Write the updated tree back, replacing the old version ---

#tTES_new.Print() #check tree

tTES_new.Write()
fTES_new.Close()
fTES.Close()




