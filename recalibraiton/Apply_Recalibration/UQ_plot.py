#!/usr/bin/env python
# coding: utf-8
import numba
import ROOT
ROOT.EnableImplicitMT(3)
import mass
from mass.calibration.energy_calibration import EnergyCalibration
from mass.calibration.energy_calibration import EnergyCalibrationMaker
import numpy as np
import matplotlib.pyplot as plt

linenames = ["es122kb", "es122ka", "es136kb", "es136ka", "g122", "g136"]
curvenames = ["linear", "linear+0", "loglog", "gain", "invgain", "loggain"]
curve = curvenames[3]

fdir = "./wkdir/"
fname = "fname.root"
path = fdir+fname

file  = ROOT.TFile(path)

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
#print(graph_names, fit_results)

def hasnan(arr): 
    return np.isnan(np.sum(arr))


cal_crvs = {}
unc_crvs = {}


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

    #if bad stats produce a >XXeV uncertainty, drop the 
    #if run_no == 99 and chan == 99: cal_maker = cal_maker.remove_cal_point_name("badpoint")

    cal_crvs[channel] = cal_maker.make_calibration(curvename = curve, approximate=True)


def getslope(x_axis, y_axis):
    xi, xf = x_axis[0], x_axis[-1]
    yi, yf = y_axis[0], y_axis[-1]
    dx = xf-xi
    dy = yf-yi 
    m = float(dy/dx)
    slope = [x*m for x in x_axis]
    slope += yi - m*xi

    return slope

f, ax = plt.subplots(1,1,figsize=(16,4))

erange  = range(95000, 140000)

for i,ch in enumerate(cal_crvs):
    if ch > 2 : continue
    cal = cal_crvs[ch].energy2ph(erange)
    unc = cal_crvs[ch].energy2uncertainty(erange)

    slope = getslope(cal, erange)
    plt.plot(cal, erange-slope, label = "ch "+str(ch))
    plt.fill_between(cal, erange-slope+unc, erange-slope-unc, alpha = .1)


plt.title("run "+str(run_no)+" ch "+" E2P Curves")
plt.xlabel(('Pulse Height Units'))
plt.ylabel(('Energy Diff. From Linear (eV)'))
ax.spines[['right', 'top']].set_visible(False)
#plt.xlim(120e3,140e3)
#plt.ylim(-.5,.5)
plt.legend()
plt.savefig("rc"+str(run_no)+"_calcurve_ex.png")
