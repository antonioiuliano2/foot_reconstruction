from __future__ import division
import pandas as pd
import ROOT as r
from scipy import stats
import fedra_utils as feutils
import numpy as np
import rootnumpy_myutils as myrootnp
import sys
'''launch with ipython -i trackinquality.py nsection'''
nsection = sys.argv[1] #1,2, 3, 4,5,6,7

print("Processing Data Frame for section {}".format(nsection))
df = pd.read_csv("GSI1_S{}.csv".format(nsection))
simdf = df.query("MCTrack>=0 or TrackID>=0") #segments tracked or from the simulation
simdf["Theta"] = np.arctan(np.sqrt(simdf["TX"] * simdf["TX"] + simdf["TY"] * simdf["TY"]))

nseg = simdf.groupby("TrackID").count()["ID"]
nsegsamemc = simdf.groupby(["TrackID","MCEvent","MCTrack"]).count()["PID"] #associated to the true MC track
#taking tracked segments
trackdf = simdf.query("TrackID>=0")
#For each track, take the first segment
trackdf = simdf.groupby("TrackID").last()
trackdf["nseg"] = nseg

trackdf = trackdf.query("MCEvent>=0")

maxnseg = int(np.max(nseg))
#group by same MCTrack and MCEvent
simdf = simdf.query("MCEvent>=0")
nsegtrue = simdf.groupby(["MCEvent","MCTrack"]).count()["ID"] #true length of the Monte Carlo track
maxnsegtrue = int(np.max(nsegtrue))
#merging this information in trackID

trackdf = trackdf.reset_index().merge(nsegtrue.reset_index(), how = "left", on = ["MCEvent","MCTrack"])
trackdf = trackdf.reset_index().merge(nsegsamemc.reset_index(), how = "left", on = ["TrackID","MCEvent","MCTrack"])

#renaming labels
trackdf["PID"] = trackdf["PID_x"]
trackdf["ID"] = trackdf["ID_x"]
trackdf["nsegsamemc"] = trackdf["PID_y"]
trackdf["nsegtrue"] = trackdf["ID_y"]

del trackdf["ID_y"]
del trackdf["ID_x"]
del trackdf["PID_y"]
del trackdf["PID_x"]

print("end preparation of dataframe, making plots")

#computing efficiency
trackdf["efficiency"] = trackdf["nsegsamemc"]/trackdf["nsegtrue"]

#how many reco tracks for each true track?
nsplit = trackdf.groupby(["MCEvent","MCTrack"]).count()["ID"]

#making histograms, drawing them
hsplit = r.TH1I("hsplit","number of reconstructed tracks for each true track;Nsplit",10,1,10)
heff = r.TH1D("heff","track efficiency;eff",11,0,1.1)

hefflength_truenseg = r.TProfile("hefflength_truenseg","Track efficiency vs true number of segments;Nseg;eff",maxnsegtrue,0,maxnsegtrue,0,1.1)
hefflength_nseg = r.TProfile("hefflength_nseg","Track efficiency vs number of segments;Nseg;eff",maxnseg,0,maxnseg,0,1.1)

heffangle = r.TProfile("heffangle","Track efficiency vs angle;#theta[rad];eff",10,0.,1,0,1.1)

csplit = r.TCanvas()
myrootnp.fillhist1D(hsplit,nsplit)
hsplit.Draw()

ceff = r.TCanvas()
myrootnp.fillhist1D(heff,trackdf["efficiency"])
heff.Draw()

cefflengthtrue = r.TCanvas()
myrootnp.fillprofile2D(hefflength_truenseg,trackdf["nsegtrue"],trackdf["efficiency"])
hefflength_truenseg.Draw()

cefflength = r.TCanvas()
myrootnp.fillprofile2D(hefflength_nseg,trackdf["nseg"],trackdf["efficiency"])
hefflength_nseg.Draw()

ceffangle = r.TCanvas()
myrootnp.fillprofile2D(heffangle,trackdf["Theta"],trackdf["efficiency"])
heffangle.Draw()
#csplitlength = r.TCanvas()
#myrootnp.fillhist2D(hsplitlength,trackdf["nsegtrue"],nsplit)
#hsplitlength.Draw("COLZ")

#printing trackdf

print(trackdf[["TrackID","MCEvent","MCTrack","nsegtrue","efficiency"]])
