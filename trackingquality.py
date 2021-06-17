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

setfile = r.TFile.Open("b000001.{}.0.0.set.root".format(nsection))
footset = setfile.Get("set")

#def getPlatebyPID(df):
# '''getting plate by pID, according to set'''
# return footset.GetID(int(df.PID)).ePlate

lastplateset = footset.GetID(0).ePlate
print("Opened set of section {}, last plate is {}".format(nsection,lastplateset))


print("Processing Data Frame for section {}".format(nsection))
df = pd.read_csv("GSI1_S{}.csv".format(nsection))
#df["Plate"] = df.apply(getPlatebyPID,axis=1) #getting plate by pID

simdf = df.query("MCTrack>=0 or TrackID>=0") #segments tracked or from the simulation
simdf["Theta"] = np.arctan(np.sqrt(simdf["TX"] * simdf["TX"] + simdf["TY"] * simdf["TY"]))

nseg = simdf.groupby("TrackID").count()["ID"]
nsegsamemc = simdf.groupby(["TrackID","MCEvent","MCTrack"]).count()["PID"] #associated to the true MC track
#taking tracked segments
trackdf = simdf.query("TrackID>=0")
#computing npl
PIDlast = trackdf.groupby("TrackID").min()["PID"]
PIDfirst = trackdf.groupby("TrackID").max()["PID"]
npl = (PIDfirst - PIDlast) + 1

#For each track, take the last segment, accept them if they have at least one segment Monte Carlo
trackdf["simulation"] = trackdf["MCEvent"]>=0 
atleastonemc = trackdf.groupby("TrackID").any()["simulation"] #at least one segment coming from simulation

#which are the most frequent MonteCarlo Event and Track for this reconstructed track?
#mostfrequentevent = trackdf.query("MCEvent>=0").groupby(['TrackID'])['MCEvent'].agg(pd.Series.mode)
#mostfrequenttrack = trackdf.query("MCTrack>=0").groupby(['TrackID'])['MCTrack'].agg(pd.Series.mode) #mode is not ok, returns only values present at least twice
mostfrequentevent = trackdf.query("MCEvent>=0").groupby(['TrackID'])['MCEvent'].agg(lambda x:x.value_counts().index[0])
mostfrequenttrack = trackdf.query("MCEvent>=0").groupby(['TrackID'])['MCTrack'].agg(lambda x:x.value_counts().index[0])

trackdf = trackdf.groupby("TrackID").last()

#adding computed information to our dataset
trackdf["nseg"] = nseg
trackdf["npl"] = npl
trackdf["fedraeff"] = nseg/npl
trackdf["MostCommonMCEvent"] = mostfrequentevent
trackdf["MostCommonMCTrack"] = mostfrequenttrack

trackdf = trackdf[atleastonemc]

maxnseg = int(np.max(nseg))
#group by same MCTrack and MCEvent
simdf = simdf.query("MCEvent>=0")
nsegtrue = simdf.groupby(["MCEvent","MCTrack"]).count()["ID"] #true length of the Monte Carlo track
maxnsegtrue = int(np.max(nsegtrue))

#let us consider now signal segments which have been tracked
simtrackdf = simdf.query("TrackID>=0")
#grouping by MCEvent and MCTrack, we can get where the track ends (even counting splitting in different reco tracks)
simtrackdf = simtrackdf.sort_values("PID")
lastplateMCTrack = simtrackdf.groupby(["MCEvent","MCTrack"]).first()["PID"]

#merging the MCTrue information in trackID, according to most common MCEvent and MCTrack
nsegtrue = nsegtrue.reset_index()
nsegsamemc = nsegsamemc.reset_index()
lastplateMCTrack = lastplateMCTrack.reset_index()

nsegtrue["MostCommonMCEvent"] = nsegtrue["MCEvent"]
nsegtrue["MostCommonMCTrack"] = nsegtrue["MCTrack"]

del nsegtrue["MCEvent"]
del nsegtrue["MCTrack"]

nsegsamemc["MostCommonMCEvent"] = nsegsamemc["MCEvent"]
nsegsamemc["MostCommonMCTrack"] = nsegsamemc["MCTrack"]

del nsegsamemc["MCEvent"]
del nsegsamemc["MCTrack"]

lastplateMCTrack["MostCommonMCEvent"] = lastplateMCTrack["MCEvent"]
lastplateMCTrack["MostCommonMCTrack"] = lastplateMCTrack["MCTrack"]

del lastplateMCTrack["MCEvent"]
del lastplateMCTrack["MCTrack"]

trackdf = trackdf.reset_index().merge(nsegtrue, how = "left", on = ["MostCommonMCEvent","MostCommonMCTrack"])
trackdf = trackdf.reset_index().merge(nsegsamemc, how = "left", on = ["TrackID","MostCommonMCEvent","MostCommonMCTrack"])

#renaming labels
trackdf["PID"] = trackdf["PID_x"]
trackdf["ID"] = trackdf["ID_x"]
trackdf["nsegsamemc"] = trackdf["PID_y"]
trackdf["nsegtrue"] = trackdf["ID_y"]

del trackdf["ID_y"]
del trackdf["ID_x"]
del trackdf["PID_y"]
del trackdf["PID_x"]

trackdf = trackdf.reset_index().merge(lastplateMCTrack, how = "left", on = ["MostCommonMCEvent","MostCommonMCTrack"])
trackdf["PID"] = trackdf["PID_x"]
trackdf["lastplateMCTrack"] = trackdf["PID_y"]

print("end preparation of dataframe, making plots")

#computing efficiency
trackdf["efficiency"] = trackdf["nsegsamemc"]/trackdf["nsegtrue"]

#how many reco tracks for each true track?
nsplit = trackdf.groupby(["MostCommonMCEvent","MostCommonMCTrack"]).count()["ID"]

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
#grouping by MCEvent for clarity
trackdf = trackdf.sort_values("MCEvent")

#printing trackdf
print (trackdf[["TrackID","MCEvent","MCTrack","nseg","nsegsamemc","nsegtrue","efficiency"]])


print(trackdf[["TrackID","MCEvent","MCTrack","nsegtrue","efficiency"]])
