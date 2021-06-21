from __future__ import division
import pandas as pd
import ROOT as r
import fedrarootlogon
from scipy import stats
import numpy as np
import rootnumpy_myutils as myrootnp
import sys
'''launch with ipython -i trackinquality.py nsection'''
nsection = sys.argv[1] #1,2, 3, 4,5,6,7

def GetSectionBorders(nsection):
   '''returns first and last plate of this section'''
   setfile = r.TFile.Open("b000001.{}.0.0.set.root".format(nsection))
   footset = setfile.Get("set")
   nplates = footset.eIDS.GetSize()
   lastplateset = footset.GetID(0).ePlate
   firstplateset = footset.GetID(nplates-1).ePlate
   print("Opened set of section {}, first plate is {}, last plate is {}".format(nsection,firstplateset, lastplateset))
   return firstplateset, lastplateset

firstplateset, lastplateset = GetSectionBorders(nsection)

print("Processing Data Frame for section {}".format(nsection))
#df = pd.read_csv("GSI1_S{}_standard.csv".format(nsection))
df = pd.read_csv("GSI1_S{}.csv".format(nsection))
df["Theta"] = np.arctan(np.sqrt(df["TX"] * df["TX"] + df["TY"] * df["TY"]))
df["Plate"] = lastplateset - df["PID"] #getting plate by pID (we assume no plates are skipped)

simdf = df.query("MCTrack>=0 or TrackID>=0") #segments tracked or from the simulation

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
mostfrequentevent = trackdf.query("MCEvent>=0").groupby(['TrackID'])['MCEvent'].agg(lambda x:x.value_counts().index[0])
mostfrequenttrack = trackdf.query("MCEvent>=0").groupby(['TrackID'])['MCTrack'].agg(lambda x:x.value_counts().index[0])
endmostfrequenttrack = trackdf.query("MCEvent>=0").groupby(['TrackID'])['LastPlate'].agg(lambda x:x.value_counts().index[0])

trackdf = trackdf.groupby("TrackID").last()

#adding computed information to our dataset
trackdf["nseg"] = nseg
trackdf["npl"] = npl
trackdf["fedraeff"] = nseg/npl

print("replacing MCEvent,MCTrack and LastPlate values in trackdf with information from most common MCEvent and MCTrack")
del trackdf["LastPlate"]
trackdf["LastPlate"] = endmostfrequenttrack

del trackdf["MCEvent"]
del trackdf["MCTrack"]
trackdf["MCEvent"] = mostfrequentevent
trackdf["MCTrack"] = mostfrequenttrack

trackdf = trackdf[atleastonemc]

maxnseg = int(np.max(nseg))
#group by same MCTrack and MCEvent
simdf = simdf.query("MCEvent>=0")
nsegtrue = simdf.groupby(["MCEvent","MCTrack"]).count()["ID"] #true length of the Monte Carlo track
maxnsegtrue = int(np.max(nsegtrue))

#let us consider now signal segments which have been tracked
simtrackdf = simdf.query("TrackID>=0")
#grouping by MCEvent and MCTrack, we can get where the track ends (even counting splitting in different reco tracks)
simtrackdf = simtrackdf.sort_values("Plate")
lastplateMCTrack = simtrackdf.groupby(["MCEvent","MCTrack"]).last()["Plate"]

#merging the MCTrue information in trackID, according to most common MCEvent and MCTrack
nsegtrue = nsegtrue.reset_index()
nsegsamemc = nsegsamemc.reset_index()
lastplateMCTrack = lastplateMCTrack.reset_index()

trackdf = trackdf.reset_index().merge(nsegtrue, how = "left", on = ["MCEvent","MCTrack"])
trackdf = trackdf.reset_index().merge(nsegsamemc, how = "left", on = ["TrackID","MCEvent","MCTrack"])

#renaming labels
trackdf["PID"] = trackdf["PID_x"]
trackdf["ID"] = trackdf["ID_x"]
trackdf["nsegsamemc"] = trackdf["PID_y"]
trackdf["nsegtrue"] = trackdf["ID_y"]

del trackdf["ID_y"]
del trackdf["ID_x"]
del trackdf["PID_y"]
del trackdf["PID_x"]

trackdf = trackdf.reset_index().merge(lastplateMCTrack, how = "left", on = ["MCEvent","MCTrack"])
trackdf["Plate"] = trackdf["Plate_x"]
trackdf["lastplateMCTrack"] = trackdf["Plate_y"]

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

#angle of tracked and not tracked Monte Carlo segments
hangletracked = r.TH1D("hangletracked","Theta of tracked MC segments",10,0,1)
hanglemissed = r.TH1D("hanglemissed","Theta of missed MC segments",10,0,1)
cthetatracking = r.TCanvas()

htxtracked = r.TH1D("htxtracked","TX of tracked MC segments",20,-1,1)
htxmissed = r.TH1D("htxmissed","TX of missed MC segments",20,-1,1)
ctxtracking = r.TCanvas()

htytracked = r.TH1D("htytracked","TY of tracked MC segments",20,-1,1)
htymissed = r.TH1D("htymissed","TY of missed MC segments",20,-1,1)
ctytracking = r.TCanvas()

def comparetrackedsegments(variable,htracked,hmissed,canvas):
 '''Tracked and not tracked segments comparison'''

 myrootnp.fillhist1D(htracked, simdf.query("TrackID>=0")[variable])
 myrootnp.fillhist1D(hmissed, simdf[simdf.isnull()["TrackID"]][variable]) #not tracked segments have TrackID NaN

 htracked.Scale(1./htracked.Integral())
 hmissed.Scale(1./hmissed.Integral())
 hmissed.SetLineColor(r.kRed)

 canvas.cd()
 htracked.Draw("histo")
 hmissed.Draw("SAMES and histo")
 canvas.BuildLegend()

comparetrackedsegments("Theta",hangletracked,hanglemissed,cthetatracking)
comparetrackedsegments("TX",htxtracked,htxmissed,ctxtracking)
comparetrackedsegments("TY",htytracked,htymissed,ctytracking)


#csplitlength = r.TCanvas()
#myrootnp.fillhist2D(hsplitlength,trackdf["nsegtrue"],nsplit)
#hsplitlength.Draw("COLZ")
#grouping by MCEvent for clarity
trackdf = trackdf.sort_values("MCEvent")

#printing trackdf and some information
print(trackdf[["TrackID","MCEvent","MCTrack","nsegtrue","efficiency"]])


print("len(trackdf) = {}".format(len(trackdf)))
print("len(trackdf.query(nsegsamemc<nseg)) = {}".format(len(trackdf.query("nsegsamemc<nseg"))))

#how many particles actually end before the end of the section?
endingtracksdf = trackdf.query("LastPlate<{}".format(lastplateset))
print("{} particles stop before the end of the section".format(len(endingtracksdf)))
#where does track reconstruction actually stop, in this subsample
endingtracksdf["missinglastplates"] = endingtracksdf.eval("LastPlate - lastplateMCTrack")

hlastplates = r.TH1D("hlastplates","Reconstructed tracks stop N plates before the actual end of the particle;NPlates",6,0,6)
clastplates = r.TCanvas()
myrootnp.fillhist1D(hlastplates, endingtracksdf["missinglastplates"])
hlastplates.Draw()

#selecteddf = df.query("MCEvent==71 and MCTrack==3")

seglist = r.TObjArray()
trackstodraw = r.TObjArray()

ds=r.EdbDisplay("Display tracked and not tracked segments",-60000.,60000.,-50000.,50000.,-4000.,80000.)
ds.SetDrawTracks(4)
#opening track file, we need it if we want to access all linked tracks information, such as fitted segments
trackfile = r.TFile.Open("b000001.{}.0.0.trk.root".format(nsection))
tracktree = trackfile.Get("tracks")

def drawsegments(selection):
 '''draw segments from dataframe according to selection'''
 #clear arrays for this event
 seglist.Clear()
 trackstodraw.Clear()
 #applyselection for this event and loop over it
 selecteddf = df.query(selection)
 selecteddf = selecteddf.fillna(-1); #not tracked are labelled as -1 
 for index, row in selecteddf.iterrows():
  myseg = r.EdbSegP(int(row["ID"]),float(row["x"]),float(row["y"]),float(row["TX"]),float(row["TY"]))
  trackID = int(row["TrackID"])
  
  if (trackID>=0):
   tracktree.GetEntry(trackID)
   origtrack = tracktree.t
   segments = tracktree.s
   fittedsegments = tracktree.sf

   copiedtrack = r.EdbTrackP()
   copiedtrack.Copy(origtrack)  
  
   #adding a copy of the track, adding segments
   for (segment, fittedsegment) in zip(segments,fittedsegments):
    segment.SetDZ(300)
    fittedsegment.SetDZ(300)
    copiedtrack.AddSegment(r.EdbSegP(segment))
    copiedtrack.AddSegmentF(r.EdbSegP(fittedsegment))
   copiedtrack.SetSegmentsTrack(copiedtrack.ID()) #track segments association
   copiedtrack.SetCounters()
   trackstodraw.Add(copiedtrack);
  else:
   #other segment information
   myseg.SetDZ(300)
   myseg.SetZ(float(row["z"]))
   myseg.SetPID(int(row["PID"]))
   #customized segment line (required to see the segments in new ROOT versions?)
   mysegline = ds.SegLine(myseg)
   mysegline.SetLineColor(r.kBlue)
   mysegline.SetLineWidth(6)
   seglist.Add(mysegline)

 ds.SetArrSegG(seglist) #not setarrsegp, since I have customized the lines
 ds.SetArrTr(trackstodraw)
 ds.Draw()

 return seglist
