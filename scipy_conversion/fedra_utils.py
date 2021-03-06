import numpy as np
import pandas as pd
import fedrarootlogon
import ROOT as r

def builddataframe(brick, path = "..", cutstring = "1", major = 0, minor = 0, newzprojection = None, footsim = False):
 """build pandas dataframe starting from couples and scanset 
    brick = Number of brick as in b0000*
    path = input path to the folder containing theb b0000* folder
    cutsring = eventual selection to couples
    newzprojection = list of projection to a new z reference system
 """
 nplate =0

 print("Reading ScanSet at path ",path)

 #reading scanset
 sproc = r.EdbScanProc()
 sproc.eProcDirClient=path
 id = r.EdbID(brick,nplate,major,minor)
 ss = sproc.ReadScanSet(id)
 ss.Brick().SetID(brick)
 
 #preparing patterns
 npl = ss.eIDS.GetEntries()

 cut = r.TCut(cutstring)

 #intial empty arrays
 IDall = np.zeros(0,dtype=int)
 PIDall = np.zeros(0,dtype=int)

 xall = np.zeros(0,dtype=np.float32)
 yall = np.zeros(0,dtype=np.float32)
 zall = np.zeros(0,dtype=np.float32)
 TXall = np.zeros(0,dtype=np.float32)
 TYall = np.zeros(0,dtype=np.float32)

 MCEvtall = np.zeros(0,dtype=int)
 MCTrackall = np.zeros(0,dtype=int)
 MCMotherIDall = np.zeros(0,dtype=int)
 Pall = np.zeros(0,dtype=np.float32)
 Flagall = np.zeros(0,dtype=int)
 Weightall = np.zeros(0,dtype=np.float32)

 FirstPlateall = np.zeros(0,dtype=int) #added by Giuliana: where does the track end?
 LastPlateall = np.zeros(0,dtype=int)

 print ("Cut on couples ")
 cut.Print()

 print("Try to open folders at path ",path+"/b00000"+str(brick))
 for i in range(npl):
  idplate = ss.GetID(i)
      
  nplate = idplate.ePlate
  plate = ss.GetPlate(idplate.ePlate)
  #read pattern information
  p = r.EdbPattern()

  ect = r.EdbCouplesTree()
  if (nplate) <10:
   ect.InitCouplesTree("couples",path+"/b00000"+str(brick)+"/p00{}/{}.{}.{}.{}.cp.root".format(nplate,brick,nplate,major,minor),"READ")
  elif (nplate)<100:
   ect.InitCouplesTree("couples",path+"/b00000"+str(brick)+"/p0{}/{}.{}.{}.{}.cp.root".format(nplate,brick,nplate,major,minor),"READ")
  else:
   ect.InitCouplesTree("couples",path+"/b00000"+str(brick)+"/p{}/{}.{}.{}.{}.cp.root".format(nplate,brick,nplate,major,minor),"READ")
  #addingcut
  ect.eCut = cut 
  cutlist = ect.InitCutList()
  
  nsegcut = cutlist.GetN()
  nseg = ect.eTree.GetEntries()

  IDarray_plate = np.zeros(nsegcut,dtype=int)
  PIDarray_plate = np.zeros(nsegcut,dtype=int)

  xarray_plate = np.zeros(nsegcut,dtype=np.float32)
  yarray_plate = np.zeros(nsegcut,dtype=np.float32)
  zarray_plate = np.zeros(nsegcut,dtype=np.float32)
  TXarray_plate = np.zeros(nsegcut,dtype=np.float32)
  TYarray_plate = np.zeros(nsegcut,dtype=np.float32)
   
  MCEvtarray_plate = np.zeros(nsegcut,dtype=int)
  MCTrackarray_plate = np.zeros(nsegcut,dtype=int)
  MCMotherIDarray_plate = np.zeros(nsegcut,dtype=int)
  Parray_plate = np.zeros(nsegcut,dtype=np.float32)
  Flagarray_plate = np.zeros(nsegcut,dtype=int)
  Weightarray_plate = np.zeros(nsegcut,dtype=np.float32)

  FirstPlatearray_plate = np.zeros(nsegcut,dtype=int)
  LastPlatearray_plate = np.zeros(nsegcut,dtype=int)

  print ("loop on {} segments over  {} for plate {}".format(nsegcut, nseg,nplate))
  for ientry in range(nsegcut):
   iseg = cutlist.GetEntry(ientry)
   ect.GetEntry(iseg)
 
   seg=ect.eS
   #//setting z and affine transformation
   seg.SetZ(plate.Z())
   seg.SetPID(i)
   seg.Transform(plate.GetAffineXY())

   if(newzprojection is not None):
    seg.PropagateTo(newzprojection[i])

   IDarray_plate[ientry] = seg.ID()
   PIDarray_plate[ientry] = seg.PID()
   
   xarray_plate[ientry] = seg.X()
   yarray_plate[ientry] = seg.Y()
   zarray_plate[ientry] = seg.Z()
   TXarray_plate[ientry] = seg.TX()
   TYarray_plate[ientry] = seg.TY()

   MCEvtarray_plate[ientry] = seg.MCEvt()
   MCTrackarray_plate[ientry] = seg.MCTrack()
   MCMotherIDarray_plate[ientry] = seg.Aid(1)
   Parray_plate[ientry] = seg.P()     
   Weightarray_plate[ientry] = seg.W()     
   if footsim: #different place where pdgcode is stored
    Flagarray_plate[ientry] = seg.Aid(0)
   else:
    Flagarray_plate[ientry] = seg.Flag()
   
   FirstPlatearray_plate[ientry] = seg.Vid(0)
   LastPlatearray_plate[ientry] = seg.Vid(1)  

  #end of loop, storing them in global arrays
  IDall = np.concatenate((IDall,IDarray_plate),axis=0)
  PIDall = np.concatenate((PIDall,PIDarray_plate),axis=0)

  xall = np.concatenate((xall,xarray_plate),axis=0)
  yall = np.concatenate((yall,yarray_plate),axis=0)
  zall = np.concatenate((zall,zarray_plate),axis=0)
  TXall = np.concatenate((TXall,TXarray_plate),axis=0)
  TYall = np.concatenate((TYall,TYarray_plate),axis=0)
  MCEvtall = np.concatenate((MCEvtall,MCEvtarray_plate),axis=0)
  MCTrackall = np.concatenate((MCTrackall,MCTrackarray_plate),axis=0)
  MCMotherIDall = np.concatenate((MCMotherIDall,MCMotherIDarray_plate),axis=0)
  Pall = np.concatenate((Pall,Parray_plate),axis=0)
  Weightall = np.concatenate((Weightall,Weightarray_plate),axis=0)
  Flagall = np.concatenate((Flagall,Flagarray_plate),axis=0)

  FirstPlateall = np.concatenate((FirstPlateall, FirstPlatearray_plate),axis=0)
  LastPlateall = np.concatenate((LastPlateall, LastPlatearray_plate),axis=0)

 data = {'ID':IDall,'PID':PIDall,'x':xall,'y':yall,'z':zall,'TX':TXall,'TY':TYall,'MCEvent':MCEvtall,'MCTrack':MCTrackall,'MCMotherID':MCMotherIDall,'P':Pall,'Flag':Flagall, 
'Weight':Weightall,  "FirstPlate":FirstPlateall,"LastPlate":LastPlateall}
 df = pd.DataFrame(data, columns = ['ID','PID','x','y','z','TX','TY','MCEvent','MCTrack','MCMotherID','P','Flag','Weight',"FirstPlate","LastPlate"] )

 return df


def addtrueMCinfo(df,simfile, ship_charm):
 '''getting additional true MC info from source file, If ship_charm is true, TM is taken into account to spread XY position'''
 import pandas as pd
 import numpy as np
 import ROOT as r
 
 simtree = simfile.Get("cbmsim")

 #position differences from FairShip2Fedra: initialized to 0
 xoffset = 0.
 yoffset = 0.
 zoffset = 0.
 #computing zoffset: in our couples, most downstream plate has always z=0
 simtree.GetEntry(0)
 emulsionhits = simtree.BoxPoint
 ihit = 0
 while (zoffset >= 0.):
  hit = emulsionhits[ihit]  
  if (hit.GetDetectorID()==29):
   zoffset = 0. - hit.GetZ()
  ihit = ihit + 1

 #virtual TM parameters (only ship-charm simulations)
 spilldy = 1
 targetmoverspeed = 2.6

 print("ZOffset between FairShip and FEDRA",zoffset) 

 df = df.sort_values(["MCEvent","MCTrack","PID"],ascending=[True,True,False]) #sorting by MCEventID allows to access each event only once
 df.reset_index()

 currentevent=-1

 nsegments = len(df)
 isegment = 0
 print("dataframe prepared, starting loop over {} segments".format(nsegments))

 #preparing arrays with new columns
 arr_MotherID = np.zeros(nsegments, dtype=int)
 arr_ProcID = np.zeros(nsegments, dtype=int)

 arr_startX = np.zeros(nsegments,dtype=np.float32)
 arr_startY = np.zeros(nsegments,dtype=np.float32)
 arr_startZ = np.zeros(nsegments,dtype=np.float32)
 arr_startT = np.zeros(nsegments,dtype=np.float32)
 
 arr_startPx = np.zeros(nsegments,dtype=np.float32)
 arr_startPy = np.zeros(nsegments,dtype=np.float32)
 arr_startPz = np.zeros(nsegments,dtype=np.float32)

 for (MCEvent, MCTrack) in zip(df['MCEvent'], df['MCTrack']):

  if (MCEvent != currentevent):
   currentevent = MCEvent
   simtree.GetEntry(currentevent)
   eventtracks = simtree.MCTrack
   if(currentevent%10000 == 0):
    print("Arrived at event ",currentevent)
   #getting virtual timestamp and replicating Target Mover
   if (ship_charm):
    eventheader = simtree.ShipEventHeader
    virtualtimestamp = eventheader.GetEventTime()
    nspill = int(virtualtimestamp/100)
    pottime = virtualtimestamp - nspill * 100
    #computing xy offset for this event
    xoffset = -12.5/2. + pottime * targetmoverspeed
    yoffset = - 9.9/2. + 0.5 + nspill * spilldy 
    #print(virtualtimestamp, pottime, nspill, xoffset, yoffset)
   

  #adding values
  mytrack = eventtracks[MCTrack]
  arr_MotherID[isegment] = mytrack.GetMotherId()
  arr_ProcID[isegment] = mytrack.GetProcID()

  arr_startX[isegment] = (mytrack.GetStartX() + xoffset) * 1e+4 + 62500 #we need also to convert cm to micron
  arr_startY[isegment] = (mytrack.GetStartY() + yoffset) * 1e+4 + 49500
  arr_startZ[isegment] = (mytrack.GetStartZ() + zoffset) * 1e+4
  arr_startT[isegment] = mytrack.GetStartT()
  
  arr_startPx[isegment] = mytrack.GetPx()
  arr_startPy[isegment] = mytrack.GetPy()
  arr_startPz[isegment] = mytrack.GetPz()
  
  isegment = isegment + 1
 
 #adding the new columns to the dataframe
 df["MotherID"] = arr_MotherID 
 df["ProcID"] = arr_ProcID

 df["StartX"] = arr_startX
 df["StartY"] = arr_startY
 df["StartZ"] = arr_startZ
 df["StartTime"] = arr_startT

 df["startPx"] = arr_startPx
 df["startPy"] = arr_startPy
 df["startPz"] = arr_startPz
 
 return df

def addvertextrackindexes(df, vertexfilename):
 '''adding track and vertex index to dataframe'''
 vertexfile = r.TFile.Open(vertexfilename)
 vrec = vertexfile.Get("EdbVertexRec")
 vtxcontainer = vrec.eVTX
 
 nvertices = vtxcontainer.GetEntries()

 #initial empty arrays, to be filled with segments from all tracks
 IDall = np.zeros(0,dtype=int)
 PIDall = np.zeros(0,dtype=int)
 TrackIDall = np.zeros(0,dtype=int)
 VertexIDall = np.zeros(0,dtype=int)
 #ID and PID do NOT identify unique segment in FOOT simulations, since signal and backgrounds are merged 
 MCEventall = np.zeros(0,dtype=int)
 MCTrackall = np.zeros(0,dtype=int)
 print("start loop on {} vertices".format(vtxcontainer.GetEntries()))
 for ivertex, vertex in enumerate (vtxcontainer):
  ntracks = vertex.N()
  for itrack in range(ntracks):
   #getting track with index itrack
   track = vertex.GetTrack(itrack)
   nseg = track.N()

   #initial arrays, length given by number of segments of this tracks
   IDarr = np.zeros(nseg,dtype=int)
   PIDarr = np.zeros(nseg,dtype=int)
   TrackIDarr = np.zeros(nseg,dtype=int)
   VertexIDarr = np.zeros(nseg,dtype=int)
   #add true MC info
   MCEventarr = np.zeros(nseg,dtype=int)
   MCTrackarr = np.zeros(nseg,dtype=int)

   #start loop on segments
   for iseg in range(nseg):
    seg = track.GetSegment(iseg)
    IDarr[iseg] = seg.ID()
    PIDarr[iseg] = seg.PID()
    TrackIDarr[iseg] = track.Track() #seg.Track() may be the same for different tracks (also when selection is applied, use track.Track() for having the original track)
    VertexIDarr[iseg] = ivertex #not vertexID (they are not the same in vertices_improved
    #add true MC info
    MCEventarr[iseg] = seg.MCEvt()
    MCTrackarr[iseg] = seg.MCTrack()
    
    if(seg.MCEvt()==1 and seg.MCTrack()==3):
     print("Found segment: ",seg.MCEvt(),seg.MCTrack(),seg.ID(), seg.PID(), track.Track(),ivertex)
   
   #concatenate, adding segments for this track to the global arrays
   IDall = np.concatenate((IDall,IDarr),axis=0)
   PIDall = np.concatenate((PIDall,PIDarr),axis=0)
   TrackIDall = np.concatenate((TrackIDall,TrackIDarr),axis=0)
   VertexIDall = np.concatenate((VertexIDall,VertexIDarr),axis=0)
   #true MC info
   MCEventall = np.concatenate((MCEventall,MCEventarr),axis=0)
   MCTrackall = np.concatenate((MCTrackall,MCTrackarr),axis=0)

 #end of loop, preparing dataframe
 labels = ["ID","PID","TrackID","VertexID","MCEvent","MCTrack"]
 dftracks = pd.DataFrame({"ID":IDall,"PID":PIDall,"TrackID":TrackIDall,"VertexID":VertexIDall,"MCEvent":MCEventall,"MCTrack":MCTrackall},columns = labels) 
 print("Track dataframe ready: merging it with all couples dataframe: not tracked segments will be labelled as NA") 
 #removing duplicates (same track can belong to more than one vertex, so the same segment can be filled multiple times)
 dftracks = dftracks.groupby(["PID","ID","MCEvent","MCTrack"]).first()
 dftracks = dftracks.reset_index()
 #Now I need to merge them, however I want to keep all the segments, not only the ones which have been tracked. Luckily, there are many ways to do a merge (default is inner)
 dfwithtracks = df.merge(dftracks,how = 'left', on=["PID","ID","MCEvent","MCTrack"])
 return dfwithtracks
 

def addtrackindex(df, trackfilename):
 ''' adding track index to dataframe, if tracking was performed'''
 trackfile = r.TFile.Open(trackfilename)
 tracktree = trackfile.Get("tracks")
 
 ntracks = tracktree.GetEntries()
 #initial empty arrays, to be filled with segments from all tracks
 IDall = np.zeros(0,dtype=int)
 PIDall = np.zeros(0,dtype=int)
 TrackIDall = np.zeros(0,dtype=int)
 #ID and PID do NOT identify unique segment in FOOT simulations, since signal and backgrounds are merged 
 MCEventall = np.zeros(0,dtype=int)
 MCTrackall = np.zeros(0,dtype=int)
 print("start loop on {} tracks".format(tracktree.GetEntries()))
 for itrack,track in enumerate(tracktree):
  nseg = track.nseg
  segments = track.s

  #initial arrays, length given by number of segments of this tracks
  IDarr = np.zeros(nseg,dtype=int)
  PIDarr = np.zeros(nseg,dtype=int)
  TrackIDarr = np.zeros(nseg,dtype=int)
  #add true MC info
  MCEventarr = np.zeros(nseg,dtype=int)
  MCTrackarr = np.zeros(nseg,dtype=int)

  #start loop on segments
  for iseg, seg in enumerate(segments):
   IDarr[iseg] = seg.ID()
   PIDarr[iseg] = seg.PID()
   TrackIDarr[iseg] = itrack #seg.Track() may be the same for different tracks
   #add true MC info
   MCEventarr[iseg] = seg.MCEvt()
   MCTrackarr[iseg] = seg.MCTrack()
   
  #concatenate, adding segments for this track to the global arrays
  IDall = np.concatenate((IDall,IDarr),axis=0)
  PIDall = np.concatenate((PIDall,PIDarr),axis=0)
  TrackIDall = np.concatenate((TrackIDall,TrackIDarr),axis=0)
  #true MC info
  MCEventall = np.concatenate((MCEventall,MCEventarr),axis=0)
  MCTrackall = np.concatenate((MCTrackall,MCTrackarr),axis=0)

 #end of loop, preparing dataframe
 labels = ["ID","PID","TrackID","MCEvent","MCTrack"]
 dftracks = pd.DataFrame({"ID":IDall,"PID":PIDall,"TrackID":TrackIDall,"MCEvent":MCEventall,"MCTrack":MCTrackall},columns = labels) 
 print("Track dataframe ready: merging it with all couples dataframe: not tracked segments will be labelled as NA") 
 #Now I need to merge them, however I want to keep all the segments, not only the ones which have been tracked. Luckily, there are many ways to do a merge (default is inner)
 dfwithtracks = df.merge(dftracks,how = 'left', on=["PID","ID","MCEvent","MCTrack"])
 return dfwithtracks

def retrievetrackinfo(df, trackfilename):
  '''from trackindex, make a dataframe with fitted position and angles of the tracks'''
  trackfile = r.TFile.Open(trackfilename)
  tracktree = trackfile.Get("tracks")

  dftracked = df.query("TrackID>=0") #we access the subset with tracks
  #we keep only the first segment of eacht track
  dftracked = dftracked.sort_values("PID",ascending=False)
  dftracked = dftracked.groupby("TrackID").first()
  dftracked = dftracked.reset_index()
  #how many tracks we have? What are their TrackIDs?
  ntracks = len(dftracked)
  trackID = dftracked["TrackID"].to_numpy(dtype=int)
  #preparing arrays to store informatino
  Xarr = np.zeros(ntracks,dtype=np.float32)
  Yarr = np.zeros(ntracks,dtype=np.float32)
  TXarr = np.zeros(ntracks,dtype=np.float32)
  TYarr = np.zeros(ntracks,dtype=np.float32)
  #loop over tracks
  for trid in trackID:
    #retrieving trackinfo
    tracktree.GetEntry(trid)
    track = tracktree.t
    
    Xarr[trid]=track.X()
    Yarr[trid]=track.Y()
    TXarr[trid]=track.TX()
    TYarr[trid]=track.TY()
  #end of loop, preparing outputdataframe and returning it
  trackdf = pd.DataFrame({"TrackID":trackID,"X":Xarr,"Y":Yarr,"TX":TXarr,"TY":TYarr},columns = ["TrackID","X","Y","TX","TY"])
  return trackdf 

  

def addvertexindex(df,vertexfilename):
  '''adding vertex index to dataframe. Requires Track Index'''
  vertexfile = r.TFile.Open(vertexfilename)
  vertextree = vertexfile.Get("vtx")
  
  nvertices = vertextree.GetEntries()
  #initial empty arrays, to be filled with segments from all tracks
  #TrackIDall = np.zeros(0,dtype=int)
  VertexS = {}
  VertexE = {}

  #start loop on vertices
  for vtx in vertextree:
  
   vID = vtx.vID
   ntracks = vtx.n
   trackIDs = vtx.TrackID
   trackedges = vtx.incoming
   #start loop on tracks
   for itrack in range(ntracks):
    trackID = trackIDs[itrack]
    #initial values, placeholder for vertex indexes
    if trackID not in VertexS:
     VertexS[trackID] = -1
    if trackID not in VertexE:
     VertexE[trackID] = -1
    #filling tracks
    if (trackedges[itrack] == 0):
     VertexE[trackID] = vID
    else:
     VertexS[trackID] = vID

  #getting array with all tracks
  labels = ["FEDRATrackID","VertexS","VertexE"]
  TrackIDarr = list(VertexS.keys())
  VertexSarr = list(VertexS.values())
  VertexEarr = list(VertexE.values())
  dfvertices = pd.DataFrame({"FEDRATrackID":TrackIDarr, "VertexS":VertexSarr, "VertexE":VertexEarr},columns = labels, dtype = int)
  
  #dfwithvertices = df.merge(dfvertices,how='left', on = ["TrackID"])

  return dfvertices  
