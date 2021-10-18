'''make the trackid a separate file, so it will be easier to adapt this command to different versions by Giuliana :) (Antonio, 1st July)'''
import fedra_utils as feutils
import pandas as pd

df = pd.read_csv("GSI1_simonly.csv")

#df["TrackID"].to_csv("GSI1_tracks.csv")

df = feutils.addvertextrackindexes(df,"vertices_improved.root")
#removing trackid from vertexid
df["Vertex_TrackID"] = df["TrackID"]
df = df.drop("TrackID",axis = 1)

df = feutils.addtrackindex(df,"b000001.0.1.7.trk.root")

#df[["TrackID"]].to_csv("GSI1_tracks_vertices_simulationonly_alltracks.csv") #not only ones from vertices
df[["TrackID","VertexID"]].to_csv("GSI1_tracks_vertices_simulationonly_tracksandvertices.csv") #vertices_improved and tracks
