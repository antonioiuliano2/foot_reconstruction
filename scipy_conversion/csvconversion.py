import fedra_utils as feutils
import sys
import os
import pandas as pd
import ROOT as r
'''
#python csvconversion.py nsection
nsection = sys.argv[1]
#please link set.root to the right file for your section
answer = input("ATTENTION: file b000001.0.0.0.set.root will be DELETED. Please confirm with 1 that this is the link file you wish to delete: ")
if (answer is not 1):
 print("exiting")
 exit()

os.system("rm b000001.0.0.0.set.root")
os.system("ln -s b000001.{}.0.0.set.root b000001.0.0.0.set.root".format(nsection))
'''
def applyconversion(nbrick,section):
 '''convert couples ROOT files into a csv'''

 df = feutils.builddataframe(nbrick,cutstring="s.eMCEvt>=0", footsim=True)
 #df = feutils.addtrackindex(df,"b00000{}.{}.0.0.trk.root".format(nbrick,nsection))

 return df 

#the two steps can now be done together, without an intermediate file
NGSI = 1
df = applyconversion(NGSI,0)

#df = df.drop(columns = ["P","Flag"])
#simfile = r.TFile.Open(sys.argv[2])
#df = desy19.addtrueMCinfo(df,simfile, True)
#df.to_csv('GSI1_S{}.csv'.format(nsection),index=False)
df.to_csv('GSI{}_simonly.csv'.format(NGSI),index=False)
