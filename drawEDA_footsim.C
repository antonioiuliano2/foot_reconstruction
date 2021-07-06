//EDA customized setup to inspect how Monte Carlo Tracks are reconstructed
void drawEDA_footsim(){
 //which tracks do you want to draw?
 TCut mycut("1");

//getting tracks set from file
 int nsection = 3;
 cout<<"Launching display for section "<<nsection<<endl;
 EdbEDA* gEDA = new EdbEDA(Form("b000001.%i.0.0.trk.root",nsection),100,mycut, kFALSE);
 
 gEDA->GetTrackSet("TS")->SetColorMode(kCOLOR_BY_MCID); //added by me, works like ID, but with MCTrack; 

 gEDA->GetTrackSet("TS")->SetDraw(kTRUE);
 //gEDA->Redraw();
 gEDA->Run();
 gEDA->GetTrackSet("TS")->SetNsegCut(2); //default is at least 3 ses
 gEDA->GetTrackSet("TS")->DoSelection();
 //gEDA->Redraw(); //i do not redraw all tracks with at least 2 seg (memory consumption), I just need them for single events

 cout<<"Display ready, please  launch drawevent(nevent) to draw tracks with at least one segment from that event"<<endl;
 cout<<"Segments from different Monte Carlo tracks are coloured differently"<<endl;
} 

void drawevent(int inputMCEvent){
 EdbEDATrackSet *mcset = gEDA->GetTrackSet("SB");
 mcset->Clear();
 bool atleast1seg = false; //track is stored if contains at least one segment from this event
 const int ntracks = gEDA->GetTrackSet("TS")->N();
 for (int itrk = 0; itrk < ntracks; itrk++){
     EdbTrackP *mytrack = gEDA->GetTrackSet("TS")->GetTrack(itrk);
     atleast1seg = false;
     for (int iseg = 0; iseg < mytrack->N(); iseg++){
         EdbSegP *seg = mytrack->GetSegment(iseg);
         if(seg->MCEvt()==inputMCEvent) atleast1seg=true;
     }
     if (atleast1seg) mcset->AddTrack(mytrack);
 }

 gEDA->GetTrackSet("TS")->SetDraw(kFALSE);
 mcset->SetTrackAttribute(-1);
 mcset->SetColorMode(kCOLOR_BY_MCID);
 mcset->SetDraw(kTRUE);
 gEDA->Redraw();

}
