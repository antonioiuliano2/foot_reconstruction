//----------------------------------------------------------------------------
//
//  Usage: $ root -l foot_vertexing.C
/*
 scp /Users/Giuliana/Desktop/Uni/FOOT/TEST_RECO_MIC/foot_vertexing.C scanner@nusrv9.na.infn.it:/home/scanner/foot/RECO_MC_TEST_giul/b000001/
 
 scp /Users/Giuliana/Desktop/Uni/FOOT/TEST_RECO_MIC/foot_vertexing.C scanner@nusrv9.na.infn.it:/home/scanner/foot/2019_GSI/GSI2/b000002
 
 
 */
//----------------------------------------------------------------------------

#define gEDBDEBUGLEVEL 1
#define LASTLAYER 120
#define DIRECTION 1 // 9 = indietro, 1 = avanti
#define BRICKID 1

//ATTENZIONE BISOGNA MODIFICARE DZmax a secnda del brick

/*
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TH1.h"
#include "TCanvas.h"
#include <stdio.h>
#include <TROOT.h>
#include "TRandom3.h"
#include "TVector.h"
#include <vector>
#include <algorithm>
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbCouplesTree.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbDataSet.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbScanSet.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbSegP.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbPVRec.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbPattern.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbTrackFitter.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbVertex.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbDisplay.h"
*/

EdbDataProc  *dproc=0;
EdbPVRec     *gAli=0;
EdbVertexRec *gEVR=0;
EdbDisplay   *ds=0;
float Z_LAYER[LASTLAYER+1]={0};

TH1D *hip = new TH1D("hip","Impact parameters",500,0,5000);
//EdbScanSet *set = file->Get("set");


void trvol( const char *def, const char *rcut = "nseg>1" );
void init( const char *def, int iopt,  const char *rcut="1" );
void set_segments_dz(float dz=300.);
void do_propagation();
void do_vertex();
void td();
void sd();
void vd( int trmin, float amin);

#include <vector>
namespace VERTEX_PAR
{
    float DZmax = 4200; //4200 (GSI1);//7200. (GSI2);
    float ProbMinV   = 0.001;  // minimum acceptable probability for chi2-distance between tracks
    //float ProbMinV   = 0.01;
    float ImpMax     = 60.; //70.;    // maximal acceptable impact parameter [microns] (for preliminary check)
    bool  UseMom     = false;  // use or not track momentum for vertex calculations
    bool  UseSegPar  = false;  // use only the nearest measured segments for vertex fit (as Neuchatel)
    int   QualityMode= 1;      // vertex quality estimation method (0:=Prob/(sigVX^2+sigVY^2); 1:= inverse average track-vertex distance)
}

//---------------------------------------------------------------------
void foot_vertexing_oneEvent_only(char *dset=0)
{
    //trvol(dset);
    //trvol(dset,"nseg>1 && TMath::Abs(s.eY-50000)<5000 && TMath::Abs(s.eX-50000)<5000");
    // SELECTION WHICH TRACKS TO USE FOR VERTEXING
    //trvol(dset,"nseg>1 && TMath::Abs(s.eY-24000)<62500 && TMath::Abs(s.eX-24000)<50000");
    trvol(dset,"nseg>1 && s.eMCEvt == 364 && s.eY>30000 && s.eY<70000 && s.eX>40000 && s.eX<85000");

    // reconstruct vertexes starting from linked_tracks.root
}

//---------------------------------------------------------------------
void trvol( const char *def, const char *rcut )
{
    // this function read volume tracks and do the vertex reconstruction
    // from linked_tracks.root
    
    init(def, 100 ,rcut);                      // read tracks (option 100)
    gAli->FillCell(30,30,0.009,0.009);
    get_z_pattern();
    do_vertex();
    //vd(4,0.01);   // draw reconstructed vertex
    //td(); //to draw tracks
    //dproc->MakeTracksTree(gAli,"linked_tracks_p.root");
}

//---------------------------------------------------------------------
void init( const char *def, int iopt,  const char *rcut)
{
    if(!def)  dproc = new EdbDataProc();
    else      dproc = new EdbDataProc(def);
    
    dproc->InitVolume(iopt, rcut);
    gAli = dproc->PVR();
    checkpatterns();
    set_segments_dz(300.);

}

void checkpatterns(){
    //check for patterns, if there is one missing add it
    int np = gAli->Npatterns();
    
    TFile *setfile = TFile::Open(Form("b%06d.0.0.0.set.root",BRICKID));
    EdbScanSet *set = (EdbScanSet*) setfile->Get("set");
    
    for(int i=0; i<np; i++) {
        EdbPattern *p = gAli->GetPattern(i);
        if(!p) {
            cout<<"missing pattern "<<i<<" now adding it "<<endl;
            float zmissingPID = set->eB.GetPlate(i)->Z(); //note, set->eB->GetPlate(i) uses the PID, set->GetPlate(i) uses the number of plate, be careful!
            EdbPattern *pat = new EdbPattern( 0., 0.,zmissingPID);
            pat->SetID(i);
            pat->SetScanID(0);
            gAli->AddPatternAt(pat,i);
        }//end if
    } //end for loop
}

//---------------------------------------------------------------------
void set_segments_dz(float dz)
{
    int np = gAli->Npatterns();
    for(int i=0; i<np; i++) {
        EdbPattern *p = gAli->GetPattern(i);
        int ns = p->N();
        for(int j=0; j<ns; j++) p->GetSegment(j)->SetDZ(dz);
    }
}

//---------------------------------------------------------------------
void do_vertex()
{
    using namespace VERTEX_PAR;
    
    
    //gAli->FitTracks(4.,0.139 );
    gEVR = new EdbVertexRec();
    gEVR->eEdbTracks = gAli->eTracks;
    gEVR->eVTX       = gAli->eVTX;
    gEVR->SetPVRec(gAli);
    
    gEVR->eDZmax=DZmax;
    gEVR->eProbMin=ProbMinV;
    gEVR->eImpMax=ImpMax;
    gEVR->eUseMom=UseMom;
    gEVR->eUseSegPar=UseSegPar;
    gEVR->eQualityMode=QualityMode;
    
    //    EdbTrackP * testtrack = ( EdbTrackP *) gEVR->eEdbTracks->At(0);
    //    cout << "testrack: " << testtrack->X() << endl;
    //    testtrack->COV().Print();
    //
    //
    
    printf("%d tracks for vertexing\n",  gEVR->eEdbTracks->GetEntries() );
    int nvtx = gEVR->FindVertex(); //QUI FA IL VERTEXING, LA FUNZIONE SI TROVA IN EdbVertex.cxx
    printf("%d 2-track vertexes was found\n",nvtx);
    
    if(nvtx == 0) return;
    int nadd =  gEVR->ProbVertexN();
    cout << "nadd: " << nadd << endl;
    int nl = gEVR->LinkedVertexes(); //should avoid same track associated to multiple vertices
    printf("%d vertices are linked\n",nl);
    
    //    EdbVertex *testvertex = (EdbVertex *) gEVR->eVTX->At(0);
    //    cout <<" v->vx() : "<<   testvertex->V()->vx() << endl;
    
    //ADD Mod Gg
    /*
    printf("Starting VertexNeighbor\n");
    gEVR->VertexNeighbor(3500,3,400);
    printf("Starting VertexTuning\n");
    int n_vert_tuning = gEVR->VertexTuning(0);
    printf("%d VertexTuning\n",n_vert_tuning);
    //fine Add
     */
    
    
    TFile *fvtx = new TFile("vertices_onlyonevent.root","RECREATE"); //OUTPUT FILE NAME
    TTree *vtx = new TTree("vtx","vtx");
    Float_t vx, vy, vz;
    Int_t vID = 0;
    Float_t maxaperture;
    Float_t probability;
    Int_t n;
    Int_t v_flag;
    Int_t vplate;
    const Int_t maxdim = 1000;
    //big arrays for containers of track variables
    Int_t IDTrack[maxdim]={0};
    Int_t maxgap[maxdim]={0};
    Int_t nseg[maxdim]={0};
    Int_t nseg_S1[maxdim]={0};
    Int_t npl[maxdim]={0};
    Int_t npl_S1[maxdim]={0};
    Int_t nholes[maxdim]={0};
    Int_t nholes_S1[maxdim]={0};
    Int_t plate[maxdim]={0};
    Float_t X[maxdim]={0};
    Float_t Y[maxdim]={0};
    Float_t Z[maxdim]={0};
    Float_t TX[maxdim]={0};
    Float_t TY[maxdim]={0};
    Float_t Theta[maxdim]={0};
    Float_t impactparameter[maxdim]={0};
    Int_t incoming[maxdim]={0}; //coming from the vertex
    Int_t Z_flag[maxdim]={0};
    Int_t MC_Charge_first[maxdim]={0};
    Int_t MC_Charge_last[maxdim]={0};
    Int_t MC_Charge_S2[maxdim]={0};
    Int_t MC_IDPart[maxdim]={0};
    Int_t MC_evID_last[maxdim]={0};
    Int_t MC_evID_first[maxdim]={0};
    Int_t MC_trackID_first[maxdim]={0};
    Int_t MC_trackID_last[maxdim]={0};
    Int_t MC_IDPart_last[maxdim]={0};
    Int_t MC_nseg[maxdim]={0};
    Int_t MC_mother_first[maxdim]={0};
    Int_t MC_mother_last[maxdim]={0};
    Int_t MC_pdgcode_first[maxdim]={0};
    Int_t MC_pdgcode_last[maxdim]={0};

    Int_t MC_firstplate[maxdim]={0};
    Int_t MC_lastplate[maxdim]={0};
    
    
    EdbVertex *vertex = new EdbVertex();
    EdbTrackP *track = new EdbTrackP();
    
    vtx->Branch("vID",&vID,"vID/I");
    vtx->Branch("vx",&vx,"vx/F");
    vtx->Branch("vy",&vy,"vy/F");
    vtx->Branch("vz",&vz,"vz/F");
    vtx->Branch("vplate",&vplate,"vplate/I");
    vtx->Branch("v_flag",&v_flag,"v_flag/I");
    vtx->Branch("maxaperture",&maxaperture,"maxaperture/F");
    //vtx->Branch("maxrmsthetaspace",&maxrmsthetaspace,"maxrmsthetaspace/F");
    vtx->Branch("probability",&probability,"probability/F");
    vtx->Branch("n",&n,"n/I");
    //track variables (they are array with the number of tracks as size)
    vtx->Branch("IDTrack",&IDTrack,"IDTrack[n]/I");
    vtx->Branch("X",&X,"X[n]/F");
    vtx->Branch("Y",&Y,"Y[n]/F");
    vtx->Branch("Z",&Z,"Z[n]/F");
    vtx->Branch("TX",&TX,"TX[n]/F");
    vtx->Branch("TY",&TY,"TY[n]/F");
    vtx->Branch("Theta",&Theta,"Theta[n]/F");
    vtx->Branch("nseg",&nseg,"nseg[n]/I");
    vtx->Branch("npl",&npl,"npl[n]/I");
    vtx->Branch("nholes",&nholes,"nholes[n]/I");
    vtx->Branch("nseg_S1",&nseg_S1,"nseg_S1[n]/I");
    vtx->Branch("nholes_S1",&nholes_S1,"nholes_S1[n]/I");
    vtx->Branch("plate",&plate," plate[n]/I");
    vtx->Branch("maxgap",&maxgap,"maxgap[n]/I");
    vtx->Branch("incoming",&incoming,"incoming[n]/I");
    //vtx->Branch("MC_IDPart",&MC_IDPart,"MC_IDPart[n]/I");
    vtx->Branch("impactparameter",&impactparameter,"impactparameter[n]/F");
    vtx->Branch("Z_flag",&Z_flag,"Z_flag[n]/I");
    vtx->Branch("MC_Charge_first",&MC_Charge_first,"MC_Charge_first[n]/I");
    vtx->Branch("MC_Charge_last",&MC_Charge_last,"MC_Charge_last[n]/I");
    vtx->Branch("MC_Charge_S2",&MC_Charge_S2,"MC_Charge_S2[n]/I");
    vtx->Branch("MC_evID_first",&MC_evID_first,"MC_evID_first[n]/I");
    vtx->Branch("MC_evID_last",&MC_evID_last,"MC_evID_last[n]/I");
    vtx->Branch("MC_trackID_first",&MC_trackID_first,"MC_trackID_first[n]/I");
    vtx->Branch("MC_trackID_last",&MC_trackID_last,"MC_trackID_last[n]/I");
    vtx->Branch("MC_mother_first",&MC_mother_first,"MC_mother_first[n]/I");
    vtx->Branch("MC_mother_last",&MC_mother_last,"MC_mother_last[n]/I");
    vtx->Branch("MC_pdgcode_first",&MC_pdgcode_first,"MC_pdgcode_first[n]/I");
    vtx->Branch("MC_pdgcode_last",&MC_pdgcode_last,"MC_pdgcode_last[n]/I");

    vtx->Branch("MC_firstplate",&MC_firstplate,"MC_firstplate[n]/I");
    vtx->Branch("MC_lastplate",&MC_lastplate,"MC_lastplate[n]/I");
    
    
    //gAli->Write();
    fvtx->cd();
    TObjArray *varr = new TObjArray(); //vertices to be saved
    //gEVR->Write();
    for(Int_t ivtx=0; ivtx<gEVR->eVTX->GetEntries(); ivtx++){
        vertex=(EdbVertex*)(gEVR->eVTX->At(ivtx));
        if(vertex->Flag()<0) continue; //saving only 'true' vertices in the tree file and in the object
        vID=vertex->ID();
        vx=vertex->VX();
        vy=vertex->VY();
        vz=vertex->VZ();
        n=vertex->N();
        vplate=Get_vtx_plate(vz);
        v_flag=vertex->Flag();
        maxaperture = vertex->MaxAperture();
        probability = vertex->V()->prob();
        //if(n < 4) continue;
        //if (maxaperture <0.01) continue;
        varr->Add(vertex); //true vertex, we can save it
        //adding vertex to list to be saved
        //loop on tracks //now it can be done offline (again)
        for (int itrk = 0; itrk < n; itrk++){
            track = vertex->GetTrack(itrk);
            
            
            IDTrack[itrk] = track->Track(); //track->GetSegment(0)->Track();//track->ID();
            //cout << track->GetSegment(0)->Track() << "\t" << track->ID() << "\t" << track->Theta() << endl;
            nseg[itrk] = track->N();
            npl[itrk] = track->Npl();
            Int_t zpos = vertex->GetVTa(itrk)->Zpos();
            incoming[itrk] = zpos;
            X[itrk] = track->GetSegment(0)->X();
            Y[itrk] = track->GetSegment(0)->Y();
            Z[itrk] = track->GetSegment(0)->Z();
            TX[itrk] = track->TX();
            TY[itrk] = track->TY();
            Theta[itrk] = track->Theta();
            plate[itrk] = track->GetSegment(0)->Plate();
            nholes[itrk] = track->N0();
            maxgap[itrk] = track->CheckMaxGap();
            
            nseg_S1[itrk]=0;
            npl_S1[itrk]=0;
            nholes_S1[itrk]=0;
            if(plate[itrk]<31){
                for (int iseg = 0; iseg < track->N(); iseg++){ //loop on segments
                    if(track->GetSegment(iseg)->Plate()<31){
                        nseg_S1[itrk]++;
                    }
                }
                npl_S1[itrk]=30-plate[itrk]+1;
                nholes_S1[itrk]=npl_S1[itrk]-nseg_S1[itrk];
            }
            impactparameter[itrk] = vertex->GetVTa(itrk)->Imp();
            Z_flag[itrk] = track->Flag();
            MC_Charge_S2[itrk] = -99;
            for(int is=0; is<track->N(); is++){
                if(track->GetSegment(is)->Plate()>=31){
                    if(track->GetSegmentLast()->Plate()-track->GetSegment(is)->Plate()>3){
                    MC_Charge_S2[itrk] = track->GetSegment(is)->W()-70;
                    if(MC_Charge_S2[itrk]>4) MC_Charge_S2[itrk] =4;
                    continue;
                }
                }
            }


            MC_evID_first[itrk] = track->GetSegmentFirst()->MCEvt();
            MC_evID_last[itrk] = track->GetSegmentLast()->MCEvt();
            MC_trackID_first[itrk] = track->GetSegmentFirst()->MCTrack();
            MC_trackID_last[itrk] = track->GetSegmentLast()->MCTrack();
            MC_pdgcode_first[itrk] = track->GetSegmentFirst()->Aid(0);
            MC_pdgcode_last[itrk] = track->GetSegmentLast()->Aid(0);
            MC_mother_first[itrk] = track->GetSegmentFirst()->Aid(1);
            MC_mother_last[itrk] = track->GetSegmentLast()->Aid(1);
            MC_firstplate[itrk] = track->GetSegmentFirst()->Vid(0);
            MC_lastplate[itrk] = track->GetSegmentLast()->Vid(0);
            MC_Charge_last[itrk] = track->GetSegmentLast())->W()-70();
            MC_Charge_first[itrk] = track->GetSegmentFirst())->W()-70();

            
            
        }
        vtx->Fill();
    }
    
    
    //trying to save the selected vertices
    EdbVertexRec *mygEVR = new EdbVertexRec();
    mygEVR->eVTX = varr;
    //using same parameters for vertexrec as the original one
    mygEVR->eDZmax=DZmax;
    mygEVR->eProbMin=ProbMinV;
    mygEVR->eImpMax=ImpMax;
    mygEVR->eUseMom=UseMom;
    mygEVR->eUseSegPar=UseSegPar;
    mygEVR->eQualityMode=QualityMode;
    mygEVR->Write();
    
    vtx->Write();
    fvtx->Close(); //close the file where vertices are saved
    
    //hip->Draw();
    //fclose(tvtx);
    cout<<"Ho salvato i vertici nel file"<<endl;
    
}

//---------------------------------------------------------------------
void td()
{
    // draw all tracks
    
    TObjArray *trarr=gAli->eTracks;
    gStyle->SetPalette(1);
    const char *dsname="display-t";
    ds = EdbDisplay::EdbDisplayExist(dsname);
    if(!ds)  ds=new EdbDisplay(dsname,-50000.,50000.,-50000.,50000.,-4000.,80000.);
    ds->SetVerRec(gEVR);
    ds->SetDrawTracks(4);
    ds->SetArrTr( trarr );
    printf("%d tracks to display\n", trarr->GetEntries() );
    ds->Draw();
}

//---------------------------------------------------------------------
void sd()
{
    // draw all tracks and segments (basetracks)
    
    TObjArray *trarr = gAli->eTracks;
    TObjArray *sarr  = new TObjArray();
    
    EdbSegP *s=0;
    for(int i=0; i<gAli->Npatterns(); i++) {
        EdbPattern *pat = gAli->GetPattern(i);
        for(int j=0; j<pat->N(); j++) {
            s = pat->GetSegment(j);
            if(s->Track()<0)                  // exclude segments already attached to tracks
                sarr->Add(s);
        }
    }
    
    printf("%d tracks to display\n",   trarr->GetEntries() );
    printf("%d segments to display\n", sarr->GetEntries()  );
    
    gStyle->SetPalette(1);
    const char *dsname="display-s";
    ds = EdbDisplay::EdbDisplayExist(dsname);
    if(!ds)  ds=new EdbDisplay(dsname,-50000.,50000.,-50000.,50000.,-4000.,80000.);
    ds->SetVerRec(gEVR);
    ds->SetDrawTracks(4);
    ds->SetArrSegP( sarr );
    ds->SetArrTr( trarr );
    ds->Draw();
}

//---------------------------------------------------------------------
void vd( int trmin, float amin)
{
    // draw vertexes with multiplicity>=trmin, and aperture >= amin
    
    TObjArray *varr = new TObjArray();
    TObjArray *tarr = new TObjArray();
    
    EdbVertex *v=0;
    EdbTrackP *t=0;
    
    int nv = gEVR->Nvtx();
    printf("nv=%d\n",nv);
    if(nv<1) return;
    
    for(int i=0; i<nv; i++) {
        v = (EdbVertex *)(gEVR->eVTX->At(i));
        if(v->Flag()<0)         continue;
        if( v->N()<trmin )                       continue;
        if( v->N()<3 && v->MaxAperture()<amin )  continue;
        
        varr->Add(v);
        for(int j=0; j<v->N(); j++) tarr->Add( v->GetTrack(j) );
    }
    
    gStyle->SetPalette(1);
    
    const char *dsname="display-v";
    ds = EdbDisplay::EdbDisplayExist(dsname);
    if(!ds)  ds=new EdbDisplay(dsname,-100000.,100000.,-100000.,100000.,-40000., 0.);
    ds->SetVerRec(gEVR);
    ds->SetArrTr( tarr );
    printf("%d tracks to display\n", tarr->GetEntries() );
    ds->SetArrV( varr );
    printf("%d vertex to display\n", varr->GetEntries() );
    ds->SetDrawTracks(4);
    ds->SetDrawVertex(1);
    ds->Draw();
}

//---------------------------------------------------------------------
// Experimental function to study segment distance to vertex
void add_proton(EdbVertex* vertex, EdbVertexRec* gEVR){
    //cout<<"start analyzing a vertex"<<endl;
    float track_ipmin = 10000.;
    int ntr  = gEVR->eEdbTracks->GetEntries(); //number of tracks
    //   EdbVTA *vta = 0;
    EdbSegP *potsegment;
    
    bool associated_track; //track near the vertex
    
    for(int itr=0; itr<ntr; itr++)   {
        associated_track = false;
        EdbTrackP *tr = (EdbTrackP*)(gEVR->eEdbTracks->At(itr));
        
        //if (tr->N() > 26){
        
        float iptrack2vertex = vertex->DistTrack(tr,0.);
        if (iptrack2vertex < track_ipmin) track_ipmin = iptrack2vertex;
        
        if (iptrack2vertex < 200.) associated_track = true; //this is not the right track
        
        if (associated_track) cout << "Let's see each segment for this track, track distance = "<<iptrack2vertex<<endl;
        
        for (int iseg = 0; iseg < tr->N(); iseg++){ //loop on segments
            
            potsegment = (EdbSegP*)tr->GetSegment(iseg);
            float ipseg2vertex = vertex->DistSeg(potsegment,0.);
            if (associated_track) cout<<ipseg2vertex<<" "<<potsegment->Z()<<endl;
        }
        //}
        if (associated_track) cout<<endl;
    }
    //hip->Fill(track_ipmin);
}

//---------------------------------------------------------------------

void get_z_pattern()
{
    int np = gAli->Npatterns();
    int layer=np;
    for(int i=0; i<np; i++) {
        EdbPattern *p = gAli->GetPattern(i);
        if(DIRECTION==9){
            Z_LAYER[layer]=p->Z();
            cout << "layer: " << layer << "\t" << Z_LAYER[layer] << endl;
            layer--;
        }
        else if(DIRECTION==1) {
            Z_LAYER[i]=p->Z();
            cout << "layer: " << i << "\t" << Z_LAYER[i] << endl;
        }
    }
}

//---------------------------------------------------------------------

int Get_vtx_plate(float vz){
    
    for(int layer=0;layer<=LASTLAYER;layer++){
        if(vz>Z_LAYER[layer]&&vz<Z_LAYER[layer+1]){
            //cout << "vz: " << vz << "\tlayer " << layer << "\tmin " << Z_LAYER[layer] << "\t" << Z_LAYER[layer+1] << endl;
            return layer;
        }
    }
    
    return LASTLAYER+1;
}

//    dset->GetPlate(21)->Z();

