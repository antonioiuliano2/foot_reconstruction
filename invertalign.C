//inverting affine transformations
void invertalign(){
 EdbScanProc *sproc = new EdbScanProc();
 sproc->eProcDirClient = "..";
 //currently available transformation
 int fromid[4] = {1,21,0,1000};
 int toid[4] = {1,20,0,1000};   
 EdbAffine2D affine;
 float dz;
 //reading file with current aff file
 sproc->GetAffZ(affine,dz,fromid,toid);
 //inverting the affine transformation, now affine contains the inverted one
 affine.Print();
 affine.Invert();
 affine.Print();
 //building layer with inverted transformation;
 EdbLayer l;
 l.SetZcorr(-1. * dz); //inverting dz
 l.SetAffXY(affine.A11(),affine.A12(),affine.A21(),affine.A22(),affine.B1(),affine.B2());
 //saving the output
 TString outputfilename;
 sproc->MakeAffName(outputfilename,toid,fromid,"aff.par");
 sproc->UpdateAFFPar(toid,fromid,l);
 cout<<outputfilename.Data()<<endl;
}
