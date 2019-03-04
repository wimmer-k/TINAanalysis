#include <TFile.h>
#include <TChain.h>
#include <TH2.h>
#include <TEnv.h>
#include <TCutG.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <signal.h>
#include <sys/time.h>

//energy loss, kinematics, and command line interface libraries
#include "Reconstruction.hh"
#include "Kinematics.hh"
#include "CommandLineInterface.hh"

using namespace std;
#define DEBUG 0

#ifndef rad2deg
#define rad2deg                       180./(TMath::Pi())
#endif
#ifndef deg2rad
#define deg2rad                       (TMath::Pi())/180.
#endif

#define NDET 6
#define NADC 5
#define RAWCSI 0

TSpline3* protDete_l2e[NDET][16];
TSpline3* deutDete_l2e[NDET][16];
TSpline3* protTarg_e2r;
TSpline3* protTarg_r2e;
TSpline3* deutTarg_e2r;
TSpline3* deutTarg_r2e;
TSpline3* protFoil_e2r;
TSpline3* protFoil_r2e;
TSpline3* deutFoil_e2r;
TSpline3* deutFoil_r2e;

TCutG* deutCut[NDET];
TCutG* protCut[NDET];
TCutG* csiSCut[NDET];
#ifndef KYUSHU
TCutG* csiLCut[NDET];
#endif

//map strip number to angle
map<int,double> strip2thetamap(char* filename);
//generate energy loss splines
double calcenergyloss(char* detectorsetup);
//pid cuts
void readpidcuts(char* filename);
//signal handling, enables ctrl-c to exit nicely
bool signal_received = false;
void signalhandler(int sig);
double get_time();
int main(int argc, char** argv){
  double time_start = get_time();  
  signal(SIGINT,signalhandler);
  vector<char*> InputFiles;
  char *OutFile = NULL;
  char* detectorsetup = NULL;
  char* pidfile = NULL;
  int LastEvent =-1;
  int Verbose =0;
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-i", "input files", &InputFiles);
  interface->Add("-o", "output file", &OutFile);    
  interface->Add("-d", "detector settings", &detectorsetup);
  interface->Add("-p", "pid cut file", &pidfile);  
  interface->Add("-le", "last event to be read", &LastEvent);  
  interface->Add("-v", "verbose level", &Verbose);  
  interface->CheckFlags(argc, argv);
  //Complain about missing mandatory arguments
  if(InputFiles.size() == 0){
    cout << "No input file given " << endl;
    return 1;
  }
  //set defaults
  if(detectorsetup==NULL)
    detectorsetup = (char*)"/home/wimmer/TINA/settings/detectorsetup_2um.dat";
  if(pidfile==NULL)
    pidfile = (char*)"/home/wimmer/TINA/settings/pidcuts.root";

  TEnv* detsetup = new TEnv(detectorsetup);
  double detectorangle = detsetup->GetValue("Detector.Angle",50.0);
  double foilthick = detsetup->GetValue("Foil.Thickness",36.0);
  double targetthick = detsetup->GetValue("Target.Thickness",2.0);
  double foildensity = detsetup->GetValue("Foil.Density",2.7);
  double targetdensity = detsetup->GetValue("Target.Density",4.5);

  //generate energy loss splines
  double emid = calcenergyloss(detectorsetup);
  //read in pid cuts
  readpidcuts(pidfile);
  //calculate kinematics
  Nucleus *prot = new Nucleus(1,0,massFile);
  Nucleus *deut = new Nucleus(1,1,massFile);
  Nucleus *carb = new Nucleus(6,6,massFile);
  Nucleus *ejec = new Nucleus(6,7,massFile);

  Kinematics* ddkine = NULL;
  Kinematics* ppkine = NULL;
  vector<Kinematics*> dpkine;
  
  ddkine = new Kinematics(carb,deut,deut,carb,emid,0);
  ppkine = new Kinematics(carb,prot,prot,carb,emid,0);
  double energies[4] = {0,3.089,3.684,3.853};
  for(int i=0;i<4;i++)
    dpkine.push_back(new Kinematics(carb,deut,prot,ejec,emid,energies[i]));

  cout << "calculated kinematics " << endl;
  
  cout<<"input file(s):"<<endl;
  for(unsigned int i=0; i<InputFiles.size(); i++){
    cout<<InputFiles[i]<<endl;
  }
  TChain* tr;
  tr = new TChain("tr");
  for(unsigned int i=0; i<InputFiles.size(); i++){
    tr->Add(InputFiles[i]);
  }
  if(tr == NULL){
    cout << "could not find tree \"tr\" in file " << endl;
    for(unsigned int i=0; i<InputFiles.size(); i++){
      cout<<InputFiles[i]<<endl;
    }
    return 3;
  }

  //histograms
  TList *hlist = new TList();
<<<<<<< HEAD

  //detector performance, csi calibration etc
  TH2F* hsienergy_csi[NDET];
  TH2F* hsicorr_csi[NDET];
  TH2F* hsibsen_csi[NDET];
  TH2F* hcal_csi[NDET];
  TH1F* hres_csi[NDET];
  TH2F* hsienergy_theta[NDET];
  TH2F* hsireco_theta[NDET];

  //physics
  TH2F* hen_theta[NDET];
  TH2F* hddexc_theta[NDET];
  TH2F* hdpexc_theta[NDET];
  TH1F* hddexc[NDET];
  TH1F* hdpexc[NDET];
  
=======
  TH2F* hsifsen_csiS[NDET];
  TH2F* hsicorr_csiS[NDET];
  TH2F* hsibsen_csiS[NDET];
  TH2F* hcal_csiS[NDET];
#ifndef KYUSHU
  TH2F* hsifsen_csiL[NDET];
  TH2F* hsicorr_csiL[NDET];
  TH2F* hsibsen_csiL[NDET];
  TH2F* hcal_csiL[NDET];
#endif
  TH2F* hsifsen_theta[NDET];
  TH2F* htotreco_theta[NDET];

  TH2F* htotreco_theta_a;
  TH2F* hsicorr_csi_a;

  //histograms
>>>>>>> 6d3290167d355101c1169de99269e4a20f9d2022
  for(unsigned short d=0;d<NDET;d++){
    if(RAWCSI){
      hsifsen_csiS[d] = new TH2F(Form("hsifsen_csiS_%d",d),Form("hsifsen_csiS_%d",d),512,0,4096,2000,0,20);hlist->Add(hsifsen_csiS[d]);
      hsicorr_csiS[d] = new TH2F(Form("hsicorr_csiS_%d",d),Form("hsicorr_csiS_%d",d),512,0,4096,2000,0,20);hlist->Add(hsicorr_csiS[d]);
      hsibsen_csiS[d] = new TH2F(Form("hsibsen_csiS_%d",d),Form("hsibsen_csiS_%d",d),512,0,4096,2000,0,20);hlist->Add(hsibsen_csiS[d]);
    }
    else{
      hsifsen_csiS[d] = new TH2F(Form("hsifsen_csiS_%d",d),Form("hsifsen_csiS_%d",d),500,0,40,2000,0,20);hlist->Add(hsifsen_csiS[d]);
      hsicorr_csiS[d] = new TH2F(Form("hsicorr_csiS_%d",d),Form("hsicorr_csiS_%d",d),500,0,40,2000,0,20);hlist->Add(hsicorr_csiS[d]);
      hsibsen_csiS[d] = new TH2F(Form("hsibsen_csiS_%d",d),Form("hsibsen_csiS_%d",d),500,0,40,2000,0,20);hlist->Add(hsibsen_csiS[d]);
    }
<<<<<<< HEAD
    hsienergy_theta[d] = new TH2F(Form("hsienergy_theta_%d",d),Form("hsienergy_theta_%d",d),55,25,80,2000,0,20);hlist->Add(hsienergy_theta[d]);
    hsireco_theta[d] = new TH2F(Form("hsireco_theta_%d",d),Form("hsireco_theta_%d",d),55,25,80,2000,0,20);hlist->Add(hsireco_theta[d]);
    hen_theta[d] = new TH2F(Form("hen_theta_%d",d),Form("hen_theta_%d",d),55,25,80,2000,0,20);hlist->Add(hen_theta[d]);
    
    hddexc_theta[d] = new TH2F(Form("hddexc_theta_%d",d),Form("hddexc_theta_%d",d),55,25,80,2000,-1,1);hlist->Add(hddexc_theta[d]);
    hdpexc_theta[d] = new TH2F(Form("hdpexc_theta_%d",d),Form("hdpexc_theta_%d",d),55,25,80,2000,-5,5);hlist->Add(hdpexc_theta[d]);
    hddexc[d] = new TH1F(Form("hddexc_%d",d),Form("hddexc_%d",d),2000,-1,1);hlist->Add(hddexc[d]);
    hdpexc[d] = new TH1F(Form("hdpexc_%d",d),Form("hdpexc_%d",d),2000,-5,5);hlist->Add(hdpexc[d]);
    if(RAWCSI)
      hcal_csi[d] = new TH2F(Form("hcal_csi_%d",d),Form("hcal_csi_%d",d),512,0,4096,2000,0,20);
    else{
      hcal_csi[d] = new TH2F(Form("hcal_csi_%d",d),Form("hcal_csi_%d",d),500,0,20,2000,0,20);
      hres_csi[d] = new TH1F(Form("hres_csi_%d",d),Form("hres_csi_%d",d),4000,-2,2);hlist->Add(hres_csi[d]);
    }
    hlist->Add(hcal_csi[d]);
=======
#ifndef KYUSHU
    if(RAWCSI){
      hsifsen_csiL[d] = new TH2F(Form("hsifsen_csiL_%d",d),Form("hsifsen_csiL_%d",d),512,0,4096,2000,0,20);hlist->Add(hsifsen_csiL[d]);
      hsicorr_csiL[d] = new TH2F(Form("hsicorr_csiL_%d",d),Form("hsicorr_csiL_%d",d),512,0,4096,2000,0,20);hlist->Add(hsicorr_csiL[d]);
      hsibsen_csiL[d] = new TH2F(Form("hsibsen_csiL_%d",d),Form("hsibsen_csiL_%d",d),512,0,4096,2000,0,20);hlist->Add(hsibsen_csiL[d]);
    }
    else{
      hsifsen_csiL[d] = new TH2F(Form("hsifsen_csiL_%d",d),Form("hsifsen_csiL_%d",d),500,0,40,2000,0,20);hlist->Add(hsifsen_csiL[d]);
      hsicorr_csiL[d] = new TH2F(Form("hsicorr_csiL_%d",d),Form("hsicorr_csiL_%d",d),500,0,40,2000,0,20);hlist->Add(hsicorr_csiL[d]);
      hsibsen_csiL[d] = new TH2F(Form("hsibsen_csiL_%d",d),Form("hsibsen_csiL_%d",d),500,0,40,2000,0,20);hlist->Add(hsibsen_csiL[d]);
    }
#endif
#ifdef KYUSHU
    hsifsen_theta[d] = new TH2F(Form("hsifsen_theta_%d",d),Form("hsifsen_theta_%d",d),55,25,80,2000,0,20);hlist->Add(hsifsen_theta[d]);
    htotreco_theta[d] = new TH2F(Form("htotreco_theta_%d",d),Form("htotreco_theta_%d",d),55,25,80,2000,0,40);hlist->Add(htotreco_theta[d]);
#else
    hsifsen_theta[d] = new TH2F(Form("hsifsen_theta_%d",d),Form("hsifsen_theta_%d",d),60,100,160,2000,0,20);hlist->Add(hsifsen_theta[d]);
    htotreco_theta[d] = new TH2F(Form("htotreco_theta_%d",d),Form("htotreco_theta_%d",d),60,100,160,2000,0,40);hlist->Add(htotreco_theta[d]);
#endif
    if(RAWCSI)
      hcal_csiS[d] = new TH2F(Form("hcal_csiS_%d",d),Form("hcal_csiS_%d",d),512,0,4096,2000,0,40);
    else
      hcal_csiS[d] = new TH2F(Form("hcal_csiS_%d",d),Form("hcal_csiS_%d",d),500,0,40,2000,0,40);
    hlist->Add(hcal_csiS[d]);
#ifndef KYUSHU
    if(RAWCSI)
      hcal_csiL[d] = new TH2F(Form("hcal_csiL_%d",d),Form("hcal_csiL_%d",d),512,0,4096,2000,0,40);
    else
      hcal_csiL[d] = new TH2F(Form("hcal_csiL_%d",d),Form("hcal_csiL_%d",d),500,0,40,2000,0,40);
    hlist->Add(hcal_csiL[d]);
#endif
>>>>>>> 6d3290167d355101c1169de99269e4a20f9d2022
  }
  htotreco_theta_a = new TH2F("htotreco_theta_a","htotreco_theta_a",60,100,160,2000,0,40);hlist->Add(htotreco_theta_a);
  hsicorr_csi_a = new TH2F("hsicorr_csi_a","hsicorr_csi_a",500,0,40,2000,0,20);hlist->Add(hsicorr_csi_a);

  TH2F* tdccorr = new TH2F("tdccorr","tdccorr",12,0,12,12,0,12);
  hlist->Add(tdccorr);

  int adc[NADC][32];
  int tdc[64];
  int siring[NDET];
  double sitheta[NDET];
  double sifsen[NDET];
  double sibsen[NDET];
  double csiSen[NDET];
#ifndef KYUSHU
  double csiLen[NDET];
#endif
  double totreco[2][NDET];
  tr->SetBranchAddress("adc",&adc);
  tr->SetBranchAddress("tdc",&tdc);
  tr->SetBranchAddress("siring",&siring);
  tr->SetBranchAddress("sitheta",&sitheta);
  tr->SetBranchAddress("sifsen",&sifsen);
  tr->SetBranchAddress("sibsen",&sibsen);
  tr->SetBranchAddress("totreco",&totreco);
  tr->SetBranchAddress("csiSen",&csiSen);
#ifndef KYUSHU
  tr->SetBranchAddress("csiLen",&csiLen);
#endif
  Double_t nentries = tr->GetEntries();
  cout << nentries << " entries in tree" << endl;
  if(nentries<1)
    return 4;
  if(LastEvent>0)
    nentries = LastEvent;


  int startcsitdc =16;
  long long int nbytes = 0;
  Int_t status;
  for(int i=0; i<nentries;i++){
    if(signal_received){
      break;
    }
    if(i%100000 == 0){
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*i)/nentries<<" % done\t" << 
	(Float_t)i/(time_end - time_start) << " events/s " <<
	(nentries-i)*(time_end - time_start)/(Float_t)i << "s to go \r" << flush;
    }
    for(int c=0;c<32;c++){
      for(int j=0;j<NADC;j++)
	adc[j][c] = 0;
    }
    for(unsigned short t=0;t<64;t++)
      tdc[t] = 0;
    for(unsigned short d=0;d<NDET;d++){
      siring[d] = -1;
      sitheta[d] = NAN;
      sifsen[d] = NAN;
      totreco[0][d] = NAN;
      totreco[1][d] = NAN;
      sibsen[d] = NAN;
      csiSen[d] = NAN;
#ifndef KYUSHU
      csiLen[d] = NAN;
#endif
    }
    if(Verbose>2)
      cout << "getting entry " << i << endl;
    status = tr->GetEvent(i);
    if(Verbose>2)
      cout << "status " << status << endl;
    if(status == -1){
      cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<endl;
      return 5;
    }
    else if(status == 0){
      cerr<<"Error occured, entry "<<i<<" in tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<" doesn't exist"<<endl;
      return 6;
    }
    nbytes += status;
    if(RAWCSI){
#ifdef KYUSHU 
      for(int chan = 0; chan<NDET;chan++)
	csiSen[chan] = adc[3][chan];
#else
      for(int chan = 16; chan<16+NDET;chan++)
	csiSen[chan-16] = adc[3][chan];
      for(int chan = 16+NDET; chan<16+2*NDET;chan++)
	csiLen[chan-16-NDET] = adc[3][chan];
#endif
    }
    for(unsigned short d=0;d<NDET;d++){
<<<<<<< HEAD
      //basic checks and csi calibration
      if(sienergy[d]>0)
	hsienergy_theta[d]->Fill(sitheta[d],sienergy[d]);
      if(sireco[1][d]>0)
	hsireco_theta[d]->Fill(sitheta[d],sireco[1][d]);    
      if(sienergy[d]>0 && csien[d]>0){
	hsienergy_csi[d]->Fill(csien[d],sienergy[d]);
	double alpha = 90.-detectorangle-sitheta[d];
	hsicorr_csi[d]->Fill(csien[d],cos(alpha*deg2rad)*sienergy[d]);
	if(pidCut[d] && pidCut[d]->IsInside(csien[d],sienergy[d])){
	  //cout << d <<",\tE = "<< sienergy[d] << ",\tth = "<< sitheta[d]<< ",\tcsiraw = "<< csien[d] << endl;
	  if(siring[d]>-1){
	    double en = protDete_l2e[d][siring[d]]->Eval(sienergy[d]);
	    //cout << siring[d] << "\t" << sienergy[d] << "\t" << en << endl;
	    hcal_csi[d]->Fill(csien[d],en-sienergy[d]);
	    if(!RAWCSI)
	      hres_csi[d]->Fill(csien[d] - (en-sienergy[d]));
=======
      //cout << sifsen[d] << "\t "<< sibsen[d]<< "\t "<< csiSen[d] << endl;
      if(sifsen[d]>0)
	hsifsen_theta[d]->Fill(sitheta[d],sifsen[d]);
      if(totreco[1][d]>0){
	htotreco_theta[d]->Fill(sitheta[d],totreco[1][d]);    
	htotreco_theta_a->Fill(sitheta[d],totreco[1][d]);    
      }
      if(sifsen[d]>0 && csiSen[d]>0){
	hsifsen_csiS[d]->Fill(csiSen[d],sifsen[d]);
	double alpha;
#ifdef KYUSHU
	alpha = 90.-detectorangle-sitheta[d];
#else
	alpha = -90.-detectorangle+sitheta[d];
#endif
	double ecorr = fabs(cos(alpha*deg2rad))*sifsen[d];
	hsicorr_csiS[d]->Fill(csiSen[d],ecorr);
	hsicorr_csi_a->Fill(csiSen[d],ecorr);
	//if(csiSCut[d] && csiSCut[d]->IsInside(csiSen[d],ecorr)){
	if(csiSCut[d] && csiSCut[d]->IsInside(csiSen[d],sifsen[d])){
	  //cout << d <<",\tE = "<< sifsen[d] << ",\tth = "<< sitheta[d]<< ",\tcsiraw = "<< csien[d] << endl;
	  if(siring[d]>-1){
	    double en = protDete_l2e[d][siring[d]]->Eval(sifsen[d]);
	    //cout << siring[d] << "\t" << sifsen[d] << "\t" << en << endl;
	    hcal_csiS[d]->Fill(csiSen[d],en-sifsen[d]);
	  }
	}
      }//si and csi en >0
      if(sibsen[d]>0 && csiSen[d]>0)
	hsibsen_csiS[d]->Fill(csiSen[d],sibsen[d]);
#ifndef KYUSHU
      if(sifsen[d]>0 && csiLen[d]>0){
	hsifsen_csiL[d]->Fill(csiLen[d],sifsen[d]);
	double alpha;
	alpha = -90.-detectorangle+sitheta[d];
	double ecorr = fabs(cos(alpha*deg2rad))*sifsen[d];
	hsicorr_csiL[d]->Fill(csiLen[d],ecorr);
	hsicorr_csi_a->Fill(csiLen[d],ecorr);
	//if(csiLCut[d] && csiLCut[d]->IsInside(csiLen[d],ecorr)){
	if(csiLCut[d] && csiLCut[d]->IsInside(csiLen[d],sifsen[d])){
	  //cout << d <<",\tE = "<< sifsen[d] << ",\tth = "<< sitheta[d]<< ",\tcsiraw = "<< csien[d] << endl;
	  if(siring[d]>-1){
	    double en = protDete_l2e[d][siring[d]]->Eval(sifsen[d]);
	    //cout << siring[d] << "\t" << sifsen[d] << "\t" << en << endl;
	    hcal_csiL[d]->Fill(csiLen[d],en-sifsen[d]);
>>>>>>> 6d3290167d355101c1169de99269e4a20f9d2022
	  }
	}
      }//si and csi en >0
      if(sibsen[d]>0 && csiLen[d]>0)
	hsibsen_csiL[d]->Fill(csiLen[d],sibsen[d]);
#endif
    }//ndet
    for(int c1=startcsitdc;c1<startcsitdc+12;c1++){
      for(int c2=startcsitdc;c2<startcsitdc+12;c2++){
	if(tdc[c1]>8e4 && tdc[c1]< 12e4 && tdc[c2]>8e4 && tdc[c2]< 12e4)
	  tdccorr->Fill(c1-startcsitdc,c2-startcsitdc);
      }
<<<<<<< HEAD
      if(sibsen[d]>0 && csien[d]>0)
	hsibsen_csi[d]->Fill(csien[d],sibsen[d]);
      //physics
      //pid cuts
      //deuteron, punch protons special
      //everything else proton
      //reconstruct energy losses in foils and target
      //calculated Qvalue
      if(isnan(sienergy[d]))
	continue;
      int pid = -1;
      double ene = sienergy[d];
      double exc = sqrt(-1);
      if(deutCut[d] && deutCut[d]->IsInside(sitheta[d],sienergy[d])){//deuteron
	pid = 2;
	//reconstruct energy loss in the al foil
	double alpha = 90.-detectorangle-sitheta[d];
	double althick = foilthick*2.70*0.1/cos(alpha*deg2rad);
	double range = deutFoil_e2r->Eval(ene);
	ene = deutFoil_r2e->Eval(range+althick);
	//reconstruct energy loss in half the target
	double tithick = targetthick/2*4.50*0.1/cos(sitheta[d]*deg2rad);
	range = deutTarg_e2r->Eval(ene);
	ene = deutTarg_r2e->Eval(range+tithick);
	
	TVector3 recodir(0,0,1);
	recodir.SetTheta(sitheta[d]*deg2rad);
	TLorentzVector recoLV(recodir, ene*1000+deut->GetMass()*1000);
	if(recoLV.Mag()>0)
	  recoLV.SetRho( sqrt( (ene+deut->GetMass())*(ene+deut->GetMass()) - deut->GetMass()*deut->GetMass() )*1000 );
	if(ddkine)
	  exc = ddkine->GetExcEnergy(recoLV)/1000; // to MeV
	hddexc_theta[d]->Fill(sitheta[d],exc);
	hddexc[d]->Fill(exc);
      }//deuteron
      else{
	if(protCut[d] && protCut[d]->IsInside(sitheta[d],sienergy[d])){//punch through proton
	  pid = 1;
	  if(siring[d]>-1){
	    ene = protDete_l2e[d][siring[d]]->Eval(sienergy[d]);
	  }
	}//punch through proton
	else
	  pid = 0;
	
	//reconstruct energy loss in the al foil
	double alpha = 90.-detectorangle-sitheta[d];
	double althick = foilthick*2.70*0.1/cos(alpha*deg2rad);
	double range = deutFoil_e2r->Eval(ene);
	ene = deutFoil_r2e->Eval(range+althick);
	//reconstruct energy loss in half the target
	double tithick = targetthick/2*4.50*0.1/cos(sitheta[d]*deg2rad);
	range = deutTarg_e2r->Eval(ene);
	ene = deutTarg_r2e->Eval(range+tithick);

	TVector3 recodir(0,0,1);
	recodir.SetTheta(sitheta[d]*deg2rad);
	TLorentzVector recoLV(recodir, ene*1000+prot->GetMass()*1000);
	if(recoLV.Mag()>0)
	  recoLV.SetRho( sqrt( (ene+prot->GetMass())*(ene+prot->GetMass()) - prot->GetMass()*prot->GetMass() )*1000 );
	if(dpkine.at(0))
	  exc = dpkine.at(0)->GetExcEnergy(recoLV)/1000; // to MeV
	hdpexc_theta[d]->Fill(sitheta[d],exc);
	hdpexc[d]->Fill(exc);
      }
      if(pid==-1)
	continue;
      hen_theta[d]->Fill(sitheta[d],ene);
=======
>>>>>>> 6d3290167d355101c1169de99269e4a20f9d2022
    }
  }//nevents
  cout << endl;
  cout << nbytes << " bytes read ("<<nbytes/1024/1024<<" MB)" << endl; 
  double time_end = get_time();
  cout << "Loop Run time: " << time_end - time_start << " s." << endl;
  cout << "creating outputfile " << OutFile << endl;
  TFile* ofile = new TFile(OutFile,"recreate");
  ofile->cd();
  TH1F* h1;
  TH2F* h2;
  TIter next(hlist);
  while( (h1 = (TH1F*)next()) ){
    if(h1->GetEntries()>0)
      h1->Write("",TObject::kOverwrite);
  }
  while( (h2 = (TH2F*)next()) ){
    if(h2->GetEntries()>0)
      h2->Write("",TObject::kOverwrite);
  }
  /*
  for(unsigned short d=0;d<NDET;d++){
    protDete_l2e[d][0]->Write(Form("ploss2energy%d_%d",d,0),TObject::kOverwrite);
    deutDete_l2e[d][0]->Write(Form("dloss2energy%d_%d",d,0),TObject::kOverwrite);
    protDete_l2e[d][15]->Write(Form("ploss2energy%d_%d",d,15),TObject::kOverwrite);
    deutDete_l2e[d][15]->Write(Form("dloss2energy%d_%d",d,15),TObject::kOverwrite);
  }
  */
  ofile->Close();
  time_end = get_time();
  cout << "Total Run time: " << time_end - time_start << " s." << endl;
  return 0;

}
void signalhandler(int sig){
  if(sig == SIGINT)
    signal_received = true;
}
double get_time(){  
  struct timeval t;  
  gettimeofday(&t, NULL);  
  double d = t.tv_sec + (double) t.tv_usec/1000000;  
  return d;  
}  
map<int,double> strip2thetamap(char* filename){
  ifstream infile;
  infile.open(filename);
  map<int,double> m;
  if(!infile.is_open()){
    cerr << "channel mapping not found, using default" << endl;
    for(int i=0;i<16;i++){
      m[i] = i;
    }
    return m;
  }
  infile.ignore(1000,'\n');
  for(int i=0;i<16;i++){
    int ch,st;
    double th;
    infile >> ch >> st >> th;
    infile.ignore(1000,'\n');
    m[st] = th;
  }
  return m;
}

double calcenergyloss(char* detectorsetup){
  double detectorthick[NDET];
  TEnv* detsetup = new TEnv(detectorsetup);
  double detectorangle = detsetup->GetValue("Detector.Angle",50.0);
  for(int i=0;i<NDET;i++){
    detectorthick[i] = detsetup->GetValue(Form("Detector.%d",i),300.0);
  }
  string channelmap = detsetup->GetValue("Channel.Mapping", "channelmapping.dat");
  //initalize mapping from ADC channel to strip number and angle
  map<int, double> strmap = strip2thetamap((char*)channelmap.c_str());
  for(int i=0;i<16;i++){
    cout << i << "\t" << strmap[i] << "\t" << endl;
  }
  Nucleus *si = new Nucleus(14,14);
  Compound *dete = new Compound(si);
  
  Nucleus *prot = new Nucleus(1,0);
  Nucleus *deut = new Nucleus(1,1);

  Reconstruction *protDete = new Reconstruction(prot, dete);
  Reconstruction *deutDete = new Reconstruction(deut, dete);
  double protDete_range;
  double deutDete_range;
  double detedensity = 2.32; //g/cm^3 
  for(int d=0;d<NDET;d++){
    detectorthick[d] *= detedensity*0.1;
    cout << "calculating energy loss of " <<  prot->GetSymbol() << ", " << deut->GetSymbol() << " in detector " << d << " d = " << detectorthick[d] << " mg/cm^2" << endl;
    for(int s=0;s<16;s++){
      double alpha;
#ifdef KYUSHU
	alpha = 90.-detectorangle-strmap[s];
#else
	alpha = -90.-detectorangle+strmap[s];
#endif
      if(DEBUG)
	cout << s << "\t" << alpha << "\t" << cos(alpha*deg2rad) << endl;
      protDete->SetTargetThickness(detectorthick[d]/fabs(cos(alpha*deg2rad)));
      protDete_l2e[d][s] = protDete->EnergyLoss2Energy(50,0.1);
      deutDete->SetTargetThickness(detectorthick[d]/fabs(cos(alpha*deg2rad)));
      deutDete_l2e[d][s] = deutDete->EnergyLoss2Energy(50,0.1);
    }
  }


  //energy loss in target and foils
  double foilthick = detsetup->GetValue("Foil.Thickness",36.0);
  double targetthick = detsetup->GetValue("Target.Thickness",2.0);
  double foildensity = detsetup->GetValue("Foil.Density",2.7);
  double targetdensity = detsetup->GetValue("Target.Density",4.5);
  Nucleus *ti = new Nucleus(22,26,massFile);
  Compound *target = new Compound(ti);
  //Compound *target = new Compound("2.0DTI");
  Nucleus *al = new Nucleus(13,14,massFile);
  Compound *foil = new Compound(al);
  Reconstruction *protTarg = new Reconstruction(prot, target);
  Reconstruction *deutTarg = new Reconstruction(deut, target);
 
  double protTarg_range;
  double deutTarg_range;
  targetthick *= targetdensity*0.1;
  protTarg->SetTargetThickness(targetthick/2);
  cout << "calculating energy loss of " <<  prot->GetSymbol() << " in half target " << targetthick/2 << " mg/cm^2 " ;
  protTarg_e2r = protTarg->Energy2Range(50,0.1);
  protTarg_r2e = protTarg->Range2Energy(50,0.1);
  protTarg_range = protTarg_e2r->Eval(10);
  cout << "10 MeV proton range " << protTarg_range << " mg/cm^2" << endl;

  deutTarg->SetTargetThickness(targetthick/2);
  cout << "calculating energy loss of " <<  deut->GetSymbol() << " in half target " << targetthick/2 << " mg/cm^2 " ;
  deutTarg_e2r = deutTarg->Energy2Range(50,0.1);
  deutTarg_r2e = deutTarg->Range2Energy(50,0.1);
  deutTarg_range = deutTarg_e2r->Eval(10);
  cout << "10 MeV deuteron range " << deutTarg_range << " mg/cm^2" << endl;

  Reconstruction *protFoil = new Reconstruction(prot, foil);
  Reconstruction *deutFoil = new Reconstruction(deut, foil);
  double protFoil_range;
  double deutFoil_range;
  foilthick *= foildensity*0.1;
  protFoil->SetTargetThickness(foilthick);
  cout << "calculating energy loss of " <<  prot->GetSymbol() << " in foil " << foilthick << " mg/cm^2 " ;
  protFoil_e2r = protFoil->Energy2Range(50,0.1);
  protFoil_r2e = protFoil->Range2Energy(50,0.1);
  protFoil_range = protFoil_e2r->Eval(10);
  cout << "10 MeV proton range " << protFoil_range << " mg/cm^2" << endl;

  deutFoil->SetTargetThickness(foilthick);
  cout << "calculating energy loss of " <<  deut->GetSymbol() << " in foil " << foilthick << " mg/cm^2 " ;
  deutFoil_e2r = deutFoil->Energy2Range(50,0.1);
  deutFoil_r2e = deutFoil->Range2Energy(50,0.1);
  deutFoil_range = deutFoil_e2r->Eval(10);
  cout << "10 MeV deuteron range " << deutFoil_range << " mg/cm^2" << endl;

  double ebeam = detsetup->GetValue("Beam.Energy",20);
  Nucleus *carb = new Nucleus(6,6,massFile);
  Reconstruction *carbTarg = new Reconstruction(carb, target);
  TSpline3* carbTarg_e2r = carbTarg->Energy2Range(50,0.1);
  TSpline3* carbTarg_r2e = carbTarg->Range2Energy(50,0.1);
  double range = carbTarg_e2r->Eval(ebeam);
  double emid = carbTarg_r2e->Eval(range-targetthick/2);
  cout << "calculating energy of " <<  carb->GetSymbol() << " in middle of " << targetthick/2 << " mg/cm^2: " << emid << endl;

  return emid;

}
void readpidcuts(char* filename){
  for(int i=0;i<NDET;i++){
    deutCut[i] = NULL;
    protCut[i] = NULL;
    csiSCut[i] = NULL;
#ifndef KYUSHU
    csiLCut[i] = NULL;
#endif
  }
  TFile* fc = new TFile(filename);
  for(int i=0;i<NDET;i++){
    deutCut[i] = (TCutG*)fc->Get(Form("deut%d",i));
    protCut[i] = (TCutG*)fc->Get(Form("prot%d",i));
    csiSCut[i] = (TCutG*)fc->Get(Form("csiS%d",i));
#ifndef KYUSHU
    csiLCut[i] = (TCutG*)fc->Get(Form("csiL%d",i));
#endif
  }
}

