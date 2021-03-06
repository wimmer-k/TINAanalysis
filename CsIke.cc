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

TCutG* deutCut[NDET];
TCutG* protCut[NDET];
TCutG* csiSCut[NDET];
#ifndef KYUSHU
TCutG* csiLCut[NDET];
#endif

//map strip number to angle
map<int,double> strip2thetamap(char* filename);
//generate energy loss splines
void calcenergyloss(char* detectorsetup);
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
  //generate energy loss splines
  calcenergyloss(detectorsetup);
  //read in pid cuts
  readpidcuts(pidfile);

  
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

  TList *hlist = new TList();
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
    if(i%10000 == 0){
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
  for(unsigned short d=0;d<NDET;d++){
    protDete_l2e[d][0]->Write(Form("ploss2energy%d_%d",d,0),TObject::kOverwrite);
    deutDete_l2e[d][0]->Write(Form("dloss2energy%d_%d",d,0),TObject::kOverwrite);
    protDete_l2e[d][15]->Write(Form("ploss2energy%d_%d",d,15),TObject::kOverwrite);
    deutDete_l2e[d][15]->Write(Form("dloss2energy%d_%d",d,15),TObject::kOverwrite);
  }
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


void calcenergyloss(char* detectorsetup){
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
