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
#define NSTRIPS 16
#define NADC 5
#define RAWCSI 0
#define S1PADS 2 // has to be same as in Turner

TSpline3* protDete_l2e[NDET][16];
TSpline3* deutDete_l2e[NDET][16];

double kine_energies[6] = {0,2,4,6,8,10}; // excitation energies for dp
Kinematics* ddkine;
Kinematics* ppkine;
vector<Kinematics*> dpkine;

TCutG* deutCut[NDET];
TCutG* protCut[NDET];
TCutG* kineCut;
TCutG* csiSCut[NDET];
#ifndef KYUSHU
TCutG* csiLCut[NDET];
#endif

//map strip number to angle
map<int,double> strip2thetamap(char* filename);
Double_t stripDist[NSTRIPS];
Double_t stripTheta[NSTRIPS];

//generate energy loss splines
void calckinematics(double midtarget, char* settingsfile);
//void calcenergyloss(char* detectorsetup);
//pid cuts
void readpidcuts(char* filename);

//signal handling, enables ctrl-c to exit nicely
bool signal_received = false;
void signalhandler(int sig);
double get_time();

Bool_t InsideEllipse(Double_t x_c, Double_t y_c, Double_t R_x, Double_t R_y, Double_t a_r, Double_t data_x, Double_t data_y);


int main(int argc, char** argv){
  double time_start = get_time();  
  signal(SIGINT,signalhandler);
  vector<char*> InputFiles;
  char *OutFile = NULL;
  char* settingsfile = NULL;
  //char* detectorsetup = NULL;
  char* pidfile = NULL;
  int LastEvent =-1;
  int Verbose =0;
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-i", "input files", &InputFiles);
  interface->Add("-o", "output file", &OutFile);    
  interface->Add("-s", "settings, beam, detector, calibration", &settingsfile);
  //interface->Add("-d", "detector settings", &detectorsetup);
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
  //if(detectorsetup==NULL)
  //  detectorsetup = (char*)"/home/oedo0/TINAroot/settings/detector_positions_day0.dat.dat";
  if(pidfile==NULL)
    pidfile = (char*)"/home/oedo0/TINAroot/settings/pidcuts_cal_CORREN.root";
  if(settingsfile==NULL){
    cout << "Error: no settings file given" << endl;
    return -99;
  }
  TEnv* settings = new TEnv(settingsfile);

  double detectorangle = settings->GetValue("Detector.Angle",50.0);

  //TEnv* detsetup = new TEnv(detectorsetup);
  //double detectorangle = detsetup->GetValue("Detector.Angle",50.0);
  double detectorphi[NDET];
  for(Int_t d = 0; d < NDET; d++){
    //detectorphi[d] = detsetup->GetValue(Form("Detector.Phi.%d", d), 0.0);
    detectorphi[d] = settings->GetValue(Form("Detector.Phi.%d", d), 0.0);
  }
  //generate energy loss splines
  //calcenergyloss(detectorsetup);
  //read in pid cuts
  readpidcuts(pidfile);

  double ebeam = settings->GetValue("Beam.Energy",20.);
  int Abeam = settings->GetValue("Beam.A",99);
  ebeam*=Abeam;
  calckinematics(ebeam,settingsfile); // 20.0 MeV/u at center of target assumed
  
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
  

  // define target: center and radius 
  Double_t bcx=settings->GetValue("Target.Pos.X", 0.0);
  Double_t bcy=settings->GetValue("Target.Pos.Y", 0.0);
  Double_t bcz=settings->GetValue("Target.Pos.Z", 0.0);
  Double_t brx=settings->GetValue("Target.Radius.X", 10.0);
  Double_t bry=settings->GetValue("Target.Radius.Y", 10.0);

  //// ic energy cut -- is done in Turner already
  //Double_t s1ic_lowercut[S1PADS], s1ic_uppercut[S1PADS];
  //for(int p=0;p<S1PADS;p++){
  //  s1ic_lowercut[p] = settings->GetValue(Form("S1.IC.PAD%d.Lower.Cut",p), 0.0);
  //  s1ic_uppercut[p] = settings->GetValue(Form("S1.IC.PAD%d.Upper.Cut",p), 999999.0);
  //}



  TList *hlist = new TList();
  TH2F* hsifsen_csiS[NDET];
  TH2F* hsicorr_csiS[NDET];
  TH2F* hsibsen_csiS[NDET];
  //TH2F* hcal_csiS[NDET];
#ifndef KYUSHU
  TH2F* hsifsen_csiL[NDET];
  TH2F* hsicorr_csiL[NDET];
  TH2F* hsibsen_csiL[NDET];
  //TH2F* hcal_csiL[NDET];
#endif

  TH2F* hsifsen_sitheta[NDET];
  TH2F* htr_sitheta[NDET];
  
  
  TH2F* hsifsen_csi_all;
  TH2F* hsicorr_csi_all;
  TH2F* hsicorr_csi_beam_pos_pad0;

  TH2F* htarget_xy_all;
  TH2F* htarget_xy_poscut;
  TH2F* htarget_xy_beam;
  TH2F* htarget_xy_beam_pos;
  TH2F* htarget_xy_beam_pos_pad0;
  
  TH2F* htr_sitheta_all;
  TH2F* htr_theta_all;
  TH2F* htr_theta_beam;
  TH2F* htr_theta_beam_pos; // cut on incoming beam and position on target
  TH2F* htr_theta_beam_pos_kc; // cut on kinematics plot, rough background suppression
  TH2F* htr_theta_beam_pos_pad0; 
  TH2F* htr_theta_beam_pos_pad1; 
  TH2F* htr_theta_beam_pos_pad01; 
  //TH2F* htr_theta_beam_pos_s1t0; // time cut 0
  //TH2F* htr_theta_beam_pos_s1t1; // time cut 1

  
  TH2F* hexcen_theta_all;
  TH2F* hexcen_totreco_all;
  TH2F* hexcen_theta_beam;
  TH2F* hexcen_theta_beam_pos;
  TH2F* hexcen_theta_beam_pos_pad0; 
  TH2F* hexcen_theta_beam_pos_pad1; 
  TH2F* hexcen_theta_beam_pos_pad01;

  TH1F* hec_all;
  TH1F* hec_beam;
  TH1F* hec_beam_pos;
  TH1F* hec_beam_pos_kc;
  TH1F* hec_beam_pos_pad0; 
  TH1F* hec_beam_pos_pad1; 
  TH1F* hec_beam_pos_pad01;
  TH1F* hec_beam_pos_tpzero;
  TH1F* hec_beam_pos_tpone;
  TH1F* hec_beam_pos_tpzeroone;
  TH1F* hec_beam_pos_tpmone;

  TH2F* hec_theta_all;
  TH2F* hec_theta_beam;
  TH2F* hec_theta_beam_pos;
  TH2F* hec_theta_beam_pos_kc;
  TH2F* hec_theta_beam_pos_pad0;
  TH2F* hec_theta_beam_pos_pad1;
  TH2F* hec_theta_beam_pos_pad01; // geometric mean of pad 0 and 1
  TH2F* hec_theta_beam_pos_tpzero;
  TH2F* hec_theta_beam_pos_tpone;
  TH2F* hec_theta_beam_pos_tpzeroone;
  TH2F* hec_theta_beam_pos_tpmone;
  TH2F* hec_x_all;
  TH2F* hec_x_beam;
  TH2F* hec_x_beam_pos;
  TH2F* hec_y_all;
  TH2F* hec_y_beam;
  TH2F* hec_y_beam_pos;
  TH2F* hec_anglea_all;
  TH2F* hec_anglea_beam;
  TH2F* hec_anglea_beam_pos;
  TH2F* hec_angleb_all;
  TH2F* hec_angleb_beam;
  TH2F* hec_angleb_beam_pos;
  TH2F* hec_totreco_all;
  TH2F* hec_totreco_beam;
  TH2F* hec_totreco_beam_pos;
  TH2F* hec_excen_all;
  TH2F* hec_excen_beam;
  TH2F* hec_excen_beam_pos;

  Int_t excen_bins=300;
  Float_t excen_binmin=-20.0, excen_binmax=40.0;
  Int_t theta_bins=60;
  Float_t theta_binmin=100.0, theta_binmax=160.0;
  Int_t targetx_bins=400;
  Float_t targetx_binmin=-100.0, targetx_binmax=100.0;
  Int_t targety_bins=200;
  Float_t targety_binmin=-50.0, targety_binmax=50.0;
  Int_t anglea_bins=400;
  Float_t anglea_binmin=-100.0, anglea_binmax=100.0;
  Int_t angleb_bins=200;
  Float_t angleb_binmin=-50.0, angleb_binmax=50.0;

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
    hsifsen_sitheta[d] = new TH2F(Form("hsifsen_sitheta_%d",d),Form("hsifsen_sitheta_%d",d),55,25,80,2000,0,20);hlist->Add(hsifsen_sitheta[d]);
    htr_sitheta[d] = new TH2F(Form("htr_sitheta_%d",d),Form("htr_sitheta_%d",d),55,25,80,2000,0,40);hlist->Add(htr_sitheta[d]);
#else
    hsifsen_sitheta[d] = new TH2F(Form("hsifsen_theta_%d",d),Form("hsifsen_theta_%d",d),theta_bins, theta_binmin, theta_binmax,2000,0,20);hlist->Add(hsifsen_sitheta[d]);
    htr_sitheta[d] = new TH2F(Form("htr_theta_%d",d),Form("htr_theta_%d",d),theta_bins, theta_binmin, theta_binmax,2000,0,40);hlist->Add(htr_sitheta[d]);
#endif
//    if(RAWCSI)
//      hcal_csiS[d] = new TH2F(Form("hcal_csiS_%d",d),Form("hcal_csiS_%d",d),512,0,4096,2000,0,40);
//    else
//      hcal_csiS[d] = new TH2F(Form("hcal_csiS_%d",d),Form("hcal_csiS_%d",d),500,0,40,2000,0,40);
//    hlist->Add(hcal_csiS[d]);
//#ifndef KYUSHU
//    if(RAWCSI)
//      hcal_csiL[d] = new TH2F(Form("hcal_csiL_%d",d),Form("hcal_csiL_%d",d),512,0,4096,2000,0,40);
//    else
//      hcal_csiL[d] = new TH2F(Form("hcal_csiL_%d",d),Form("hcal_csiL_%d",d),500,0,40,2000,0,40);
//    hlist->Add(hcal_csiL[d]);
//#endif
  }
  hsifsen_csi_all = new TH2F("hsifsen_csi_all","hsifsen_csi_all",500,0,40,2000,0,20);hlist->Add(hsifsen_csi_all);
  hsicorr_csi_all = new TH2F("hsicorr_csi_all","hsicorr_csi_all",500,0,40,2000,0,20);hlist->Add(hsicorr_csi_all);
  hsicorr_csi_beam_pos_pad0 = new TH2F("hsicorr_csi_beam_pos_pad0","hsicorr_csi_beam_pos_pad0",500,0,40,2000,0,20);hlist->Add(hsicorr_csi_beam_pos_pad0);

#ifndef KYUSHU
  htarget_xy_all = new TH2F("htarget_xy_all","htarget_xy_all",targetx_bins, targetx_binmin, targetx_binmax,targety_bins, targety_binmin, targety_binmax);hlist->Add(htarget_xy_all);
  htarget_xy_poscut = new TH2F("htarget_xy_poscut","htarget_xy_poscut",targetx_bins, targetx_binmin, targetx_binmax,targety_bins, targety_binmin, targety_binmax);hlist->Add(htarget_xy_poscut);
  htarget_xy_beam = new TH2F("htarget_xy_beam","htarget_xy_beam",targetx_bins, targetx_binmin, targetx_binmax,targety_bins, targety_binmin, targety_binmax);hlist->Add(htarget_xy_beam);
  htarget_xy_beam_pos = new TH2F("htarget_xy_beam_pos","htarget_xy_beam_pos",targetx_bins, targetx_binmin, targetx_binmax,targety_bins, targety_binmin, targety_binmax);hlist->Add(htarget_xy_beam_pos);
  htarget_xy_beam_pos_pad0 = new TH2F("htarget_xy_beam_pos_pad0","htarget_xy_beam_pos_pad0",targetx_bins, targetx_binmin, targetx_binmax,targety_bins, targety_binmin, targety_binmax);hlist->Add(htarget_xy_beam_pos_pad0);
#endif

#ifdef KYUSHU
  
  //htr_theta_all = new TH2F("htr_theta_all","htr_theta_all",60,100,160,200,0,40);hlist->Add(htr_theta_all);
  //hexcen_theta_all = new TH2F("hexcen_theta_all","hexcen_theta_all",60,100,160,200,-20,20);hlist->Add(hexcen_theta_all);

#else

  htr_sitheta_all = new TH2F("htr_sitheta_all","htr_sitheta_all",theta_bins,theta_binmin,theta_binmax,200,0,40);hlist->Add(htr_sitheta_all);
  
  htr_theta_all = new TH2F("htr_theta_all","htr_theta_all",theta_bins,theta_binmin,theta_binmax,200,0,40);hlist->Add(htr_theta_all);
  htr_theta_beam = new TH2F("htr_theta_beam","htr_theta_beam",theta_bins,theta_binmin,theta_binmax,200,0,40);hlist->Add(htr_theta_beam);
  htr_theta_beam_pos = new TH2F("htr_theta_beam_pos","htr_theta_beam_pos",theta_bins,theta_binmin,theta_binmax,200,0,40);hlist->Add(htr_theta_beam_pos);
  htr_theta_beam_pos_kc = new TH2F("htr_theta_beam_pos_kc","htr_theta_beam_pos_kc",theta_bins,theta_binmin,theta_binmax,200,0,40);hlist->Add(htr_theta_beam_pos_kc);
  htr_theta_beam_pos_pad0 = new TH2F("htr_theta_beam_pos_pad0","htr_theta_beam_pos_pad0",theta_bins,theta_binmin,theta_binmax,200,0,40);hlist->Add(htr_theta_beam_pos_pad0);
  htr_theta_beam_pos_pad1 = new TH2F("htr_theta_beam_pos_pad1","htr_theta_beam_pos_pad1",theta_bins,theta_binmin,theta_binmax,200,0,40);hlist->Add(htr_theta_beam_pos_pad1);
  htr_theta_beam_pos_pad01 = new TH2F("htr_theta_beam_pos_pad01","htr_theta_beam_pos_pad01",theta_bins,theta_binmin,theta_binmax,200,0,40);hlist->Add(htr_theta_beam_pos_pad01);

  //htr_theta_beam_pos_s1t0 = new TH2F("htr_theta_beam_pos_s1t0","htr_theta_beam_pos_s1t0",theta_bins,theta_binmin,theta_binmax,200,0,40);hlist->Add(htr_theta_beam_pos_s1t0);
  //htr_theta_beam_pos_s1t1 = new TH2F("htr_theta_beam_pos_s1t1","htr_theta_beam_pos_s1t1",theta_bins,theta_binmin,theta_binmax,200,0,40);hlist->Add(htr_theta_beam_pos_s1t1);

  hexcen_theta_all = new TH2F("hexcen_theta_all","hexcen_theta_all",theta_bins,theta_binmin,theta_binmax,excen_bins,excen_binmin,excen_binmax);hlist->Add(hexcen_theta_all);
  hexcen_theta_beam = new TH2F("hexcen_theta_beam","hexcen_theta_beam",theta_bins,theta_binmin,theta_binmax,excen_bins,excen_binmin,excen_binmax);hlist->Add(hexcen_theta_beam);
  hexcen_theta_beam_pos = new TH2F("hexcen_theta_beam_pos","hexcen_theta_beam_pos",theta_bins,theta_binmin,theta_binmax,excen_bins,excen_binmin,excen_binmax);hlist->Add(hexcen_theta_beam_pos);
  hexcen_theta_beam_pos_pad0 = new TH2F("hexcen_theta_beam_pos_pad0","hexcen_theta_beam_pos_pad0",theta_bins,theta_binmin,theta_binmax,excen_bins,excen_binmin,excen_binmax);hlist->Add(hexcen_theta_beam_pos_pad0);
  hexcen_theta_beam_pos_pad1 = new TH2F("hexcen_theta_beam_pos_pad1","hexcen_theta_beam_pos_pad1",theta_bins,theta_binmin,theta_binmax,excen_bins,excen_binmin,excen_binmax);hlist->Add(hexcen_theta_beam_pos_pad1);
  hexcen_theta_beam_pos_pad01 = new TH2F("hexcen_theta_beam_pos_pad01","hexcen_theta_beam_pos_pad01",theta_bins,theta_binmin,theta_binmax,excen_bins,excen_binmin,excen_binmax);hlist->Add(hexcen_theta_beam_pos_pad01);
  
  //hexcen_corr = new TH1F("hexcen_corr","hexcen_corr",600,-20,40);hlist->Add(hexcen_corr);
#endif

  hexcen_totreco_all = new TH2F("hexcen_totreco_all","hexcen_totreco_all", 200,0,40,excen_bins,excen_binmin,excen_binmax);hlist->Add(hexcen_totreco_all);
  
  hec_all = new TH1F("hec_all","hec_all",excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_all);	
  hec_beam = new TH1F("hec_beam","hec_beam",excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_beam);	
  hec_beam_pos = new TH1F("hec_beam_pos","hec_beam_pos",excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_beam_pos);	
  hec_beam_pos_kc = new TH1F("hec_beam_pos_kc","hec_beam_pos_kc",excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_beam_pos_kc);	
  hec_beam_pos_pad0 = new TH1F("hec_beam_pos_pad0","hec_beam_pos_pad0",excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_beam_pos_pad0);	
  hec_beam_pos_pad1 = new TH1F("hec_beam_pos_pad1","hec_beam_pos_pad1",excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_beam_pos_pad1);	
  hec_beam_pos_pad01 = new TH1F("hec_beam_pos_pad01","hec_beam_pos_pad01",excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_beam_pos_pad01);	
  hec_beam_pos_tpzero = new TH1F("hec_beam_pos_tpzero","hec_beam_pos_tpzero",excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_beam_pos_tpzero);	
  hec_beam_pos_tpone = new TH1F("hec_beam_pos_tpone","hec_beam_pos_tpone",excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_beam_pos_tpone);	
  hec_beam_pos_tpzeroone	= new TH1F("hec_beam_pos_tpzeroone","hec_beam_pos_tpzeroone",excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_beam_pos_tpzeroone);	
  hec_beam_pos_tpmone = new TH1F("hec_beam_pos_tpmone","hec_beam_pos_tpmone",excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_beam_pos_tpmone);	
  
  hec_theta_all = new TH2F("hec_theta_all","hec_theta_all",theta_bins,theta_binmin,theta_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_theta_all);
  hec_theta_beam	= new TH2F("hec_theta_beam","hec_theta_beam",theta_bins,theta_binmin,theta_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_theta_beam);
  hec_theta_beam_pos = new TH2F("hec_theta_beam_pos","hec_theta_beam_pos",theta_bins,theta_binmin,theta_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_theta_beam_pos);
  hec_theta_beam_pos_kc = new TH2F("hec_theta_beam_pos_kc","hec_theta_beam_pos_kc",theta_bins,theta_binmin,theta_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_theta_beam_pos_kc);
  hec_theta_beam_pos_pad0 = new TH2F("hec_theta_beam_pos_pad0","hec_theta_beam_pos_pad0",theta_bins,theta_binmin,theta_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_theta_beam_pos_pad0);
  hec_theta_beam_pos_pad1 = new TH2F("hec_theta_beam_pos_pad1","hec_theta_beam_pos_pad1",theta_bins,theta_binmin,theta_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_theta_beam_pos_pad1);
  hec_theta_beam_pos_pad01 = new TH2F("hec_theta_beam_pos_pad01","hec_theta_beam_pos_pad01",theta_bins,theta_binmin,theta_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_theta_beam_pos_pad01);

  hec_theta_beam_pos_tpzero = new TH2F("hec_theta_beam_pos_tpzero","hec_theta_beam_pos_tpzero",theta_bins,theta_binmin,theta_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_theta_beam_pos_tpzero);
  hec_theta_beam_pos_tpone = new TH2F("hec_theta_beam_pos_tpone","hec_theta_beam_pos_tpone",theta_bins,theta_binmin,theta_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_theta_beam_pos_tpone);
  hec_theta_beam_pos_tpzeroone = new TH2F("hec_theta_beam_pos_tpzeroone","hec_theta_beam_pos_tpzeroone",theta_bins,theta_binmin,theta_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_theta_beam_pos_tpzeroone);
  hec_theta_beam_pos_tpmone = new TH2F("hec_theta_beam_pos_tpmone","hec_theta_beam_pos_tpmone",theta_bins,theta_binmin,theta_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_theta_beam_pos_tpmone);

  hec_theta_beam_pos_tpmone->Draw();

  hec_x_all = new TH2F("hec_x_all","hec_x_all",targetx_bins, targetx_binmin, targetx_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_x_all);	
  hec_x_beam = new TH2F("hec_x_beam","hec_x_beam",targetx_bins, targetx_binmin, targetx_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_x_beam);
  hec_x_beam_pos	= new TH2F("hec_x_beam_pos","hec_x_beam_pos",targetx_bins, targetx_binmin, targetx_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_x_beam_pos);
  hec_y_all = new TH2F("hec_y_all","hec_y_all",targety_bins, targety_binmin, targety_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_y_all);	
  hec_y_beam = new TH2F("hec_y_beam","hec_y_beam",targety_bins, targety_binmin, targety_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_y_beam);
  hec_y_beam_pos	= new TH2F("hec_y_beam_pos","hec_y_beam_pos",targety_bins, targety_binmin, targety_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_y_beam_pos);	
  hec_anglea_all = new TH2F("hec_anglea_all","hec_anglea_all",anglea_bins,anglea_binmin,anglea_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_anglea_all);
  hec_anglea_beam = new TH2F("hec_anglea_beam","hec_anglea_beam",anglea_bins,anglea_binmin,anglea_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_anglea_beam);
  hec_anglea_beam_pos = new TH2F("hec_anglea_beam_pos","hec_anglea_beam_pos",anglea_bins,anglea_binmin,anglea_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_anglea_beam_pos);
  hec_angleb_all = new TH2F("hec_angleb_all","hec_angleb_all",angleb_bins,angleb_binmin,angleb_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_angleb_all);
  hec_angleb_beam = new TH2F("hec_angleb_beam","hec_angleb_beam",angleb_bins,angleb_binmin,angleb_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_angleb_beam);
  hec_angleb_beam_pos = new TH2F("hec_angleb_beam_pos","hec_angleb_beam_pos",angleb_bins,angleb_binmin,angleb_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_angleb_beam_pos);
  
  hec_totreco_all = new TH2F("hec_totreco_all","hec_totreco_all",200,0,40,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_totreco_all);
  hec_totreco_beam = new TH2F("hec_totreco_beam","hec_totreco_beam",200,0,40,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_totreco_beam);
  hec_totreco_beam_pos = new TH2F("hec_totreco_beam_pos","hec_totreco_beam_pos",200,0,40,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_totreco_beam_pos);
  hec_excen_all = new TH2F("hec_excen_all","hec_excen_all",excen_bins, excen_binmin, excen_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_excen_all);	
  hec_excen_beam	= new TH2F("hec_excen_beam","hec_excen_beam",excen_bins, excen_binmin, excen_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_excen_beam);
  hec_excen_beam_pos = new TH2F("hec_excen_beam_pos","hec_excen_beam_pos",excen_bins, excen_binmin, excen_binmax,excen_bins, excen_binmin, excen_binmax);hlist->Add(hec_excen_beam_pos);
 

  
  TH2F* tdccorr = new TH2F("tdccorr","tdccorr",12,0,12,12,0,12);
  hlist->Add(tdccorr);

  int adc[NADC][32];
  int tdc[64];
  int siring[NDET];
  double sitheta[NDET];
  double sifsen[NDET];
  double sicorr[NDET];
  double sibsen[NDET];
  double csiSen[NDET];
#ifndef KYUSHU
  Int_t beampid;
  double s1paden[2];
  double csiLen[NDET];
  double theta[2][NDET];
  double excencorr[2][NDET]; // 0 is position, 1 is position and angle on tharget
#endif
  double totreco[2][NDET];
  double excen[NDET];
  //Double_t excencorr[NDET];
  Int_t tinapid[NDET];

  double targetx;
  double targety;
  double targeta;
  double targetb;



  tr->SetBranchAddress("adc",&adc);
  tr->SetBranchAddress("tdc",&tdc);
  tr->SetBranchAddress("siring",&siring);
  tr->SetBranchAddress("sitheta",&sitheta);
  tr->SetBranchAddress("sifsen",&sifsen);
  tr->SetBranchAddress("sicorr",&sicorr);
  tr->SetBranchAddress("sibsen",&sibsen);
  tr->SetBranchAddress("theta",&theta);
  tr->SetBranchAddress("totreco",&totreco);
  tr->SetBranchAddress("excen",&excen);
  tr->SetBranchAddress("excencorr",&excencorr);
  tr->SetBranchAddress("beampid",&beampid);
  tr->SetBranchAddress("tinapid",&tinapid);
  tr->SetBranchAddress("csiSen",&csiSen);
#ifndef KYUSHU
  tr->SetBranchAddress("csiLen",&csiLen);
  tr->SetBranchAddress("targetx",&targetx);
  tr->SetBranchAddress("targety",&targety);
  tr->SetBranchAddress("targeta",&targeta);
  tr->SetBranchAddress("targetb",&targetb);
  tr->SetBranchAddress("s1paden",&s1paden);
  //tr->SetBranchAddress("s1ppactime",&s1ppactime);
#endif

  Nucleus *prot = new Nucleus(1,0);

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
      sicorr[d] = NAN;
      theta[0][d] = NAN;
      theta[1][d] = NAN;
      totreco[0][d] = NAN;
      totreco[1][d] = NAN;
      tinapid[d]=-1;
      beampid=0;
      excen[d] = NAN;
      excencorr[0][d] = NAN;
      excencorr[1][d] = NAN;
      sibsen[d] = NAN;
      csiSen[d] = NAN;
#ifndef KYUSHU
      csiLen[d] = NAN;
      targetx = NAN;
      targety = NAN;
      targeta = NAN;
      targetb = NAN;
      s1paden[0]=NAN;
      s1paden[1]=NAN;
      //s1ppactime[0]=NAN;
      //s1ppactime[1]=NAN;
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

    Bool_t isontarget=false;
#ifndef KYUSHU

    htarget_xy_all->Fill(targetx, targety);
    if(beampid){
      htarget_xy_beam->Fill(targetx, targety);
    }
    if(InsideEllipse(bcx, bcy, brx, bry, 0.0, targetx, targety)){
      isontarget=true;
      htarget_xy_poscut->Fill(targetx, targety);
      if(beampid){
        htarget_xy_beam_pos->Fill(targetx, targety);
	if(s1paden[0]>0){
	  htarget_xy_beam_pos_pad0->Fill(targetx, targety);
	}
      }
    }
#endif


    for(unsigned short d=0;d<NDET;d++){
      //cout << sifsen[d] << "\t "<< sibsen[d]<< "\t "<< csiSen[d] << endl;
      if(sifsen[d]>0)
	hsifsen_sitheta[d]->Fill(sitheta[d],sifsen[d]);
      if(sifsen[d]>0 && csiSen[d]>0){
	hsifsen_csiS[d]->Fill(csiSen[d],sifsen[d]);
	//double alpha;
//#ifdef KYUSHU
//	//alpha = 90.-detectorangle-sitheta[d];
//#else
//	//alpha = -90.-detectorangle+sitheta[d];
//#endif
	//double ecorr = fabs(cos(alpha*deg2rad))*sifsen[d];
	hsicorr_csiS[d]->Fill(csiSen[d],sicorr[d]);
	hsicorr_csi_all->Fill(csiSen[d],sicorr[d]);
	hsifsen_csi_all->Fill(csiSen[d],sifsen[d]);
	//if(csiSCut[d] && csiSCut[d]->IsInside(csiSen[d],ecorr)){
	if(csiSCut[d] && csiSCut[d]->IsInside(csiSen[d],sifsen[d])){
	  //cout << d <<",\tE = "<< sifsen[d] << ",\tth = "<< sitheta[d]<< ",\tcsiraw = "<< csien[d] << endl;
	  //if(siring[d]>-1){
	  //  double en = protDete_l2e[d][siring[d]]->Eval(sifsen[d]);
	  //  //cout << siring[d] << "\t" << sifsen[d] << "\t" << en << endl;
	  //  hcal_csiS[d]->Fill(csiSen[d],en-sifsen[d]);
	  //}
	}
      }//si and csi en >0
      if(sibsen[d]>0 && csiSen[d]>0)
	hsibsen_csiS[d]->Fill(csiSen[d],sibsen[d]);
#ifndef KYUSHU
      if(sifsen[d]>0 && csiLen[d]>0){
	hsifsen_csiL[d]->Fill(csiLen[d],sifsen[d]);
	hsifsen_csi_all->Fill(csiLen[d],sifsen[d]);
	//double alpha;
	//alpha = -90.-detectorangle+sitheta[d];
	//double ecorr = fabs(cos(alpha*deg2rad))*sifsen[d];
	hsicorr_csiL[d]->Fill(csiLen[d],sicorr[d]);
	hsicorr_csi_all->Fill(csiLen[d],sicorr[d]);
	//if(csiLCut[d] && csiLCut[d]->IsInside(csiLen[d],ecorr)){
	//if(csiLCut[d] && csiLCut[d]->IsInside(csiLen[d],sifsen[d])){
	//  //cout << d <<",\tE = "<< sifsen[d] << ",\tth = "<< sitheta[d]<< ",\tcsiraw = "<< csien[d] << endl;
	//  //if(siring[d]>-1){
	//  //  double en = protDete_l2e[d][siring[d]]->Eval(sifsen[d]);
	//  //  //cout << siring[d] << "\t" << sifsen[d] << "\t" << en << endl;
	//  //  hcal_csiL[d]->Fill(csiLen[d],en-sifsen[d]);
	//  //}
	//}
      }//si and csi en >0
      if(sibsen[d]>0 && csiLen[d]>0)
	hsibsen_csiL[d]->Fill(csiLen[d],sibsen[d]);
#endif    
    
      
      if(totreco[1][d]>0){
	htr_sitheta[d]->Fill(sitheta[d],totreco[1][d]);    
	htr_sitheta_all->Fill(sitheta[d],totreco[1][d]);    
	
	htr_theta_all->Fill(theta[0][d],totreco[1][d]);
	
	hexcen_theta_all->Fill(theta[0][d],excen[d]);    
	hexcen_totreco_all->Fill(totreco[1][d],excen[d]);    
	
	hec_all->Fill(excencorr[0][d]);	
	hec_theta_all->Fill(theta[0][d],excencorr[0][d]);	
	hec_x_all->Fill(targetx,excencorr[0][d]);
	hec_y_all->Fill(targety,excencorr[0][d]);	
	hec_anglea_all->Fill(targeta,excencorr[0][d]);	
	hec_angleb_all->Fill(targetb,excencorr[0][d]);	
	
	hec_totreco_all->Fill(totreco[1][d],excencorr[0][d]);    
	hec_excen_all->Fill(excen[d],excencorr[0][d]);	
	
	if(beampid){
	  htr_theta_beam->Fill(theta[0][d],totreco[1][d]);
	  hexcen_theta_beam->Fill(theta[0][d],excen[d]);    
          hec_beam->Fill(excencorr[0][d]);
          hec_theta_beam->Fill(theta[0][d],excencorr[0][d]);
	  hec_x_beam->Fill(targetx,excencorr[0][d]);
	  hec_y_beam->Fill(targety,excencorr[0][d]);
	  hec_anglea_beam->Fill(targeta,excencorr[0][d]);	
	  hec_angleb_beam->Fill(targetb,excencorr[0][d]);	
	  
	  hec_totreco_beam->Fill(totreco[1][d],excencorr[0][d]);    
	  hec_excen_beam->Fill(excen[d],excencorr[0][d]);

	  //if(InsideEllipse(bcx, bcy, brx, bry, 0.0, targetx, targety)){
	  if(isontarget){
	    htr_theta_beam_pos->Fill(theta[0][d],totreco[1][d]);
	    hexcen_theta_beam_pos->Fill(theta[0][d],excen[d]);
	    //hexcen_corr->Fill(excencorr[d]);
	    hec_beam_pos->Fill(excencorr[0][d]);		
	    hec_theta_beam_pos->Fill(theta[0][d],excencorr[0][d]);		

	    if(kineCut->IsInside(theta[0][d], totreco[1][d])){
	      htr_theta_beam_pos_kc->Fill(theta[0][d],totreco[1][d]);
	      hec_beam_pos_kc->Fill(excencorr[0][d]);		
	      hec_theta_beam_pos_kc->Fill(theta[0][d],excencorr[0][d]);		
	    }


	    if(tinapid[d]==0){
	      hec_beam_pos_tpzero->Fill(excencorr[0][d]);
	      hec_theta_beam_pos_tpzero->Fill(theta[0][d],excencorr[0][d]);
	    }
	    if(tinapid[d]==1){
	      hec_beam_pos_tpone->Fill(excencorr[0][d]);
	      hec_theta_beam_pos_tpone->Fill(theta[0][d],excencorr[0][d]);
	    }
	    if(tinapid[d]==0 || tinapid[d]==1){
	      hec_beam_pos_tpzeroone->Fill(excencorr[0][d]);
	      hec_theta_beam_pos_tpzeroone->Fill(theta[0][d],excencorr[0][d]);
	    }
	    if(tinapid[d]==-1){
	      hec_beam_pos_tpmone->Fill(excencorr[0][d]);
	      hec_theta_beam_pos_tpmone->Fill(theta[0][d],excencorr[0][d]);
	    }
	    hec_x_beam_pos->Fill(targetx,excencorr[0][d]);		
	    hec_y_beam_pos->Fill(targety,excencorr[0][d]);		
	    hec_anglea_beam_pos->Fill(targeta,excencorr[0][d]);		
	    hec_angleb_beam_pos->Fill(targetb,excencorr[0][d]);		
	    
	    hec_totreco_beam_pos->Fill(totreco[1][d],excencorr[0][d]);    
	    hec_excen_beam_pos->Fill(excen[d],excencorr[0][d]);		
	    //printf("filling hec_excen_beam_pos: %f %f\n", excen[d],excencorr[0][d]);

	    // fill with cuts on s1 pads
	    // s1ic_lowercut
	    if(s1paden[0]>0){
	            hec_beam_pos_pad0->Fill(excencorr[0][d]);
	      hec_theta_beam_pos_pad0->Fill(theta[0][d],excencorr[0][d]);
	      htr_theta_beam_pos_pad0->Fill(theta[0][d],totreco[1][d]);
	       hexcen_theta_beam_pos_pad0->Fill(theta[0][d],excen[d]);
	       if(csiLen[d]>0 && sifsen[d]>0){
	         hsicorr_csi_beam_pos_pad0->Fill(csiLen[d], sicorr[d]);
	       }
	       if(csiSen[d]>0 && sifsen[d]>0){
	         hsicorr_csi_beam_pos_pad0->Fill(csiSen[d], sicorr[d]);
	       }
	    }
	    if(s1paden[1]>0){
	            hec_beam_pos_pad1->Fill(excencorr[0][d]);
	      hec_theta_beam_pos_pad1->Fill(theta[0][d],excencorr[0][d]);
	      htr_theta_beam_pos_pad1->Fill(theta[0][d],totreco[1][d]);
	       hexcen_theta_beam_pos_pad1->Fill(theta[0][d],excen[d]);
	    }
	    double s1padenmean = TMath::Sqrt(s1paden[0]*s1paden[1]);
	    if(s1padenmean>0){
	            hec_beam_pos_pad01->Fill(excencorr[0][d]);
	      hec_theta_beam_pos_pad01->Fill(theta[0][d],excencorr[0][d]);
	      htr_theta_beam_pos_pad01->Fill(theta[0][d],totreco[1][d]);
	       hexcen_theta_beam_pos_pad01->Fill(theta[0][d],excen[d]);
	    }
	    //if((s1ppactime[0]-tdc[0]>-14200) && (s1ppactime[0]-tdc[0]<-13500)){
	    //  htr_theta_beam_pos_s1t0->Fill(theta[0][d],totreco[1][d]);
	    //}
	    //if((s1ppactime[0]-tdc[0]>-12700) && (s1ppactime[0]-tdc[0]<-11950)){
	    //  htr_theta_beam_pos_s1t1->Fill(theta[0][d],totreco[1][d]);
	    //}

	  }
	}
      }


//#endif
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
  //for(unsigned short d=0;d<NDET;d++){
  //  protDete_l2e[d][0]->Write(Form("ploss2energy%d_%d",d,0),TObject::kOverwrite);
  //  deutDete_l2e[d][0]->Write(Form("dloss2energy%d_%d",d,0),TObject::kOverwrite);
  //  protDete_l2e[d][15]->Write(Form("ploss2energy%d_%d",d,15),TObject::kOverwrite);
  //  deutDete_l2e[d][15]->Write(Form("dloss2energy%d_%d",d,15),TObject::kOverwrite);
  //}
  
  vector<TSpline3*> dp;
  TSpline3 *dd;
  TSpline3 *pp;
  for(unsigned int i =0;i< dpkine.size();i++){
    dp.push_back(dpkine.at(i)->EvslabMeV(0,180,1,2));
    dp.back()->SetName(Form("dp%d",i));
    dp.back()->Write();
  }
  dd = ddkine->EvslabMeV(0,90,1,2);
  dd->SetName("dd");
  dd->Write();
  pp = ppkine->EvslabMeV(0,90,1,2);
  pp->SetName("pp");
  pp->Write();


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
    for(int i=0;i<NSTRIPS;i++){
      m[i] = i;
    }
    return m;
  }
  infile.ignore(1000,'\n');
  for(int i=0;i<NSTRIPS;i++){
    int ch,st;
    double th, dth, dist;
    infile >> ch >> st >> th >> dth >> dist;
    infile.ignore(1000,'\n');
    m[st] = th;
    stripTheta[st]=th;
    stripDist[st]=dist;
  }
  return m;
}

void calckinematics(double midtarget, char* settingsfile){
  TEnv* settings = new TEnv(settingsfile);
  int Abeam = settings->GetValue("Beam.A",99);
  int Zbeam = settings->GetValue("Beam.Z",12);

  Nucleus *prot = new Nucleus(1,0);
  Nucleus *deut = new Nucleus(1,1);
  Nucleus *proj = new Nucleus(Zbeam,Abeam-Zbeam);
  Nucleus *ejec = new Nucleus(Zbeam,Abeam-Zbeam+1);

  ddkine = new Kinematics(proj,deut,deut,proj,midtarget,0);
  ppkine = new Kinematics(proj,prot,prot,proj,midtarget,0);
  for(int i=0;i<6;i++)
    dpkine.push_back(new Kinematics(proj,deut,prot,ejec,midtarget,kine_energies[i]));

  cout << "calculated kinematics " << endl;
}
//void calcenergyloss(char* detectorsetup){
//  double detectorthick[NDET];
//  TEnv* detsetup = new TEnv(detectorsetup);
//  double detectorangle = detsetup->GetValue("Detector.Angle",50.0);
//  for(int i=0;i<NDET;i++){
//    detectorthick[i] = detsetup->GetValue(Form("Detector.%d",i),300.0);
//  }
//  string channelmap = detsetup->GetValue("Channel.Mapping", "channelmapping.dat");
//  //initalize mapping from ADC channel to strip number and angle
//  map<int, double> strmap = strip2thetamap((char*)channelmap.c_str());
//  for(int i=0;i<16;i++){
//    cout << i << "\t" << strmap[i] << "\t" << endl;
//  }
//  Nucleus *si = new Nucleus(14,14);
//  Compound *dete = new Compound(si);
//  
//  Nucleus *prot = new Nucleus(1,0);
//  Nucleus *deut = new Nucleus(1,1);
//
//  Reconstruction *protDete = new Reconstruction(prot, dete);
//  Reconstruction *deutDete = new Reconstruction(deut, dete);
//  double protDete_range;
//  double deutDete_range;
//  double detedensity = 2.32; //g/cm^3 
//  for(int d=0;d<NDET;d++){
//    detectorthick[d] *= detedensity*0.1;
//    cout << "calculating energy loss of " <<  prot->GetSymbol() << ", " << deut->GetSymbol() << " in detector " << d << " d = " << detectorthick[d] << " mg/cm^2" << endl;
//    for(int s=0;s<16;s++){
//      double alpha;
//#ifdef KYUSHU
//	alpha = 90.-detectorangle-strmap[s];
//#else
//	alpha = -90.-detectorangle+strmap[s];
//#endif
//      if(DEBUG)
//	cout << s << "\t" << alpha << "\t" << cos(alpha*deg2rad) << endl;
//      protDete->SetTargetThickness(detectorthick[d]/fabs(cos(alpha*deg2rad)));
//      protDete_l2e[d][s] = protDete->EnergyLoss2Energy(50,0.1);
//      deutDete->SetTargetThickness(detectorthick[d]/fabs(cos(alpha*deg2rad)));
//      deutDete_l2e[d][s] = deutDete->EnergyLoss2Energy(50,0.1);
//    }
//  }
//
//}
void readpidcuts(char* filename){
  for(int i=0;i<NDET;i++){
    deutCut[i] = NULL;
    protCut[i] = NULL;
    kineCut = NULL;
    csiSCut[i] = NULL;
#ifndef KYUSHU
    csiLCut[i] = NULL;
#endif
  }
  TFile* fc = new TFile(filename);
  for(int i=0;i<NDET;i++){
    deutCut[i] = (TCutG*)fc->Get(Form("deut%d",i));
    protCut[i] = (TCutG*)fc->Get(Form("prot%d",i));
    kineCut = (TCutG*)fc->Get(Form("kineCut"));
    csiSCut[i] = (TCutG*)fc->Get(Form("csiS%d",i));
#ifndef KYUSHU
    csiLCut[i] = (TCutG*)fc->Get(Form("csiL%d",i));
#endif
  }
}



Bool_t InsideEllipse(Double_t x_c, Double_t y_c, Double_t R_x, Double_t R_y, Double_t a_r, Double_t data_x, Double_t data_y){
  
  //exclude nan values from input data
  if(TMath::IsNaN(data_x) || TMath::IsNaN(data_y)){
    return false;
  }


  // parameter which indicates that event is inside:
  // if it is smaller than 1.0 the data point is inside the ellipse
  Double_t d_check=9999.;

  Double_t a_r_rad        = a_r*TMath::Pi()/180.;
  Double_t x_p    = data_x-x_c;
  Double_t y_p    = data_y-y_c;
  Double_t x_ir   = x_p*TMath::Cos(-a_r_rad)-y_p*TMath::Sin(-a_r_rad);
  Double_t y_ir   = x_p*TMath::Sin(-a_r_rad)+y_p*TMath::Cos(-a_r_rad);
  Double_t delta_x        = x_ir/R_x;
  Double_t delta_y        = y_ir/R_y;
  
  //squared values!!
  d_check = delta_x*delta_x+delta_y*delta_y;

  if(d_check<=1.0){
    return true;
  }
  else{
    return false;
  }

}

