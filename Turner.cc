#include "TArtEventStore.hh"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCutG.h>
#include <TTree.h>
#include <TObject.h>
#include <TEnv.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
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

#ifndef rad2deg
#define rad2deg                       180./(TMath::Pi())
#endif
#ifndef deg2rad
#define deg2rad                       (TMath::Pi())/180.
#endif

#define ADC0 11  //ch 0-15 YY1 No1, ch 16-31 YY1 No2
#define ADC1 12  //ch 0-15 YY1 No3, ch 16-31 YY1 No4
#define ADC2 13  //ch 0-15 YY1 No5, ch 16-31 YY1 No5
#define ADC3 14  //ch 0-15 CsI,     ch 16-31 YY1 backsides
#define ADC4 15  //spare
#define TDC0 0
#define TIMEREF 0

#define NDET 6
#define NSTRIPS 16
#define NADC 5

#define BEAMDEVICE 60
#define TINADEVICE 0
#define S1DEVICE 11
#define BEAMFOCALPLANE 0
#define TINAFOCALPLANE 0
#define S1FOCALPLANE 21
#define TOFDET 13
#define PPACDET 15
#define S1DET 31
#define S1PADS 2 // up to 12

#define OVERFLOWBINS 4080

int vlevel =0;
TSpline3* protTarg_e2r;
TSpline3* protTarg_r2e;
TSpline3* deutTarg_e2r;
TSpline3* deutTarg_r2e;
TSpline3* protFoil_e2r;
TSpline3* protFoil_r2e;
TSpline3* deutFoil_e2r;
TSpline3* deutFoil_r2e;
TSpline3* protDete_e2r[NDET];
TSpline3* protDete_r2e[NDET];
TSpline3* deutDete_e2r[NDET];
TSpline3* deutDete_r2e[NDET];
Kinematics* ddkine;
Kinematics* ppkine;
vector<Kinematics*> dpkine;

TCutG* deutCut[NDET];
TCutG* protCut[NDET];
TCutG* kineCut;
TCutG* pidSCut[NDET];
TCutG* pidLCut[NDET];

//geometric address to module number
int geo2num(int geo);
//map adc ch to strip number and angle
map<int,pair<int,double> > ch2stripmap(char* filename);
Double_t stripDist[NSTRIPS];
Double_t stripTheta[NSTRIPS];

//energy losses
double calcenergyloss(char* filename);
//kinematics
void calckinematics(double midtarget, char* settingsfile);
//pid cuts
void readpidcuts(const char* filename);

//signal handling, enables ctrl-c to exit nicely
bool signal_received = false;
void signalhandler(int sig);
double get_time();
int main(int argc, char** argv){
  double time_start = get_time();  
  signal(SIGINT,signalhandler);
#ifdef KYUSHU
  cout << "Turner compiled for Kyushu data!" << endl;
#else
  cout << "Turner compiled for OEDO day0!" << endl;
#endif
  Int_t RunNumber = -1;
  char* settingsfile = NULL;
  char* prefix = "./ridf/phys";//"./data/run";
  vlevel = 0;
  int LastEvent =-1;
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-r", "run number", &RunNumber);
  interface->Add("-s", "settings, beam, detector, calibration", &settingsfile);
  interface->Add("-le", "last event to be read", &LastEvent);  
  interface->Add("-v", "verbose level", &vlevel);  
  interface->Add("-p", "prefix for data files", &prefix);  

  //check flags and complain if somethings is missing
  interface->CheckFlags(argc, argv);
  if(RunNumber<0){
    cout << "invalid runnumber: " << RunNumber << endl;
    return 1;
  }
  if(settingsfile==NULL){
    cout << "Error: no settings file given" << endl;
    return -99;
  }

  //input and output files 
  cout <<" analyzing run " << RunNumber << endl;
  char* inputfilename  = (char*)Form("%s%04d.ridf",prefix,RunNumber);
  char* outputfilename = (char*)Form("root/run%04d.root",RunNumber);
  TArtEventStore *estore = new TArtEventStore();
  estore->Open(inputfilename);
  TArtRawEventObject *rawevent = estore->GetRawEventObject();
  TFile *fout = new TFile(outputfilename,"RECREATE");

  //output tree
  TTree *tr = new TTree("tr","events");
  Int_t adc[NADC][32];
  Int_t tdc[64];
  Int_t tdcchmult[64];
  Int_t adcmult;
  Int_t tdcmult;

  //all rings
  Double_t sirien[NDET][16];
 
  Double_t sifsti[NDET];
  //highest per detector (for now)
  Double_t sifsen[NDET];
  Int_t siring[NDET];
  Double_t sitheta[NDET];
  //strip mult per detector
  Int_t simult[NDET];
  //energyloss correction incident angle alpha on detector (just for pid), from sifsen
  Double_t sicorr[NDET];
    
  //backside
  Double_t sibsen[NDET];
  Double_t sibsti[NDET];
  
  //csi
  Double_t csiSen[NDET];
  Double_t csiSti[NDET];
#ifndef KYUSHU 
  Double_t csiLen[NDET];
  Double_t csiLti[NDET];
#endif

  //total energy for sifsen and csi
  Double_t toten[NDET];
  //reconstructed energy at target center, from sifsen + csi
  Double_t totreco[2][NDET]; //0 only foils, 1 foils and target energy loss reconstructed
  //pid flag: 0 identified (proton), 1 E/theta pid proton, 2 E/theta pid deuteron
  Int_t tinapid[NDET]; //check still kyushu type, id only for protons
  //excitation energy
  Double_t theta[2][NDET]; // for excen correction
  Double_t excen[NDET];
  Double_t excencorr[2][NDET]; // 0 position on target, 1 position and angle

#ifndef KYUSHU 
  //tof
  Double_t f3tof;
  Double_t f5tof;
  Double_t f3f5tof;
  Int_t beampid;

  //ppac position
  Int_t fe12tref;
  Int_t fe12traw[2];
  Double_t fe12t[2];
  Double_t fe12x[2];
  Double_t fe12y[2];
  Double_t targetx;
  Double_t targety;
  Double_t targeta;
  Double_t targetb;

  Int_t s1tref;
  Int_t s1traw[2]; // anodes
  Double_t s1t[2]; 
  Double_t s1x[2];
  Double_t s1y[2];
  Double_t s1paden[S1PADS]; 
  
  Int_t registr;
#endif
  
  Double_t targetz; // defined in settings

  //raw
  tr->Branch("adc",&adc,Form("adc[%d][32]/I",NADC));
  tr->Branch("tdc",&tdc,"tdc[64]/I");
  tr->Branch("tdcchmult",&tdcchmult,"tdcchmult[64]/I");
  tr->Branch("adcmult",&adcmult,"adcmult/I");
  tr->Branch("tdcmult",&tdcmult,"tdcmult/I");
  //cal
  tr->Branch("sirien",  &sirien,  Form("sirien[%d][16]/D",NDET));
  tr->Branch("sifsti",  &sifsti,  Form("sifsti[%d]/D",NDET));
  tr->Branch("sifsen",  &sifsen,  Form("sifsen[%d]/D",NDET));
  tr->Branch("sitheta", &sitheta, Form("sitheta[%d]/D",NDET));
  tr->Branch("siring",  &siring,  Form("siring[%d]/I",NDET));
  tr->Branch("simult",  &simult,  Form("simult[%d]/I",NDET));

  tr->Branch("sicorr",  &sicorr,  Form("sicorr[%d]/D",NDET));
  
  tr->Branch("toten",   &toten,   Form("toten[%d]/D",NDET));
  tr->Branch("tinapid", &tinapid, Form("tinapid[%d]/I",NDET));
  tr->Branch("totreco", &totreco, Form("totreco[2][%d]/D",NDET));
  tr->Branch("theta", &theta, Form("theta[2][%d]/D",NDET));
  tr->Branch("excen",   &excen,   Form("excen[%d]/D",NDET));
  tr->Branch("excencorr",   &excencorr,   Form("excencorr[2][%d]/D",NDET));

  tr->Branch("sibsen",  &sibsen,  Form("sibsen[%d]/D",NDET));
  tr->Branch("sibsti",  &sibsti,  Form("sibsti[%d]/D",NDET));
  tr->Branch("csiSen",  &csiSen,  Form("csiSen[%d]/D",NDET));
  tr->Branch("csiSti",  &csiSti,  Form("csiSti[%d]/D",NDET));
#ifndef KYUSHU 
  tr->Branch("csiLen",  &csiLen,  Form("csiLen[%d]/D",NDET));
  tr->Branch("csiLti",  &csiLti,  Form("csiLti[%d]/D",NDET));

  //trigger register
  tr->Branch("registr", &registr, "registr/I");

  //tof
  tr->Branch("f3tof",   &f3tof,   "f3tof/D");
  tr->Branch("f5tof",   &f5tof,   "f5tof/D");
  tr->Branch("f3f5tof", &f3f5tof, "f3f5tof/D");
  tr->Branch("beampid", &beampid, "beampid/I");

  //ppac positions
  tr->Branch("fe12tref",   &fe12tref,   "fe12tref/I");
  tr->Branch("fe12traw",   &fe12traw,   "fe12traw[2]/I");
  tr->Branch("fe12t",   &fe12t,   "fe12t[2]/D");
  tr->Branch("fe12x",   &fe12x,   "fe12x[2]/D");
  tr->Branch("fe12y",   &fe12y,   "fe12y[2]/D");
  tr->Branch("targeta", &targeta, "targeta/D");
  tr->Branch("targetb", &targetb, "targetb/D");
  tr->Branch("targetx", &targetx, "targetx/D");
  tr->Branch("targety", &targety, "targety/D");
  
  //s1 stuff
  tr->Branch("s1tref",   &s1tref,   "s1tref/I");
  tr->Branch("s1traw",   &s1traw,   "s1traw[2]/I");
  tr->Branch("s1t",   &s1t,   "s1t[2]/D");
  tr->Branch("s1x",   &s1x,   "s1x[2]/D");
  tr->Branch("s1y",   &s1y,   "s1y[2]/D");
  //tr->Branch("s1padraw",   &s1padraw,   "s1padraw[2]/I");
  tr->Branch("s1paden",   &s1paden,   "s1paden[2]/D");
#endif


  TEnv* settings = new TEnv(settingsfile);
  const char* calfile = settings->GetValue("Calibration.File","");
  const char* pedfile = settings->GetValue("Pedestal.File","");
  const char* pidfile = settings->GetValue("PID.Cut.File","");
  targetz = settings->GetValue("Target.Pos.Z", 0.0);
  if(strlen(calfile)==0)
    cout << "Warning: no calibration parameters given!" << endl;
  if(strlen(pedfile)==0)
    cout << "Warning: no pedestal parameters given!" << endl;
  if(strlen(pidfile)==0)
    cout << "Warning: no PID cuts given!" << endl;


  //define gain and offset parameters for calibration, (NADC) ADC, 1 TDC
  //pedestals == software thresholds for (NADC) adcs
  //check: only 32 channels of tdc here, if we need more this needs to be adjusted
  double ped[NDET][32];
  double gain[NDET][32];
  double offs[NDET][32];
  TEnv *cal = new TEnv(calfile);
  TEnv *pedestals = new TEnv(pedfile);
  for(int m=0;m<NDET;m++){
    for(int c=0;c<32;c++){
      if(strlen(calfile)==0){
	gain[m][c] = 1.0;
	offs[m][c] = 0.0;
      }
      else{
	gain[m][c] = cal->GetValue(Form("Gain.Module%d.Ch%d",m,c),0.0);
	offs[m][c] = cal->GetValue(Form("Offset.Module%d.Ch%d",m,c),0.0);
      }
      if(m<NADC){
	if(strlen(pedfile)==0)
	  ped[m][c] = 0.0;  
	else
	  ped[m][c] = pedestals->GetValue(Form("Pedestal.Module%d.Ch%d",m,c),0.0);
      }
      //cout << "gain["<<m<<"]["<<c<<"] = " << gain[m][c] << endl;
    } 
    
  }

  //reading in detector setup
  double detectorthick[NDET];
  double foilthick = settings->GetValue("Foil.Thickness",36.0);
  double targetthick = settings->GetValue("Target.Thickness",2.0);
  double foildensity = settings->GetValue("Foil.Density",2.7);
  double targetdensity = settings->GetValue("Target.Density",4.5);
  double detectorangle = settings->GetValue("Detector.Angle",50.0);
  double detectorphi[NDET];
  for(int i=0;i<NDET;i++){
    detectorthick[i] = settings->GetValue(Form("Detector.%d",i),300.0);
    detectorphi[i] = settings->GetValue(Form("Detector.Phi.%d", i), 0.0);
  }

  string channelmap = settings->GetValue("Channel.Mapping", "channelmapping.dat");
  //initalize mapping from ADC channel to strip number and angle
  map<int,pair<int,double> > chmap = ch2stripmap((char*)channelmap.c_str());
  if(vlevel>0){
    for(int i=0;i<16;i++)
      cout << i << "\t" << chmap[i].first << "\t" << chmap[i].second << "\t" << stripTheta[i] << "\t" << stripDist[i] << endl;
  }

  //prepare the energy reconstruction
  double emid = calcenergyloss(settingsfile);

  //calculate kinematics
  calckinematics(emid,settingsfile);
  fout->cd();
  vector<TSpline3*> dp;
  for(unsigned int i =0;i< dpkine.size();i++){
    dp.push_back(dpkine.at(i)->EvslabMeV(0,180,1,2));
    dp.back()->SetName(Form("dp%d",i));
    dp.back()->Write();
  }
  TSpline3 *dd = ddkine->EvslabMeV(0,90,1,2);
  dd->SetName("dd");
  dd->Write();
  TSpline3 *pp = ppkine->EvslabMeV(0,90,1,2);
  pp->SetName("pp");
  pp->Write();
  Nucleus *prot = new Nucleus(1,0);


  //read in pid cuts
  readpidcuts(pidfile);
  fout->cd();
  for(int i=0;i<NDET;i++){
    if(deutCut[i]!=NULL)
      deutCut[i]->Write();
    if(protCut[i]!=NULL)
      protCut[i]->Write();
    if(kineCut!=NULL)
      kineCut->Write();
    if(pidSCut[i]!=NULL)
      pidSCut[i]->Write();
    if(pidLCut[i]!=NULL)
      pidLCut[i]->Write();
  }
#ifndef KYUSHU  
  //tof cuts
  int f3refcut[3][2];
  int f5refcut[3][2];
  double tofoffset;
  double tofpidcut[2];
  
  int triggercut[3][2];

  char* namecut[2] = {"Lower","Upper"};
  char* trigname[3] = {"TINA","F3DS","PPAC"};
  for(int i=0;i<2;i++){
    tofpidcut[i] = settings->GetValue(Form("F3F5.PID.%s.Cut",namecut[i]),0.0);
    for(int t=0;t<3;t++){
      triggercut[t][i] = settings->GetValue(Form("Trigger.%s.%s.Cut",trigname[t],namecut[i]),0);
      f3refcut[t][i] = settings->GetValue(Form("F3.Ref.%s.%s.Cut",trigname[t],namecut[i]),0);
      f5refcut[t][i] = settings->GetValue(Form("F5.Ref.%s.%s.Cut",trigname[t],namecut[i]),0);
      if(f3refcut[t][i]==0)
	f3refcut[t][i] = settings->GetValue(Form("F3.Ref.%s.Cut",namecut[i]),0);
      if(f5refcut[t][i]==0)
	f5refcut[t][i] = settings->GetValue(Form("F5.Ref.%s.Cut",namecut[i]),0);
    }
  }
  cout << " beam pid cut " << tofpidcut[0] << " to " << tofpidcut[1] << " ns " << endl;
  tofoffset = settings->GetValue("F3F5.TOF.Offset",0.0);

  //fe12 ppac parameters
  double fe12delayoffset[2][2];
  double fe12linecalib[2][2];
  double fe12ns2mm[2][2];
  double fe12ppacpos[2][3];
  //double targetz;
  char* cxy[2] = {"X","Y"};

  for(int p=0;p<2;p++){
    for(int xy=0;xy<2;xy++){
      fe12delayoffset[p][xy] = settings->GetValue(Form("FE12.PPAC%d.%s.DelayOffset",p+1,cxy[xy]),0.0);
      fe12linecalib[p][xy] = settings->GetValue(Form("FE12.PPAC%d.%s.LineCalib",p+1,cxy[xy]),0.0);
      fe12ns2mm[p][xy] = settings->GetValue(Form("FE12.PPAC%d.%s.ns2mm",p+1,cxy[xy]),0.0);
      fe12ppacpos[p][xy] = settings->GetValue(Form("FE12.PPAC%d.Pos.%s",p+1,cxy[xy]),0.0);
    }
    fe12ppacpos[p][2] = settings->GetValue(Form("FE12.PPAC%d.Pos.Z",p+1),0.0);
  }
  //targetz = settings->GetValue("Target.Pos.Z",0.0);

  // s1
  double s1delayoffset[2][2];
  double s1linecalib[2][2];
  double s1ns2mm[2][2];
  double s1ppacpos[2][3];

  for(int p=0;p<2;p++){
    for(int xy=0;xy<2;xy++){
      s1delayoffset[p][xy] = settings->GetValue(Form("S1.PPAC%d.%s.DelayOffset",p+1,cxy[xy]),0.0);
      s1linecalib[p][xy] = settings->GetValue(Form("S1.PPAC%d.%s.LineCalib",p+1,cxy[xy]),0.0);
      s1ns2mm[p][xy] = settings->GetValue(Form("S1.PPAC%d.%s.ns2mm",p+1,cxy[xy]),0.0);
      s1ppacpos[p][xy] = settings->GetValue(Form("S1.PPAC%d.Pos.%s",p+1,cxy[xy]),0.0);
    }
    s1ppacpos[p][2] = settings->GetValue(Form("S1.PPAC%d.Pos.Z",p+1),0.0);
  }

  double s1ic_lowercut[S1PADS], s1ic_uppercut[S1PADS];
  for(int p=0;p<S1PADS;p++){
    s1ic_lowercut[p] = settings->GetValue(Form("S1.IC.PAD%d.Lower.Cut",p), 0.0);
    s1ic_uppercut[p] = settings->GetValue(Form("S1.IC.PAD%d.Upper.Cut",p), 999999.0);
  }

  // time cuts
  double fe12tcut[1][2]; // [detector][lower/upper], use only one global cut so far
  double s1tcut[1][2];
  fe12tcut[0][0] = settings->GetValue("FE12.Time.Lower.Cut", -999999);
  fe12tcut[0][1] = settings->GetValue("FE12.Time.Upper.Cut",  999999);
  s1tcut[0][0] = settings->GetValue("S1.Time.Lower.Cut", -999999);
  s1tcut[0][1] = settings->GetValue("S1.Time.Upper.Cut",  999999);

#endif

  int neve = 0;  // event number
  TRandom3 *rand = new TRandom3();


  int tdcevt = 0;
  int adcevt[NADC] = {};
  int bugevt = 0;
  TH2F* hraw[NADC+2];
  for(int i=0;i<NADC;i++){
    hraw[i] = new TH2F(Form("adc_%d",i),Form("adc_%d",i),32,0,32,4096,0,4096);
  }
  hraw[NADC] = new TH2F("tdc_0","tdc_0",64,0,64,4096,0,32*4096);
  hraw[NADC+1] = new TH2F("tdcchmult_0","tdcchmult_0",64,0,64,16,0,16);

  // TH1F* reff3_tof = new TH1F("reff3_tof","reff3_tof",1000,64000,65000);
  // TH1F* reff5_tof = new TH1F("reff5_tof","reff5_tof",2000,00,168000);
#ifndef KYUSHU
  TH1F* reff3_tof[3];
  TH1F* reff5_tof[3];
  for(int t=0;t<3;t++){
    reff3_tof[t] = new TH1F(Form("reff3_tof_%s",trigname[t]),Form("reff3_tof_%s",trigname[t]),2000,00,200000);
    reff5_tof[t] = new TH1F(Form("reff5_tof_%s",trigname[t]),Form("reff5_tof_%s",trigname[t]),2000,00,200000);
  }
  TH1F* f3f5_tof_all = new TH1F("f3f5_tof_all","f3f5_tof_all",6000,-300,300);
  TH1F* f3f5_tof = new TH1F("f3f5_tof","f3f5_tof",6000,-300,300);
  
  TH1F* f3mult = new TH1F("f3mult","f3mult",100,0,100);
  TH1F* f5mult = new TH1F("f5mult","f5mult",100,0,100);

  TH1F* fe12_ppac1x = new TH1F("fe12_ppac1x","fe12_ppac1x",4000,-200,200);
  TH1F* fe12_ppac1y = new TH1F("fe12_ppac1y","fe12_ppac1y",4000,-200,200);
  TH1F* fe12_ppac2x = new TH1F("fe12_ppac2x","fe12_ppac2x",4000,-200,200);
  TH1F* fe12_ppac2y = new TH1F("fe12_ppac2y","fe12_ppac2y",4000,-200,200);

#endif

  while(estore->GetNextEvent()){ 
    if(signal_received){
      break;
    }
    if(vlevel >0){
      cout << "-------------------------------------next event segments = "<<rawevent->GetNumSeg()<< endl;
    }
    if(neve == LastEvent)
      break;
    if(neve%10000==0){
      double time_end = get_time();
      cout << neve << " events analyzed, " << (Float_t)neve/(time_end - time_start) << " events/s\r" << flush;
    }
    //clear

    //tina
    for(int i=0;i<32;i++){
      for(int j=0;j<NADC;j++)
	adc[j][i] = 0;
    }
    for(int i=0;i<64;i++){
      tdc[i] = 0;
      tdcchmult[i] = 0;
    }
    adcmult = 0;
    tdcmult = 0;
    for(int i=0;i<NDET;i++){
      for(int j=0;j<16;j++)
	sirien[i][j] = sqrt(-1);
      simult[i] = 0;
      siring[i] = -1;
      sifsti[i] = sqrt(-1);
      sifsen[i] = sqrt(-1);
      sitheta[i] = sqrt(-1);
      sicorr[i] = sqrt(-1);
      toten[i] = sqrt(-1);
      totreco[0][i] = sqrt(-1);
      totreco[1][i] = sqrt(-1);
      tinapid[i] = -1;
      theta[0][i] = sqrt(-1);
      theta[1][i] = sqrt(-1);
      excen[i] = sqrt(-1);
      excencorr[0][i] = sqrt(-1);
      excencorr[1][i] = sqrt(-1);
      sibsen[i] = sqrt(-1);
      sibsti[i] = sqrt(-1);
      csiSen[i] = sqrt(-1);
      csiSti[i] = sqrt(-1);
#ifndef KYUSHU 
      csiLen[i] = sqrt(-1);
      csiLti[i] = sqrt(-1);
#endif
    }
#ifndef KYUSHU 
    //trigger
    registr =0;

    //tof
    f3tof = sqrt(-1);
    f5tof = sqrt(-1);
    f3f5tof = sqrt(-1);
    beampid = 0;

    //ppac position
    fe12tref = 999999;
    fe12traw[0] = 999999;
    fe12traw[1] = 999999;
    fe12x[0] = sqrt(-1);;
    fe12y[0] = sqrt(-1);;
    fe12x[1] = sqrt(-1);;
    fe12y[1] = sqrt(-1);;
    targetx = sqrt(-1);;
    targety = sqrt(-1);;
    targeta = sqrt(-1);;
    targetb = sqrt(-1);;

    //s1 stuff
    s1tref=999999;
    s1traw[0]=999999;
    s1traw[1]=999999;
    s1x[0] = sqrt(-1);;
    s1y[0] = sqrt(-1);;
    s1x[1] = sqrt(-1);;
    s1y[1] = sqrt(-1);;
    s1paden[0]=sqrt(-1);
    s1paden[1]=sqrt(-1);
#endif

    //temp data, not written to tree
    bool hadTINA = false;
    bool hadadc[NADC] = {};
    bool hadtdc = false;
    bool hadbug = false;

    int reftime = 0;
    vector<int> f3times;
    vector<int> f5times;
    f3times.clear();
    f5times.clear();

    bool refleading = true;
    bool f3leading = true;
    bool f5leading = true;

    int fe12xtime[2][2] = {{0,0},{0,0}};
    int fe12ytime[2][2] = {{0,0},{0,0}};

    int s1xtime[2][2] = {{0,0},{0,0}};
    int s1ytime[2][2] = {{0,0},{0,0}};
    int s1padraw[S1PADS] = {0,0};

    //loop over the segments
    for(int i=0;i<rawevent->GetNumSeg();i++){
      TArtRawSegmentObject *seg = rawevent->GetSegment(i);
      if(vlevel>1){
	cout << "event " << neve << "\tsegment " << i << "\tAddress " << seg->GetAddress() << "\tDevice " << seg->GetDevice() << "\tFocalPlane " << seg->GetFP() << "\tDetector " << seg->GetDetector() << "\tModule "<< seg->GetModule() << "\tnum data = " << seg->GetNumData() << endl;
      }
#ifndef KYUSHU
      if(seg->GetDevice()==BEAMDEVICE && seg->GetFP()==BEAMFOCALPLANE){
	//for PID:
	//dev 60, fp 0, det 12, geo 1 (moco + mtdc)  (13 instead of 12 for v1290) 
	//ch 10 F3dia9Pad
	//ch 26 F5dia9Pad
	if(vlevel>1)
	  cout << "BEAM data! " << endl;
	if(seg->GetDetector()==TOFDET){
	  if(vlevel>1)
	    cout << "TOF detector" << endl;
	  for(int j=0;j<seg->GetNumData();j++){
	    TArtRawDataObject *d = seg->GetData(j);
	    //get geometric address, channel and value
	    int geo = d->GetGeo(); 
	    int chan = d->GetCh();
	    int val = d->GetVal(); 
	    if(vlevel>2)
	      cout << "geo: " << geo  << "\tchan: " << chan  << "\tval: " << val << endl;  
	    if(geo==0){
	      if(chan==0){
		//cout << "Reference = " << val << endl;
		if(refleading){
		  reftime = val;
		  refleading = false;
		}
	      }
	      if(chan==10){
		//cout << "F3diaPad = " << val << endl;
		if(f3leading){
		  f3times.push_back(val);
		  f3leading = false;
		}
		else
		  f3leading = true;
	      }
	      if(chan==26){
		//cout << "F5diaPad = " << val << endl;
		if(f5leading){
		  f5times.push_back(val);
		  f5leading = false;
		}
		else
		  f5leading = true;
	      }
	    }
<<<<<<< HEAD
	    if(DEBUG)
	      cout << "det: " << det << "\tchan: " << chan  << "\tstrip: " << siring[det] << "\ten: " << en << endl;  
	    //cout << "num: " << num << "\tdet: " << det  << "\tchan: " << chan  << "\tchan%16: " << chan%16  << "\tchmap[chan%16].first: " << chmap[chan%16].first << "\tstrip: " << siring[det] << "\ten: " << en << endl;  
	  }
	  else if(geo==ADC0){
	    if(DEBUG)
	      cout << "geo: " << geo  << "\tchan: " << chan  << "\tval: " << val << endl; 
	    if(chan>15 && chan<15+NDET+1)//check
	      sibsen[chan-16] = en;
	    else if(chan<NDET)
	      csien[chan] = en;
	    else{
	      //cout << "unknown channel" << endl;
	      //cout << "geo: " << geo  << "\tchan: " << chan  << "\tval: " << val << endl; 
=======
	  }//data
	}//tof detector
	if(seg->GetDetector()==PPACDET){
	  if(vlevel>1)
	    cout << "PPAC detector" << endl;
	  for(int j=0;j<seg->GetNumData();j++){
	    TArtRawDataObject *d = seg->GetData(j);
	    //get geometric address, channel and value
	    int geo = d->GetGeo(); 
	    int chan = d->GetCh();
	    int val = d->GetVal(); 
	    if(vlevel>2)
	      cout << "geo: " << geo  << "\tchan: " << chan  << "\tval: " << val << endl;  
	    //PPAC 1
	    if(chan==0 && fe12tref==999999)//time ref
              fe12tref=val;
	    if(chan==16 && fe12traw[0]==999999)//PPAC1 A
              fe12traw[0]=val;
 	    if(chan==17 && fe12xtime[0][0]==0)//PPAC1 X1
	      fe12xtime[0][0]=val;
 	    if(chan==19 && fe12xtime[0][1]==0)//PPAC1 X2
	      fe12xtime[0][1]=val;
 	    if(chan==18 && fe12ytime[0][0]==0)//PPAC1 Y1
	      fe12ytime[0][0]=val;
 	    if(chan==20 && fe12ytime[0][1]==0)//PPAC1 Y2
	      fe12ytime[0][1]=val;
	    //PPAC2
	    if(chan==21 && fe12traw[1]==999999)//PPAC2 A
              fe12traw[1]=val;
 	    if(chan==22 && fe12xtime[1][0]==0)//PPAC2 X1
	      fe12xtime[1][0]=val;
 	    if(chan==24 && fe12xtime[1][1]==0)//PPAC2 X2
	      fe12xtime[1][1]=val;
 	    if(chan==23 && fe12ytime[1][0]==0)//PPAC2 Y1
	      fe12ytime[1][0]=val;
 	    if(chan==25 && fe12ytime[1][1]==0)//PPAC2 Y2
	      fe12ytime[1][1]=val;
	  }	  
	}//ppac detector	
      }//beam data
      
      // s1 stuff
      if(seg->GetDevice()==S1DEVICE && seg->GetFP()==S1FOCALPLANE){
	if(seg->GetDetector()==S1DET){
	  if(vlevel>1)
	    cout << "S1 detector" << endl;
	  for(int j=0;j<seg->GetNumData();j++){
	    TArtRawDataObject *d = seg->GetData(j);
	    //get geometric address, channel and value
	    int geo = d->GetGeo(); 
	    int chan = d->GetCh();
	    int val = d->GetVal(); 
	    if(vlevel>2)
	      cout << "geo: " << geo  << "\tchan: " << chan  << "\tval: " << val << endl;  
	    if(geo==0){// PPACs
	      //PPAC 1
	      if(chan==0 && s1tref==999999)//time ref
                s1tref=val;
	      if(chan==6 && s1traw[0]==999999)//PPAC1 A
                s1traw[0]=val;
 	      if(chan==2 && s1xtime[0][0]==0)//PPAC1 X1
	        s1xtime[0][0]=val;
 	      if(chan==4 && s1xtime[0][1]==0)//PPAC1 X2
	        s1xtime[0][1]=val;
 	      if(chan==3 && s1ytime[0][0]==0)//PPAC1 Y1
	        s1ytime[0][0]=val;
 	      if(chan==5 && s1ytime[0][1]==0)//PPAC1 Y2
	        s1ytime[0][1]=val;
	      //PPAC2
	      if(chan==11 && s1traw[1]==999999)//PPAC2 A
                s1traw[1]=val;
 	      if(chan==7 && s1xtime[1][0]==0)//PPAC2 X1
	        s1xtime[1][0]=val;
 	      if(chan==9 && s1xtime[1][1]==0)//PPAC2 X2
	        s1xtime[1][1]=val;
 	      if(chan==8 && s1ytime[1][0]==0)//PPAC2 Y1
	        s1ytime[1][0]=val;
 	      if(chan==10 && s1ytime[1][1]==0)//PPAC2 Y2
	        s1ytime[1][1]=val;
	    } // geo0
	    ////PPACS 
 	    //if(geo==0 && chan==6){//PPAC1 A
	    //  s1ppactime[0]=val;
	    // }
 	    //if(geo==0 && chan==11){//PPAC2 A
	    //  s1ppactime[1]=val;
	    //}
	    if(geo==1){ //IC pads
 	      if(chan==30 && s1padraw[0]==0) //IC Pad 0
	        s1padraw[0]=val;
 	      if(chan==29 && s1padraw[1]==0)//IC Pad 1
	        s1padraw[1]=val;
 	      //if(chan==28 && s1padraw[2]==0)//IC Pad 2
	      //  s1padraw[2]=val;
 	      //if(chan==27 && s1padraw[3]==0)//IC Pad 3
	      //  s1padraw[3]=val;
 	      //if(chan==26 && s1padraw[4]==0)//IC Pad 4
	      //  s1padraw[4]=val;
 	      //if(chan==25 && s1padraw[5]==0)//IC Pad 5
	      //  s1padraw[5]=val;
 	      //if(chan==24 && s1padraw[6]==0)//IC Pad 6
	      //  s1padraw[6]=val;
 	      //if(chan==23 && s1padraw[7]==0)//IC Pad 7
	      //  s1padraw[7]=val;
 	      //if(chan==22 && s1padraw[8]==0)//IC Pad 8
	      //  s1padraw[8]=val;
 	      //if(chan==21 && s1padraw[9]==0)//IC Pad 9
	      //  s1padraw[9]=val;
 	      //if(chan==20 && s1padraw[10]==0)//IC Pad 10
	      //  s1padraw[10]=val;
 	      //if(chan==19 && s1padraw[11]==0)//IC Pad 11
	      //  s1padraw[11]=val;
>>>>>>> 6d3290167d355101c1169de99269e4a20f9d2022
	    }
	  }	  
	}// s1 detector	
      }//s1 data
#endif


      if(seg->GetDevice()==TINADEVICE && seg->GetFP()==TINAFOCALPLANE){
	if(vlevel>1)
	  cout << "TINA data! " << endl;
	for(int j=0;j<seg->GetNumData();j++){
	  TArtRawDataObject *d = seg->GetData(j);
	  //get geometric address, channel and value
	  int geo = d->GetGeo(); 
	  int chan = d->GetCh();
	  int val = d->GetVal(); 
	  if(vlevel>1)
	    cout << "geo: " << geo  << "\tchan: " << chan  << "\tval: " << val << endl;  
	  int num, det;
	  double en;
	  double ti;
	  switch(geo){
	  case ADC0:
	  case ADC1:
	  case ADC2:
#ifndef KYUSHU
	  case ADC3:
	  case ADC4:
#endif
	    num = geo2num(geo); 
	    hadadc[num] = true;
	    adc[num][chan] = val;
	    hraw[num]->Fill(chan,val);
	    //if(val>0 && val<)
	    adcmult++;
	    if(val < ped[num][chan])
	      continue;
	    if(val > OVERFLOWBINS)
	      continue;
	    en = gain[num][chan]*(val+rand->Uniform(0,1)) + offs[num][chan];
	    det = -1;
#ifdef KYUSHU 
	    //kyushu: adc1 and 2 for 4 YY1, adc0 for backsides and csi
	    if(geo == ADC1)
	      det = chan/16;
	    if(geo == ADC2)
	      det = chan/16+2;
#else
	    //OEDO: adc0,1 and 2 for 6 YY1, adc3 for backsides and csi
	    if(geo == ADC0)
	      det = chan/16;
	    if(geo == ADC1)
	      det = chan/16+2;
	    if(geo == ADC2)
	      det = chan/16+4;
#endif
	    if(det>-1&&en>0){
	      simult[det]++;
	      int ring = chmap[chan%16].first; // mapping to real strip number
	      sirien[det][ring] = en;
	      if(siring[det]<0 || en > sifsen[det]){
		siring[det] = ring;
		sitheta[det] = chmap[chan%16].second; // mapping to theta
		sifsen[det] = en;
		double alpha;
#ifdef KYUSHU
		alpha = 90.-detectorangle-sitheta[det];
#else
		alpha = -90.-detectorangle+sitheta[det];
#endif
		sicorr[det] = fabs(cos(alpha*deg2rad))*en;
	      }
	      if(vlevel>1)
		cout << "det: " << det  << "\tchan: " << chan  << "\tstrip: " << siring[det] << "\ten: " << en << endl;  
	    }
#ifdef KYUSHU 
	    else if(geo==ADC0){
	      if(vlevel>1)
		cout << "geo: " << geo  << "\tchan: " << chan  << "\tval: " << val << endl; 
	      if(chan>15 && chan<15+NDET+1)//check
		sibsen[chan-16] = en;
	      else if(chan<NDET)
		csiSen[chan] = en;
	    }
#else
	    else if(geo==ADC3){
	      if(chan<NDET)
		sibsen[chan] = en;
	      else if(chan>15 && chan<16+NDET)
		csiSen[chan-16] = en;
	      else if(chan>15+NDET && chan<16+2*NDET)
		csiLen[chan-16-NDET] = en;
	    }
	    else if(geo==ADC4){
	      // cout << "ADC4 hit, should be empty/spare" << endl;
	      // cout << "geo: " << geo  << "\tchan: " << chan  << "\tval: " << val << endl; 
	    }
#endif	  
	    hadTINA = true;
	    if(vlevel>1)
	      cout << neve <<"\tADC\t" << num << "\t" << chan << "\t" << val <<"\tadc[num][chan] = " << adc[num][chan]<< endl;
	    break;
	  case TDC0:
	    hadtdc = true;
	    num = geo2num(geo); 
	    //take first hit only, maybe add more later
	    tdcchmult[chan]++;//since we read leading and trailing edge, this should be a multiple of 2
	    if(tdc[chan]==0){
	      tdc[chan] = val;
	      hraw[NADC]->Fill(chan,val);
	    }
	    tdcmult++;
	    if(vlevel>1)
	      cout << neve <<"\tTDC\t"<< chan << "\t" << val <<"\ttdc[chan] = " << tdc[chan]<<"\ttdcchmult[chan] = " << tdcchmult[chan]<< endl;
	    hadTINA = true;
	    break;
	  case 99:
	    hadbug = true;
	    //if()
	    cout << neve <<"\tBUG\t"<< chan << "\t" << val << endl;
	  default:
	    cout << "unknown geo = " << geo << " in TINAdevice and focal plane" << endl; 
	    break;
	  }//geo
	}//data
      }//tina data
    }//segments
    // if(hadTINA){
    //   cout << "----------------------" << endl;
    // }
    if(hadtdc){
      for(int t=0;t<3;t++){
	if(tdc[TIMEREF]-tdc[1+t] > triggercut[t][0] && tdc[TIMEREF]-tdc[1+t] < triggercut[t][1])
	  registr |= 1U << t;
      }
      if(vlevel>1){
	cout << "registr = " << registr << " = ";
	for(int b=0;b<3;b++){
	  int bb = (registr >> b) & 1U;
	  cout << bb;
	}
	cout << endl;
      }
      for(int ch=0;ch<64;ch++){
	if(tdcchmult[ch]>0)
	  hraw[NADC+1]->Fill(ch,tdcchmult[ch]);
      }
      tdcevt++;
      for(int i=0;i<NDET;i++){
 
#ifdef KYUSHU 
	int chan = i+2;
	if(tdc[chan]>0)
	  sibsti[i] = gain[geo2num(TDC0)][chan]*(tdc[chan]-tdc[TIMEREF]+rand->Uniform(0,1)) + offs[geo2num(TDC0)][chan];
	chan = i+6;
	if(tdc[chan]>0)
	  csiSti[i] = gain[geo2num(TDC0)][chan]*(tdc[chan]-tdc[TIMEREF]+rand->Uniform(0,1)) + offs[geo2num(TDC0)][chan];	
#else
	int chan = i+4;
	if(tdc[chan]>0)
	  sifsti[i] = gain[geo2num(TDC0)][chan]*(tdc[chan]-tdc[TIMEREF]+rand->Uniform(0,1)) + offs[geo2num(TDC0)][chan];
	chan = i+10;
	if(tdc[chan]>0)
	  sibsti[i] = gain[geo2num(TDC0)][chan]*(tdc[chan]-tdc[TIMEREF]+rand->Uniform(0,1)) + offs[geo2num(TDC0)][chan];
	chan = i+16;
	if(tdc[chan]>0)
	  csiSti[i] = gain[geo2num(TDC0)][chan]*(tdc[chan]-tdc[TIMEREF]+rand->Uniform(0,1)) + offs[geo2num(TDC0)][chan];
	chan = i+22;
	if(tdc[chan]>0)
	  csiLti[i] = gain[geo2num(TDC0)][chan]*(tdc[chan]-tdc[TIMEREF]+rand->Uniform(0,1)) + offs[geo2num(TDC0)][chan];
#endif
      }
    }
    for(int i=0;i<NADC;i++){
      if(hadadc[i])
	adcevt[i]++;
    }
    if(hadbug)
      bugevt++;
 
#ifndef KYUSHU   
    //tof
    // cout << "time fo flights " << endl;
    // cout << f3times.size() <<"\t" << f5times.size() << endl;
    f3mult->Fill(f3times.size());
    f5mult->Fill(f5times.size());
    for(vector<int>::iterator f3 = f3times.begin(); f3 != f3times.end(); f3++){
      //cout << *f3 << endl;
      for(int b=0;b<3;b++){
	if( (registr >> b) & 1U)
	  reff3_tof[b]->Fill(reftime - *f3);
	if( (reftime - *f3) > f3refcut[b][0] && (reftime - *f3) < f3refcut[b][1] )
	  f3tof = (*f3+rand->Uniform(0,1))*100./4096;
      }
      for(vector<int>::iterator f5 = f5times.begin(); f5 != f5times.end(); f5++){
     	f3f5_tof_all->Fill((*f5+rand->Uniform(0,1) - (*f3+rand->Uniform(0,1)))*100./4096 + tofoffset);
      } 
    }
    for(vector<int>::iterator f5 = f5times.begin(); f5 != f5times.end(); f5++){
      //cout << *f5 << endl;
      for(int b=0;b<3;b++){
	if( (registr >> b) & 1U)
	  reff5_tof[b]->Fill(reftime - *f5);
	if( (reftime - *f5) > f5refcut[b][0] && (reftime - *f5) < f5refcut[b][1] )
	  f5tof = (*f5+rand->Uniform(0,1))*100./4096;
      }
    } 
    
    f3f5tof = f5tof - f3tof + tofoffset; 
    f3tof = (reftime+rand->Uniform(0,1))*100./4096 - f3tof;
    f5tof = (reftime+rand->Uniform(0,1))*100./4096 - f5tof;

    f3f5_tof->Fill(f3f5tof);
    if(f3f5tof > tofpidcut[0] && f3f5tof < tofpidcut[1])
      beampid = 1;

    //ppac
    if(vlevel>2){
      cout << "PPAC1 X " << fe12xtime[0][0] << "\t" << fe12xtime[0][1] << endl;
      cout << "PPAC1 Y " << fe12ytime[0][0] << "\t" << fe12ytime[0][1] << endl;
      cout << "PPAC2 X " << fe12xtime[1][0] << "\t" << fe12xtime[1][1] << endl;
      cout << "PPAC2 Y " << fe12ytime[1][0] << "\t" << fe12ytime[1][1] << endl;
    }
    /*
    if(fe12x[0][0]>0&&fe12x[0][1]>0)
      fe12_ppac1x->Fill((fe12x[0][0]+rand->Uniform(0,1)-(fe12x[0][1]+rand->Uniform(0,1)))*100./1024);
    if(fe12y[0][0]>0&&fe12y[0][1]>0)
      fe12_ppac1y->Fill((fe12y[0][0]+rand->Uniform(0,1)-(fe12y[0][1]+rand->Uniform(0,1)))*100./1024);
    if(fe12x[1][0]>0&&fe12x[1][1]>0)
      fe12_ppac2x->Fill((fe12x[1][0]+rand->Uniform(0,1)-(fe12x[1][1]+rand->Uniform(0,1)))*100./1024);
    if(fe12y[1][0]>0&&fe12y[1][1]>0)
      fe12_ppac2y->Fill((fe12y[1][0]+rand->Uniform(0,1)-(fe12y[1][1]+rand->Uniform(0,1)))*100./1024);
    */
    
    // FE12 calibration
    for(int p=0;p<2;p++){
      if(fe12xtime[p][0]>0&&fe12xtime[p][1]>0){
	fe12x[p] = (fe12xtime[p][0]+rand->Uniform(0,1)-(fe12xtime[p][1]+rand->Uniform(0,1)))*100./1024; //in ns
	fe12x[p] += fe12delayoffset[p][0] - fe12linecalib[p][0];
	fe12x[p] *= fe12ns2mm[p][0]*0.5;
	fe12x[p] -= fe12ppacpos[p][0];	
      }
      if(fe12ytime[p][0]>0&&fe12ytime[p][1]>0){
	fe12y[p] = (fe12ytime[p][0]+rand->Uniform(0,1)-(fe12ytime[p][1]+rand->Uniform(0,1)))*100./1024; //in ns
	fe12y[p] += fe12delayoffset[p][1] - fe12linecalib[p][1];
	fe12y[p] *= fe12ns2mm[p][1]*0.5;
	fe12y[p] -= fe12ppacpos[p][1]; 
      }
      if((fe12traw[p]-fe12tref)>fe12tcut[0][0] && (fe12traw[p]-fe12tref)<fe12tcut[0][1]){
        fe12t[p] = (fe12traw[p]-fe12tref)*100.0/1024.0;
      }
    }
    fe12_ppac1x->Fill(fe12x[0]);
    fe12_ppac1y->Fill(fe12y[0]);
    fe12_ppac2x->Fill(fe12x[1]);
    fe12_ppac2y->Fill(fe12y[1]);

    //extrapolate to target
    if(!isnan(fe12x[0]) && !isnan(fe12x[1])){
      targeta = atan((fe12x[1] - fe12x[0])/(fe12ppacpos[1][2]-fe12ppacpos[0][2]))*1000.; //in mrad
      targetx = fe12x[1] + (targetz - fe12ppacpos[1][2])*tan(targeta/1000);
    }
    if(!isnan(fe12y[0]) && !isnan(fe12y[1])){
      targetb = atan((fe12y[1] - fe12y[0])/(fe12ppacpos[1][2]-fe12ppacpos[0][2]))*1000.; //in mrad
      targety = fe12y[1] + (targetz - fe12ppacpos[1][2])*tan(targetb/1000);
    }
    
    
    // S1 calibration
    for(int p=0;p<2;p++){
      if(s1xtime[p][0]>0&&s1xtime[p][1]>0){
	s1x[p] = (s1xtime[p][0]+rand->Uniform(0,1)-(s1xtime[p][1]+rand->Uniform(0,1)))*100./1024; //in ns
	s1x[p] += s1delayoffset[p][0] - s1linecalib[p][0];
	s1x[p] *= s1ns2mm[p][0]*0.5;
	s1x[p] -= s1ppacpos[p][0];	
      }
      if(s1ytime[p][0]>0&&s1ytime[p][1]>0){
	s1y[p] = (s1ytime[p][0]+rand->Uniform(0,1)-(s1ytime[p][1]+rand->Uniform(0,1)))*100./1024; //in ns
	s1y[p] += s1delayoffset[p][1] - s1linecalib[p][1];
	s1y[p] *= s1ns2mm[p][1]*0.5;
	s1y[p] -= s1ppacpos[p][1]; 
      }
      if(s1padraw[p]>0 && s1ic_lowercut[p]<s1padraw[p] && s1ic_uppercut[p]>s1padraw[p])
        s1paden[p] = (Double_t)s1padraw[p]; // put here calibration parameter, if needed
      if((s1traw[p]-s1tref)>s1tcut[0][0] && (s1traw[p]-s1tref)<s1tcut[0][1]){
        s1t[p] = (s1traw[p]-s1tref)*100.0/1024.0;
      }
    }
#endif

    //reconstruc energy loss
    if(hadTINA){
      //calculate total energy
      for(int i=0;i<NDET;i++){
<<<<<<< HEAD
	if(pidCut[i]!=NULL && pidCut[i]->IsInside(csien[i],sicorr[i])){
	  pid[i] = 10;
=======
#ifndef KYUSHU 
	if(csiLen[i] > 0 && csiSen[i] >0){
	  if(csiLen[i]>csiSen[i])
	    toten[i] = sifsen[i] + csiLen[i];
	  else
	    toten[i] = sifsen[i] + csiSen[i];
>>>>>>> 6d3290167d355101c1169de99269e4a20f9d2022
	}
	else if(csiLen[i] > 0){
	  toten[i] = sifsen[i] + csiLen[i];
	}
	else if(csiSen[i] > 0)
	  toten[i] = sifsen[i] + csiSen[i];
	else
	  toten[i] = sifsen[i];
#else
	if(csiSen[i] > 0)
	  toten[i] = sifsen[i] + csiSen[i];
	else
	  toten[i] = sifsen[i];
#endif
      }
      
      for(int i=0;i<NDET;i++){
        //if(sicorr[i]>0)
	//printf("i = %d, sitheta[i] = %f, sifsen[i] = %f, csiLen[i] = %f\n", i, sitheta[i], sifsen[i], csiLen[i]);
#ifdef KYUSHU
	if(csiSen[i] > 0 && pidSCut[i]!=NULL && pidSCut[i]->IsInside(csiSen[i],sicorr[i]))
	  tinapid[i] = 0;
#else
 	if(csiLen[i] > 0 && pidLCut[i]!=NULL && pidLCut[i]->IsInside(csiLen[i],sicorr[i]))
	  tinapid[i] = 0;
 	else if(csiSen[i] > 0 && pidSCut[i]!=NULL && pidSCut[i]->IsInside(csiSen[i],sicorr[i]))
	  tinapid[i] = 0;
	
	// if only Si fire, assume proton
	else if(sicorr[i]>0 && !(pidLCut[i]->IsInside(csiLen[i],sicorr[i])) && (TMath::IsNaN(csiLen[i]) || csiLen[i]<3.0) ){
	  tinapid[i] = 1;
	  toten[i]=sifsen[i];
	}
#endif
	else if(deutCut[i]!=NULL && deutCut[i]->IsInside(sitheta[i],sifsen[i]))
	  tinapid[i] = 2;
	else if(protCut[i]!=NULL && protCut[i]->IsInside(sitheta[i],sifsen[i]))
	  tinapid[i] = 1;
	// if only Si fire, assume proton
	//else if(!protCut[i]->IsInside(sitheta[i],sifsen[i]) && !deutCut[i]->IsInside(sitheta[i],sifsen[i]) && (TMath::IsNaN(csiSen[i]) || csiSen[i]<2.0) )
	//  tinapid[i] = 1;
	else if(sicorr[i]>0 && !(pidSCut[i]->IsInside(csiSen[i],sicorr[i])) && (TMath::IsNaN(csiSen[i]) || csiSen[i]<2.0) ){
	  tinapid[i] = 1;
	  toten[i]=sifsen[i];
	}

	double ene = toten[i];
	double range;
	double alpha;
	if(tinapid[i]==0||tinapid[i]==1){
#ifdef KYUHSU
	  alpha = 90.-detectorangle-sitheta[i];
	  double althick = foilthick*foildensity*0.1/fabs(cos(alpha*deg2rad));
	  range = protFoil_e2r->Eval(ene);
	  ene = protFoil_r2e->Eval(range+althick);
	  totreco[0][i] = ene;
#else
	  alpha = -90.-detectorangle+sitheta[i];
#endif
	  //reconstruct energy loss in half the target
	  double tathick = targetthick/2*targetdensity*0.1/fabs(cos(sitheta[i]*deg2rad));
	  range = protTarg_e2r->Eval(ene);
<<<<<<< HEAD
	  ene = protTarg_r2e->Eval(range+tithick);
	  sireco[1][i] = ene;
	}
	else{ //everything is a proton
	  if(pid[i]<0)
	    pid[i] = 0;	  
	  double ene = sienergy[i];
	  //reconstruct energy loss in the al foil
	  double alpha = 90.-detectorangle-sitheta[i];
	  double althick = foilthick*2.70*0.1/cos(alpha*deg2rad);
	  double range = protFoil_e2r->Eval(ene);
	  ene = protFoil_r2e->Eval(range+althick);
	  sireco[0][i] = ene;
	  //reconstruct energy loss in half the target
	  double tithick = targetthick/2*4.50*0.1/cos(sitheta[i]*deg2rad);
	  range = protTarg_e2r->Eval(ene);
	  ene = protTarg_r2e->Eval(range+tithick);
	  sireco[1][i] = ene;
	}
      }
=======
	  ene = protTarg_r2e->Eval(range+tathick);
	  totreco[1][i] = ene;
	}//proton
	if(tinapid[i]==2){
	  //reconstruct energy loss in the al foil
#ifdef KYUHSU
	  alpha = 90.-detectorangle-sitheta[i];
	  double althick = foilthick*foildensity*0.1/fabs(cos(alpha*deg2rad));
	  range = deutFoil_e2r->Eval(ene);
	  ene = deutFoil_r2e->Eval(range+althick);
	  totreco[0][i] = ene;
#else
	  alpha = -90.-detectorangle+sitheta[i];
#endif
	  //reconstruct energy loss in half the target
	  double tathick = targetthick/2*targetdensity*0.1/fabs(cos(sitheta[i]*deg2rad));
	  range = deutTarg_e2r->Eval(ene);
	  ene = deutTarg_r2e->Eval(range+tathick);
	  totreco[1][i] = ene;
	}//deuteron
	//kinematics excitation energy
	//if(tinapid[i]==0){//identified proton
	if(tinapid[i]==0 || tinapid[i]==1){//identified or stopped proton
	  TVector3 recodir(0,0,1);
	  recodir.SetTheta(sitheta[i]*deg2rad);
	  TLorentzVector recoLV(recodir, totreco[1][i]*1000+prot->GetMass()*1000);
	  if(recoLV.Mag()>0)
	    recoLV.SetRho( sqrt( (totreco[1][i]+prot->GetMass())*(totreco[1][i]+prot->GetMass()) - prot->GetMass()*prot->GetMass() )*1000 );
	  if(dpkine.at(0))
	    excen[i] = dpkine.at(0)->GetExcEnergy(recoLV)/1000; // to MeV
>>>>>>> 6d3290167d355101c1169de99269e4a20f9d2022

	  
	  // corr[0]: correct for position on target
	  recodir.SetMagThetaPhi(stripDist[siring[i]], sitheta[i]*deg2rad, detectorphi[i]*deg2rad);
	  TVector3 beampos(targetx, targety, targetz);
	  recodir -= beampos;
          
	  theta[0][i] = recodir.Theta()*rad2deg;
	  recodir.SetMag(1.0);

          recoLV.SetVectM(recodir, totreco[1][i]*1000+prot->GetMass()*1000);
          if(recoLV.Mag()>0)
            recoLV.SetRho( sqrt( (totreco[1][i]+prot->GetMass())*(totreco[1][i]+prot->GetMass()) - prot->GetMass()*prot->GetMass() )*1000 );
          if(dpkine.at(0))
            excencorr[0][i] = dpkine.at(0)->GetExcEnergy(recoLV)/1000; // to MeV


	  // corr[1]: correcto for angle of beam on target
	  // todo: check rotations
	  recodir.RotateY(targeta/1000.0);
	  recodir.RotateX(-targetb/1000.0);

	  theta[1][i] = recodir.Theta()*rad2deg;

          recoLV.SetVectM(recodir, totreco[1][i]*1000+prot->GetMass()*1000);
          if(recoLV.Mag()>0)
            recoLV.SetRho( sqrt( (totreco[1][i]+prot->GetMass())*(totreco[1][i]+prot->GetMass()) - prot->GetMass()*prot->GetMass() )*1000 );
          if(dpkine.at(0))
            excencorr[1][i] = dpkine.at(0)->GetExcEnergy(recoLV)/1000; // to MeV



	}//identified proton
      }//ndet
    }//hadtinadata
    tr->Fill();
    estore->ClearData();
    neve++;
  }//events
  fout->Write();
  fout->Close();
  cout << endl;
  cout << "tdc events " << tdcevt << endl;
  for(int i=0;i<NADC;i++)
    cout << "adc"<<i<<" events " << adcevt[i] << endl;
  cout << "bug events " << bugevt << endl;
  double time_end = get_time();
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
int geo2num(int geo){
  switch(geo){
  case ADC0:
    return 0;
  case ADC1:
    return 1;
  case ADC2:
    return 2;
#ifdef KYUSHU
  case TDC0:
    return 3;
#else
  case ADC3:
    return 3;
  case ADC4:
    return 4;
  case TDC0:
    return 5;
#endif
  default:
    cout <<"disaster, wrong geometrical address "<<endl;
    return -1;
  }
}
map<int,pair<int,double> > ch2stripmap(char* filename){
  ifstream infile;
  infile.open(filename);
  map<int,pair<int,double> > m;
  if(!infile.is_open()){
    cerr << "channel mapping not found, using default" << endl;
    for(int i=0;i<16;i++){
      m[i].first = i;
      m[i].second = 0;
    }
    return m;
  }
  infile.ignore(1000,'\n');
  for(int i=0;i<16;i++){
    int ch,st;
    double th, dth, dist;
    infile >> ch >> st >> th >> dth >> dist;
    infile.ignore(1000,'\n');
    m[ch].first = st;
    m[ch].second = th;
    stripTheta[st]=th;
    stripDist[st]=dist;
  }
  return m;
}
double calcenergyloss(char *settingsfile){
  TEnv* settings = new TEnv(settingsfile);
  double foilthick = settings->GetValue("Foil.Thickness",36.0);
  double targetthick = settings->GetValue("Target.Thickness",2.0);
  double foildensity = settings->GetValue("Foil.Density",2.7);
  double targetdensity = settings->GetValue("Target.Density",4.5);
  double detectordensity = settings->GetValue("Detector.Density",2.33);
  double detectorthick[NDET];
  for(int i=0;i<NDET;i++){
    detectorthick[i] = settings->GetValue(Form("Detector.%d",i),300.0);
  }
  double ebeam = settings->GetValue("Beam.Energy",20.);
  int Abeam = settings->GetValue("Beam.A",99);
  int Zbeam = settings->GetValue("Beam.Z",12);
  ebeam*=Abeam;

  Compound *target;
#ifdef KYUSHU
  Nucleus *ti = new Nucleus(22,26);
  target = new Compound(ti);
#else
  target = new Compound("DPE");
#endif
  if(vlevel>0)
    cout << " target material " << target->GetSymbol() << endl;
  Nucleus *al = new Nucleus(13,14);
  Compound *foil = new Compound(al);
  Nucleus *si = new Nucleus(14,14);
  Compound *dete = new Compound(si);

  Nucleus *prot = new Nucleus(1,0);
  Nucleus *deut = new Nucleus(1,1);

  Reconstruction *protTarg = new Reconstruction(prot, target);
  Reconstruction *deutTarg = new Reconstruction(deut, target);
  double protTarg_range;
  double deutTarg_range;
  targetthick *= targetdensity*0.1;
  protTarg->SetTargetThickness(targetthick/2);
  if(vlevel>0)
    cout << "calculating energy loss of " <<  prot->GetSymbol() << " in half target " << targetthick/2 << " mg/cm^2 " ;
  protTarg_e2r = protTarg->Energy2Range(50,0.1);
  protTarg_r2e = protTarg->Range2Energy(50,0.1);
  protTarg_range = protTarg_e2r->Eval(10);
  if(vlevel>0)
    cout << "10 MeV proton range " << protTarg_range << " mg/cm^2" << endl;

  deutTarg->SetTargetThickness(targetthick/2);
  if(vlevel>0)
    cout << "calculating energy loss of " <<  deut->GetSymbol() << " in half target " << targetthick/2 << " mg/cm^2 " ;
  deutTarg_e2r = deutTarg->Energy2Range(50,0.1);
  deutTarg_r2e = deutTarg->Range2Energy(50,0.1);
  deutTarg_range = deutTarg_e2r->Eval(10);
  if(vlevel>0)
    cout << "10 MeV deuteron range " << deutTarg_range << " mg/cm^2" << endl;

  Reconstruction *protFoil = new Reconstruction(prot, foil);
  Reconstruction *deutFoil = new Reconstruction(deut, foil);
  double protFoil_range;
  double deutFoil_range;
  foilthick *= foildensity*0.1;
  protFoil->SetTargetThickness(foilthick);
  if(vlevel>0)
    cout << "calculating energy loss of " <<  prot->GetSymbol() << " in foil " << foilthick << " mg/cm^2 " ;
  protFoil_e2r = protFoil->Energy2Range(50,0.1);
  protFoil_r2e = protFoil->Range2Energy(50,0.1);
  protFoil_range = protFoil_e2r->Eval(10);
  if(vlevel>0)
    cout << "10 MeV proton range " << protFoil_range << " mg/cm^2" << endl;

  deutFoil->SetTargetThickness(foilthick);
  if(vlevel>0)
    cout << "calculating energy loss of " <<  deut->GetSymbol() << " in foil " << foilthick << " mg/cm^2 " ;
  deutFoil_e2r = deutFoil->Energy2Range(50,0.1);
  deutFoil_r2e = deutFoil->Range2Energy(50,0.1);
  deutFoil_range = deutFoil_e2r->Eval(10);
  if(vlevel>0)
    cout << "10 MeV deuteron range " << deutFoil_range << " mg/cm^2" << endl;

  Reconstruction *protDete = new Reconstruction(prot, dete);
  Reconstruction *deutDete = new Reconstruction(deut, dete);
  double protDete_range;
  double deutDete_range;
  for(int i=0;i<NDET;i++){
    detectorthick[i] *= detectordensity*0.1;
    protDete->SetTargetThickness(detectorthick[i]);
    if(vlevel>0)
      cout << "calculating energy loss of " <<  prot->GetSymbol() << " in detector " << i << " d = " << detectorthick[i] << " mg/cm^2 " ;
    protDete_e2r[i] = protDete->Energy2Range(50,0.1);
    protDete_r2e[i] = protDete->Range2Energy(50,0.1);
    protDete_range = protDete_e2r[i]->Eval(10);
    if(vlevel>0)
      cout << "10 MeV proton range " << protDete_range << " mg/cm^2" << endl;
    
    deutDete->SetTargetThickness(detectorthick[i]);
    if(vlevel>0)
      cout << "calculating energy loss of " <<  deut->GetSymbol() << " in detector " << i << " d = " << detectorthick[i] << " mg/cm^2 " ;
    deutDete_e2r[i] = deutDete->Energy2Range(50,0.1);
    deutDete_r2e[i] = deutDete->Range2Energy(50,0.1);
    deutDete_range = deutDete_e2r[i]->Eval(10);
    if(vlevel>0)
      cout << "10 MeV deuteron range " << deutDete_range << " mg/cm^2" << endl;
  }
  if(vlevel>0)
    cout << "beam A = " << Abeam <<", Z = " << Zbeam <<", E = " << ebeam << endl;
  Nucleus *beam = new Nucleus(Zbeam,Abeam-Zbeam);
  Reconstruction *beamTarg = new Reconstruction(beam, target);
  TSpline3* beamTarg_e2r = beamTarg->Energy2Range(5000,1);
  TSpline3* beamTarg_r2e = beamTarg->Range2Energy(5000,1);
  double range = beamTarg_e2r->Eval(ebeam);
  if(vlevel>0)
    cout << "range " << range << endl;
  double emid = beamTarg_r2e->Eval(range-targetthick/2);
  if(vlevel>0)
    cout << "calculating energy of " <<  beam->GetSymbol() << " in middle of " << targetthick/2 << " mg/cm^2: " << emid << endl;
  
  cout << "calculated energy losses " << endl;
  return emid;
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
  double energies[6] = {0,2,4,6,8,10};
  for(int i=0;i<6;i++)
    dpkine.push_back(new Kinematics(proj,deut,prot,ejec,midtarget,energies[i]));

  cout << "calculated kinematics " << endl;
}
void readpidcuts(const char* filename){
  for(int i=0;i<NDET;i++){
    deutCut[i] = NULL;
    protCut[i] = NULL;
    kineCut = NULL;
    pidSCut[i] = NULL;
    pidLCut[i] = NULL;
  }
  if(strlen(filename)==0){
    cout << "no PID cuts to read"<<endl;
    return;
  }
  TFile* fc = new TFile(filename);
  for(int i=0;i<NDET;i++){
    deutCut[i] = (TCutG*)fc->Get(Form("deut%d",i));
    protCut[i] = (TCutG*)fc->Get(Form("prot%d",i));
    kineCut = (TCutG*)fc->Get(Form("kineCut"));
    pidSCut[i] = (TCutG*)fc->Get(Form("csiS%d",i));
    pidLCut[i] = (TCutG*)fc->Get(Form("csiL%d",i));
  }
  cout << "read PID cuts"<<endl;
}
