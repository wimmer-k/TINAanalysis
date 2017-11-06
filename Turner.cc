#include "TArtEventStore.hh"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCutG.h>
#include <TTree.h>
#include <TObject.h>
#include <TEnv.h>
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
#define NADC 5

#define BEAMDEVICE 60
#define TINADEVICE 0
#define TOFDET 13
#define PPACDet 15

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
TCutG* pidCut[NDET];

//geometric address to module number
int geo2num(int geo);
//map adc ch to strip number and angle
map<int,pair<int,double> > ch2stripmap(char* filename);
//energy losses
double calcenergyloss(char* filename);
//kinematics
void calckinematics(double midtarget, char* detectorsetup);
//pid cuts
void readpidcuts(char* filename);

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
  char* detectorsetup = NULL;
  char* calfile = NULL;
  char* pedfile = NULL;
  char* pidfile = NULL;
  vlevel = 0;
  int LastEvent =-1;
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-r", "run number", &RunNumber);
  interface->Add("-d", "detector settings", &detectorsetup);
  interface->Add("-c", "calibration parameters", &calfile);  
  interface->Add("-p", "pedestal parameters", &pedfile);  
  interface->Add("-i", "pid cut file", &pidfile);  
  interface->Add("-le", "last event to be read", &LastEvent);  
  interface->Add("-v", "verbose level", &vlevel);  

  //check flags and complain if somethings is missing
  interface->CheckFlags(argc, argv);
  if(RunNumber<0){
    cout << "invalid runnumber: " << RunNumber << endl;
    return 1;
  }
  if(detectorsetup==NULL){
    cout << "Error: no detectorsetup file given" << endl;
    return -99;
  }
  if(calfile==NULL)
    cout << "Warning: no calibration parameters given!" << endl;
  if(pedfile==NULL)
    cout << "Warning: no pedestal parameters given!" << endl;
  if(pidfile==NULL)
    cout << "Warning: no PID cuts given!" << endl;

  //input and output files 
  cout <<" analyzing run " << RunNumber << endl;
  char* inputfilename  = (char*)Form("data/run%04d.ridf",RunNumber);
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
  Double_t csiLen[NDET];
  Double_t csiLti[NDET];

  //total energy for sifsen and csi
  Double_t toten[NDET];
  //reconstructed energy at target center, from sifsen + csi
  Double_t totreco[2][NDET]; //0 only foils, 1 foils and target energy loss reconstructed
  //pid flag: 0 identified (proton), 1 E/theta pid proton, 2 E/theta pid deuteron
  Int_t pid[NDET]; //check still kyushu type, id only for protons


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
  tr->Branch("pid",     &pid,     Form("pid[%d]/I",NDET));
  tr->Branch("totreco", &totreco, Form("totreco[2][%d]/D",NDET));

  tr->Branch("sibsen",  &sibsen,  Form("sibsen[%d]/D",NDET));
  tr->Branch("sibsti",  &sibsti,  Form("sibsti[%d]/D",NDET));
  tr->Branch("csiSen",  &csiSen,  Form("csiSen[%d]/D",NDET));
  tr->Branch("csiSti",  &csiSti,  Form("csiSti[%d]/D",NDET));
  tr->Branch("csiLen",  &csiLen,  Form("csiLen[%d]/D",NDET));
  tr->Branch("csiLti",  &csiLti,  Form("csiLti[%d]/D",NDET));


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
      if(calfile==NULL){
	gain[m][c] = 1.0;
	offs[m][c] = 0.0;
      }
      else{
	gain[m][c] = cal->GetValue(Form("Gain.Module%d.Ch%d",m,c),0.0);
	offs[m][c] = cal->GetValue(Form("Offset.Module%d.Ch%d",m,c),0.0);
      }
      if(m<NADC){
	if(pedfile==NULL)
	  ped[m][c] = 0.0;  
	else
	  ped[m][c] = pedestals->GetValue(Form("Pedestal.Module%d.Ch%d",m,c),0.0);
      }
      //cout << "gain["<<m<<"]["<<c<<"] = " << gain[m][c] << endl;
    } 
    
  }
  
  //reading in detector setup
  double detectorthick[NDET];
  TEnv* detsetup = new TEnv(detectorsetup);
  double foilthick = detsetup->GetValue("Foil.Thickness",36.0);
  double targetthick = detsetup->GetValue("Target.Thickness",2.0);
  double foildensity = detsetup->GetValue("Foil.Density",2.7);
  double targetdensity = detsetup->GetValue("Target.Density",4.5);
  double detectorangle = detsetup->GetValue("Detector.Angle",50.0);
  for(int i=0;i<NDET;i++){
    detectorthick[i] = detsetup->GetValue(Form("Detector.%d",i),300.0);
  }

  string channelmap = detsetup->GetValue("Channel.Mapping", "channelmapping.dat");
  //initalize mapping from ADC channel to strip number and angle
  map<int,pair<int,double> > chmap = ch2stripmap((char*)channelmap.c_str());
  if(vlevel>0){
    for(int i=0;i<16;i++)
      cout << i << "\t" << chmap[i].first<< "\t" << chmap[i].second << endl;
  }

  //prepare the energy reconstruction
  double emid = calcenergyloss(detectorsetup);

  //calculate kinematics
  calckinematics(emid,detectorsetup);
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

  //read in pid cuts
  readpidcuts(pidfile);
  fout->cd();
  for(int i=0;i<NDET;i++){
    if(deutCut[i]!=NULL)
      deutCut[i]->Write();
    if(protCut[i]!=NULL)
      protCut[i]->Write();
    if(pidCut[i]!=NULL)
      pidCut[i]->Write();
  }
 


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

  while(estore->GetNextEvent()){ 
    if(signal_received){
      break;
    }
    if(vlevel >1){
      cout << "-------------------------------------next event segments = "<<rawevent->GetNumSeg()<< endl;
    }
    if(neve == LastEvent)
      break;
    if(neve%10000==0){
      double time_end = get_time();
      cout << neve << " events analyzed, " << (Float_t)neve/(time_end - time_start) << " events/s\r" << flush;
    }
    //clear 
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
      totreco[0][i] = sqrt(-1);
      totreco[1][i] = sqrt(-1);
      pid[i] = -1;
      sibsen[i] = -sqrt(-1);
      sibsti[i] = sqrt(-1);
      csiSen[i] = sqrt(-1);
      csiSti[i] = sqrt(-1);
      csiLen[i] = sqrt(-1);
      csiLti[i] = sqrt(-1);
    }

    bool haddata = false;
    bool hadadc[NADC] = {};
    bool hadtdc = false;
    bool hadbug = false;

    //loop over the segments
    for(int i=0;i<rawevent->GetNumSeg();i++){
      TArtRawSegmentObject *seg = rawevent->GetSegment(i);
      if(vlevel>1){
	cout << "event " << neve << "\tsegment " << i << "\tAddress " << seg->GetAddress() << "\tDevice " << seg->GetDevice() << "\tFocalPlane " << seg->GetFP() << "\tDetector " << seg->GetDetector() << "\tModule "<< seg->GetModule() << "\tnum data = " << seg->GetNumData() << endl;
      }
      if(seg->GetDevice()==BEAMDEVICE){
	//for PID:
	//dev 60, fp 0, det 12, geo 1 (moco + mtdc)  (13 instead of 12 for v1290) 
	//ch 10 F3dia9Pad
	//ch 26 F5dia9Pad
	cout << "BEAM data! " << endl;
	if(seg->GetDetector()==TOFDET){
	  cout << "TOF detector" << endl;
	  for(int j=0;j<seg->GetNumData();j++){
	    TArtRawDataObject *d = seg->GetData(j);
	    //get geometric address, channel and value
	    int geo = d->GetGeo(); 
	    int chan = d->GetCh();
	    int val = d->GetVal(); 
	    if(vlevel>1)
	      cout << "geo: " << geo  << "\tchan: " << chan  << "\tval: " << val << endl;  
	    if(geo==0){
	      if(chan==10)
		cout << "F3diaPad = " << val << endl;
	      if(chan==26)
		cout << "F5diaPad = " << val << endl;
	    }
	  }
	}//tof detector
      }//beam data
      if(seg->GetDevice()==TINADEVICE){
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
		double alpha = 90.-detectorangle-sitheta[det];
		sicorr[det] = cos(alpha*deg2rad)*en;
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
	    haddata = true;
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
	    haddata = true;
	    break;
	  case 99:
	    hadbug = true;
	    //if()
	    cout << neve <<"\tBUG\t"<< chan << "\t" << val << endl;
	  default:
	    cout << "unknown geo = " << geo << endl; 
	    break;
	  }//geo
	}//data
      }//tina data
    }//segments
    // if(haddata){
    //   cout << "----------------------" << endl;
    // }
    if(hadtdc){
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
    //reconstruc energy loss
    if(haddata){
      //calculate total energy
      for(int i=0;i<NDET;i++){
	if(csiLen[i] > csiSen[i])
	  toten[i] = sifsen[i] + csiLen[i];
	else if(csiSen[i] > 0)
	  toten[i] = sifsen[i] + csiSen[i];
      }
      
      for(int i=0;i<NDET;i++){
	if(pidCut[i]!=NULL && pidCut[i]->IsInside(csiLen[i],sicorr[i])){//check change
	  pid[i] = 0;
	  continue;
	}
	double range;
	double alpha;
	if(deutCut[i]!=NULL && deutCut[i]->IsInside(sitheta[i],sifsen[i])){
	  pid[i] = 2;
	  double ene = sifsen[i];
	  //reconstruct energy loss in the al foil
#ifdef KYUHSU
	  alpha = 90.-detectorangle-sitheta[i];
	  double althick = foilthick*foildensity*0.1/cos(alpha*deg2rad);
	  range = deutFoil_e2r->Eval(ene);
	  ene = deutFoil_r2e->Eval(range+althick);
	  totreco[0][i] = ene;
#else
	  alpha = -90.-detectorangle+sitheta[i];
#endif
	  //reconstruct energy loss in half the target
	  double tathick = targetthick/2*targetdensity*0.1/cos(sitheta[i]*deg2rad);
	  range = deutTarg_e2r->Eval(ene);
	  ene = deutTarg_r2e->Eval(range+tathick);
	  totreco[1][i] = ene;
	}
	else if(protCut[i]!=NULL && protCut[i]->IsInside(sitheta[i],sifsen[i])){
	  pid[i] = 1;
	  double ene = sifsen[i];
#ifdef KYUHSU
	  alpha = 90.-detectorangle-sitheta[i];
	  double althick = foilthick*foildensity*0.1/cos(alpha*deg2rad);
	  range = protFoil_e2r->Eval(ene);
	  ene = protFoil_r2e->Eval(range+althick);
	  totreco[0][i] = ene;
#else
	  alpha = -90.-detectorangle+sitheta[i];
#endif
	  //reconstruct energy loss in half the target
	  double tathick = targetthick/2*targetdensity*0.1/cos(sitheta[i]*deg2rad);
	  range = protTarg_e2r->Eval(ene);
	  ene = protTarg_r2e->Eval(range+tathick);
	  totreco[1][i] = ene;
	}
      }//energy loss reconstruction
      tr->Fill();
    }//all segments
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
    double th;
    infile >> ch >> st >> th;
    infile.ignore(1000,'\n');
    m[ch].first = st;
    m[ch].second = th;
  }
  return m;
}
double calcenergyloss(char *detectorsetup){
  TEnv* detsetup = new TEnv(detectorsetup);
  double foilthick = detsetup->GetValue("Foil.Thickness",36.0);
  double targetthick = detsetup->GetValue("Target.Thickness",2.0);
  double foildensity = detsetup->GetValue("Foil.Density",2.7);
  double targetdensity = detsetup->GetValue("Target.Density",4.5);
  double detectordensity = detsetup->GetValue("Detector.Density",2.33);
  double detectorthick[NDET];
  for(int i=0;i<NDET;i++){
    detectorthick[i] = detsetup->GetValue(Form("Detector.%d",i),300.0);
  }
  double ebeam = detsetup->GetValue("Beam.Energy",20.);
  int Abeam = detsetup->GetValue("Beam.A",99);
  int Zbeam = detsetup->GetValue("Beam.Z",12);
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
void calckinematics(double midtarget, char* detectorsetup){
  TEnv* detsetup = new TEnv(detectorsetup);
  int Abeam = detsetup->GetValue("Beam.A",99);
  int Zbeam = detsetup->GetValue("Beam.Z",12);

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
void readpidcuts(char* filename){
  for(int i=0;i<NDET;i++){
    deutCut[i] = NULL;
    protCut[i] = NULL;
    pidCut[i] = NULL;
  }
  if(filename ==NULL){
    cout << "no PID cuts to read"<<endl;
    return;
  }
  TFile* fc = new TFile(filename);
  for(int i=0;i<NDET;i++){
    deutCut[i] = (TCutG*)fc->Get(Form("deut%d",i));
    protCut[i] = (TCutG*)fc->Get(Form("prot%d",i));
    pidCut[i] = (TCutG*)fc->Get(Form("csi%d",i));
  }
  cout << "read PID cuts"<<endl;
}
