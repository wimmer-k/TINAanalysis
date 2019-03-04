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

//energy loss lib
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

#define ADC0 10
#define ADC1 11
#define ADC2 12
#define TDC0 0
#define TIMEREF 0
#define DEBUG 0

#define NDET 4
char* massFile = (char*)"/home/wimmer/progs/eloss/mass.dat";
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
double calcenergyloss(double targetthick, double foilthick, double detthick[NDET], double ebeam);
//kinemaicts
void calckinematics(double midtarget);
//pid cuts
void readpidcuts(char* filename);

int main(int argc, char** argv){
  Int_t RunNumber = -1;
  char* detectorsetup = NULL;
  char* calfile = NULL;
  char* pedfile = NULL;
  char* pidfile = NULL;

  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-r", "RunNumber", &RunNumber);
  interface->Add("-d", "detector settings", &detectorsetup);
  interface->Add("-c", "calibration parameters", &calfile);  
  interface->Add("-p", "pedestal parameters", &pedfile);  
  interface->Add("-i", "pid cut file", &pidfile);  
  interface->CheckFlags(argc, argv);
  if(RunNumber<0){
    cout << "invalid runnumber: " << RunNumber << endl;
    return 1;
  }
  if(detectorsetup==NULL)
    detectorsetup = (char*)"/home/wimmer/TINA/settings/detectorsetup_2um.dat";
  if(calfile==NULL)
    calfile = (char*)"/home/wimmer/TINA/settings/calparams_jul26.dat";
  if(pedfile==NULL)
    pedfile = (char*)"/home/wimmer/TINA/settings/pedestals.dat";
  if(pidfile==NULL)
    pidfile = (char*)"/home/wimmer/TINA/settings/pidcuts.root";
    


  // ========= Define input and output files
  cout << endl << "---> Start analyzing Run# = " << RunNumber << endl;
  char* inputfilename  = (char*)Form("data/run%04d.ridf",RunNumber);
  char* outputfilename = (char*)Form("root/run%04d.root",RunNumber);


  // ========= Define input and output files, trees, histograms etc...
  // -- open input file
  TArtEventStore *estore = new TArtEventStore();
  estore->Open(inputfilename);
  TArtRawEventObject *rawevent = estore->GetRawEventObject();

  // -- open output file
  TFile *fout = new TFile(outputfilename,"RECREATE");

  TTree *tr = new TTree("tr","events");
  Int_t adc[3][32];
  Int_t tdc[64];
  Int_t adcmult;
  Int_t tdcmult;

  //all rings
  Double_t sirien[NDET][16];
 
  //highest per detector (for now)
  Double_t sienergy[NDET];
  Int_t siring[NDET];
  Double_t sitheta[NDET];
  Int_t simult[NDET];
  //energyloss correction for pid
  Double_t sicorr[NDET];
  //reconstructed energy at target center
  Double_t sireco[2][NDET];
  Int_t pid[NDET];

  //todo
  Double_t sibsen[NDET];
  Double_t sibsti[NDET];
  Double_t csien[NDET];
  Double_t csiti[NDET];

  //raw
  tr->Branch("adc",&adc,"adc[3][32]/I");
  tr->Branch("tdc",&tdc,"tdc[64]/I");
  tr->Branch("adcmult",&adcmult,"adcmult/I");
  tr->Branch("tdcmult",&tdcmult,"tdcmult/I");
  //cal
  tr->Branch("sirien",  &sirien,  Form("sirien[%d][16]/D",NDET));
  tr->Branch("sienergy",&sienergy,Form("sienergy[%d]/D",NDET));
  tr->Branch("sitheta", &sitheta, Form("sitheta[%d]/D",NDET));
  tr->Branch("siring",  &siring,  Form("siring[%d]/I",NDET));
  tr->Branch("simult",  &simult,  Form("simult[%d]/I",NDET));

  tr->Branch("sicorr",  &sicorr,  Form("sicorr[%d]/D",NDET));
  tr->Branch("pid",     &pid,     Form("pid[%d]/I",NDET));
  tr->Branch("sireco",  &sireco,  Form("sireco[2][%d]/D",NDET));

  tr->Branch("sibsen",  &sibsen,  Form("sibsen[%d]/D",NDET));
  tr->Branch("sibsti",  &sibsti,  Form("sibsti[%d]/D",NDET));
  tr->Branch("csien",   &csien,   Form("csien[%d]/D",NDET));
  tr->Branch("csiti",   &csiti,   Form("csiti[%d]/D",NDET));


  //define gain and offset parameters for calibration, 3 ADC, 1 TDC
  //pedestals == software thresholds for 3 adcs
  double ped[NDET][32];
  double gain[NDET][32];
  double offs[NDET][32];
  TEnv *cal = new TEnv(calfile);
  TEnv *pedestals = new TEnv(pedfile);
  for(int m=0;m<NDET;m++){
    for(int c=0;c<32;c++){
      gain[m][c] = cal->GetValue(Form("Gain.Module%d.Ch%d",m,c),0.0);
      offs[m][c] = cal->GetValue(Form("Offset.Module%d.Ch%d",m,c),0.0);
      if(m<3)
	ped[m][c] = pedestals->GetValue(Form("Pedestal.Module%d.Ch%d",m,c),0.0);
      //cout << "gain["<<m<<"]["<<c<<"] = " << gain[m][c] << endl;
    } 
    
  }
  
  double detectorthick[NDET];
  TEnv* detsetup = new TEnv(detectorsetup);
  double foilthick = detsetup->GetValue("Foil.Thickness",36.0);
  double targetthick = detsetup->GetValue("Target.Thickness",2.0);
  double detectorangle = detsetup->GetValue("Detector.Angle",50.0);
  for(int i=0;i<NDET;i++){
    detectorthick[i] = detsetup->GetValue(Form("Detector.%d",i),300.0);
  }
  double ebeam = detsetup->GetValue("Beam.Energy",20);
  string channelmap = detsetup->GetValue("Channel.Mapping", "channelmapping.dat");
  //initalize mapping from ADC channel to strip number and angle
  map<int,pair<int,double> > chmap = ch2stripmap((char*)channelmap.c_str());
  for(int i=0;i<16;i++){
    cout << i << "\t" << chmap[i].first<< "\t" << chmap[i].second << endl;
  }

  //prepare the energy reconstruction
  double emid = calcenergyloss(targetthick, foilthick, detectorthick,ebeam);

  //calculate kinematics
  calckinematics(emid);
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
 

  // ========= LOOP for each event

  int neve = 0;  // event number
  TRandom3 *rand = new TRandom3();


  int tdcevt = 0;
  int adcevt[3] = {0,0,0};
  int bugevt = 0;
  TH2F* hraw[3];
  for(int i=0;i<3;i++){
    hraw[i] = new TH2F(Form("adc_%d",i),Form("adc_%d",i),32,0,32,4096,0,4096);
  }

  while(estore->GetNextEvent()){ // from here, sorting the data
    if(DEBUG){
      cout << "-------------------------------------next event segments = "<<rawevent->GetNumSeg()<< endl;
      if(neve == 10)
	break;
    }
    // --- print sorted event number
    if (neve%10000==0)
      cout << neve << " events sorted\r" << flush;

    //clear 
    for(int i=0;i<32;i++){
      for(int j=0;j<3;j++)
	adc[j][i] = 0;
    }
    for(int i=0;i<64;i++){
      tdc[i] = 0;
    }
    adcmult = 0;
    tdcmult = 0;
    for(int i=0;i<NDET;i++){
      for(int j=0;j<16;j++)
	sirien[i][j] = sqrt(-1);
      simult[i] = 0;
      siring[i] = -1;
      sienergy[i] = sqrt(-1);
      sitheta[i] = sqrt(-1);
      sicorr[i] = sqrt(-1);
      sireco[0][i] = sqrt(-1);
      sireco[1][i] = sqrt(-1);
      pid[i] = -1;
      sibsen[i] = -sqrt(-1);
      sibsti[i] = sqrt(-1);
      csien[i] = sqrt(-1);
      csiti[i] = sqrt(-1);
    }

    bool haddata = false;
    bool hadadc[3] = {false,false,false};
    bool hadtdc = false;
    bool hadbug = false;



    // --- reading data
    for(int i=0;i<rawevent->GetNumSeg();i++){
      TArtRawSegmentObject *seg = rawevent->GetSegment(i);
      if(DEBUG){
	cout << "segment " << i << "\tseg->GetAddress " << seg->GetAddress() << "\tset->GetModule "<< seg->GetModule() << "\tnum data = " << seg->GetNumData() << endl;
      }
      for(int j=0;j<seg->GetNumData();j++){
	TArtRawDataObject *d = seg->GetData(j);
	//get geometric address, channel and value
	int geo = d->GetGeo(); 
	int chan = d->GetCh();
	int val = d->GetVal(); 
	if(DEBUG)
	  cout << "geo: " << geo  << "\tchan: " << chan  << "\tval: " << val << endl;  
	int num, det;
	double en;
	double ti;
	switch(geo){
	case ADC0:
	case ADC1:
	case ADC2:
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
	  // adc1 and 2 for 4 YY1, adc0 for backsides and csi
	  if(geo == ADC1)
	    det = chan/16;
	  if(geo == ADC2)
	    det = chan/16+2;
	  if(det>-1&&en>0){
	    simult[det]++;
	    int ring = chmap[chan%16].first; // mapping to real strip number
	    sirien[det][ring] = en;
	    if(siring[det]<0 || en > sienergy[det]){
	      siring[det] = ring;
	      sitheta[det] = chmap[chan%16].second; // mapping to theta
	      sienergy[det] = en;
	      double alpha = 90.-detectorangle-sitheta[det];
	      sicorr[det] = cos(alpha*deg2rad)*en;
	    }
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
	    }
	  }
	  
	  haddata = true;
	  if(DEBUG)
	    cout << neve <<"\tADC\t" << num << "\t" << chan << "\t" << val <<"\tadc[num][chan] = " << adc[num][chan]<< endl;
	  break;
	case TDC0:
	  hadtdc = true;
	  num = geo2num(geo); 
	  //take first hit only, maybe add more later
	  if(tdc[chan]==0)
	    tdc[chan] = val;
	  else
	    continue;
	  tdcmult++;
	  if(DEBUG)
	    cout << neve <<"\tTDC\t"<< chan << "\t" << val <<"\ttdc[chan] = " << tdc[chan]<< endl;
	  haddata = true;
	  break;
	case 99:
	  hadbug = true;
	  if(DEBUG)
	    cout << neve <<"\tBUG\t"<< chan << "\t" << val << endl;
	default:
	  cout << "unknown geo = " << geo << endl; 
	  break;
	}//geo
      }//data
    }//segments
    // if(haddata){
    //   cout << "----------------------" << endl;
    // }
    if(hadtdc){
      tdcevt++;
      for(int i=0;i<NDET;i++){
	int chan = i+2;
	sibsti[i] = gain[geo2num(TDC0)][chan]*(tdc[chan]-tdc[TIMEREF]+rand->Uniform(0,1)) + offs[geo2num(TDC0)][chan];
	chan = i+6;
	csiti[i] = gain[geo2num(TDC0)][chan]*(tdc[chan]-tdc[TIMEREF]+rand->Uniform(0,1)) + offs[geo2num(TDC0)][chan];	
      }
    }
    for(int i=0;i<3;i++){
      if(hadadc[i])
	adcevt[i]++;
    }
    if(hadbug)
      bugevt++;
    
    if(haddata){
      for(int i=0;i<NDET;i++){
	if(pidCut[i]!=NULL && pidCut[i]->IsInside(csien[i],sicorr[i])){
	  pid[i] = 10;
	}
	if(deutCut[i]!=NULL && deutCut[i]->IsInside(sitheta[i],sienergy[i])){
	  pid[i] = 2;
	  double ene = sienergy[i];
	  //reconstruct energy loss in the al foil
	  double alpha = 90.-detectorangle-sitheta[i];
	  double althick = foilthick*2.70*0.1/cos(alpha*deg2rad);
	  double range = deutFoil_e2r->Eval(ene);
	  ene = deutFoil_r2e->Eval(range+althick);
	  sireco[0][i] = ene;
	  //reconstruct energy loss in half the target
	  double tithick = targetthick/2*4.50*0.1/cos(sitheta[i]*deg2rad);
	  range = deutTarg_e2r->Eval(ene);
	  ene = deutTarg_r2e->Eval(range+tithick);
	  sireco[1][i] = ene;
	}
	else if(protCut[i]!=NULL && protCut[i]->IsInside(sitheta[i],sienergy[i])){
	  pid[i] = 1;
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

      tr->Fill();
    }
    estore->ClearData();
    neve++;
  }
  fout->Write();
  fout->Close();
  cout << endl;
  cout << "tdc events " << tdcevt << endl;
  for(int i=0;i<3;i++)
    cout << "adc"<<i<<" events " << adcevt[i] << endl;
  cout << "bug events " << bugevt << endl;
}

int geo2num(int geo){
  switch(geo){
  case ADC0:
    return 0;
  case ADC1:
    return 1;
  case ADC2:
    return 2;
  case TDC0:
    return 3;
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
double calcenergyloss(double targetthick, double foilthick, double detectorthick[NDET], double ebeam){
  Nucleus *ti = new Nucleus(22,26,massFile);
  Compound *target = new Compound(ti);
  //Compound *target = new Compound("2.0DTI");
  Nucleus *al = new Nucleus(13,14,massFile);
  Compound *foil = new Compound(al);
  Nucleus *si = new Nucleus(14,14,massFile);
  Compound *dete = new Compound(si);

  Nucleus *prot = new Nucleus(1,0,massFile);
  Nucleus *deut = new Nucleus(1,1,massFile);

  Reconstruction *protTarg = new Reconstruction(prot, target);
  Reconstruction *deutTarg = new Reconstruction(deut, target);
  double targetdensity = 4.50; //g/cm^3 for Ti 
  //double targetdensity = 3.75; //g/cm^3 for TiH2 
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
  double foildensity = 2.70; //g/cm^3 
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

  Reconstruction *protDete = new Reconstruction(prot, dete);
  Reconstruction *deutDete = new Reconstruction(deut, dete);
  double detedensity = 2.32; //g/cm^3 
  double protDete_range;
  double deutDete_range;
  for(int i=0;i<NDET;i++){
    detectorthick[i] *= detedensity*0.1;
    protDete->SetTargetThickness(detectorthick[i]);
    cout << "calculating energy loss of " <<  prot->GetSymbol() << " in detector " << i << " d = " << detectorthick[i] << " mg/cm^2 " ;
    protDete_e2r[i] = protDete->Energy2Range(50,0.1);
    protDete_r2e[i] = protDete->Range2Energy(50,0.1);
    protDete_range = protDete_e2r[i]->Eval(10);
    cout << "10 MeV proton range " << protDete_range << " mg/cm^2" << endl;
    
    deutDete->SetTargetThickness(detectorthick[i]);
    cout << "calculating energy loss of " <<  deut->GetSymbol() << " in detector " << i << " d = " << detectorthick[i] << " mg/cm^2 " ;
    deutDete_e2r[i] = deutDete->Energy2Range(50,0.1);
    deutDete_r2e[i] = deutDete->Range2Energy(50,0.1);
    deutDete_range = deutDete_e2r[i]->Eval(10);
    cout << "10 MeV deuteron range " << deutDete_range << " mg/cm^2" << endl;
  }

  Nucleus *carb = new Nucleus(6,6,massFile);
  Reconstruction *carbTarg = new Reconstruction(carb, target);
  TSpline3* carbTarg_e2r = carbTarg->Energy2Range(50,0.1);
  TSpline3* carbTarg_r2e = carbTarg->Range2Energy(50,0.1);
  double range = carbTarg_e2r->Eval(ebeam);
  double emid = carbTarg_r2e->Eval(range-targetthick/2);
  cout << "calculating energy of " <<  carb->GetSymbol() << " in middle of " << targetthick/2 << " mg/cm^2: " << emid << endl;
  return emid;
}
void calckinematics(double midtarget){
  Nucleus *prot = new Nucleus(1,0,massFile);
  Nucleus *deut = new Nucleus(1,1,massFile);
  Nucleus *carb = new Nucleus(6,6,massFile);
  Nucleus *ejec = new Nucleus(6,7,massFile);

  ddkine = new Kinematics(carb,deut,deut,carb,midtarget,0);
  ppkine = new Kinematics(carb,prot,prot,carb,midtarget,0);
  double energies[4] = {0,3.089,3.684,3.853};
  for(int i=0;i<4;i++)
    dpkine.push_back(new Kinematics(carb,deut,prot,ejec,midtarget,energies[i]));

  cout << "calculated kinematics " << endl;
}
void readpidcuts(char* filename){
  for(int i=0;i<NDET;i++){
    deutCut[i] = NULL;
    protCut[i] = NULL;
    pidCut[i] = NULL;
  }
  TFile* fc = new TFile(filename);
  for(int i=0;i<NDET;i++){
    deutCut[i] = (TCutG*)fc->Get(Form("deut%d",i));
    protCut[i] = (TCutG*)fc->Get(Form("prot%d",i));
    pidCut[i] = (TCutG*)fc->Get(Form("csi%d",i));
  }
}
