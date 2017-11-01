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

//energy loss lib
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

#define NDET 4
#define RAWCSI 0

char* massFile = (char*)"/home/wimmer/progs/eloss/mass.dat";
TSpline3* protDete_l2e[NDET][16];
TSpline3* deutDete_l2e[NDET][16];

TCutG* deutCut[NDET];
TCutG* protCut[NDET];
TCutG* pidCut[NDET];

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
  char* calfile = NULL;
  char* pidfile = NULL;
  int LastEvent =-1;
  int Verbose =0;
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-i", "input files", &InputFiles);
  interface->Add("-o", "output file", &OutFile);    
  interface->Add("-d", "detector settings", &detectorsetup);
  interface->Add("-c", "calibration parameters", &calfile);  
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
  if(calfile==NULL)
    calfile = (char*)"/home/wimmer/TINA/settings/calparams_jul26.dat";
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
  TH2F* hsienergy_csi[NDET];
  TH2F* hsicorr_csi[NDET];
  TH2F* hsibsen_csi[NDET];
  TH2F* hcal_csi[NDET];
  TH2F* hsienergy_theta[NDET];
  TH2F* hsireco_theta[NDET];
  //histograms
  for(unsigned short d=0;d<NDET;d++){
    
    if(RAWCSI){
      hsienergy_csi[d] = new TH2F(Form("hsienergy_csi_%d",d),Form("hsienergy_csi_%d",d),512,0,4096,2000,0,20);hlist->Add(hsienergy_csi[d]);
      hsicorr_csi[d] = new TH2F(Form("hsicorr_csi_%d",d),Form("hsicorr_csi_%d",d),512,0,4096,2000,0,20);hlist->Add(hsicorr_csi[d]);
      hsibsen_csi[d] = new TH2F(Form("hsibsen_csi_%d",d),Form("hsibsen_csi_%d",d),512,0,4096,2000,0,20);hlist->Add(hsibsen_csi[d]);
    }
    else{
      hsienergy_csi[d] = new TH2F(Form("hsienergy_csi_%d",d),Form("hsienergy_csi_%d",d),500,0,20,2000,0,20);hlist->Add(hsienergy_csi[d]);
      hsicorr_csi[d] = new TH2F(Form("hsicorr_csi_%d",d),Form("hsicorr_csi_%d",d),500,0,20,2000,0,20);hlist->Add(hsicorr_csi[d]);
      hsibsen_csi[d] = new TH2F(Form("hsibsen_csi_%d",d),Form("hsibsen_csi_%d",d),500,0,20,2000,0,20);hlist->Add(hsibsen_csi[d]);
    }
    hsienergy_theta[d] = new TH2F(Form("hsienergy_theta_%d",d),Form("hsienergy_theta_%d",d),55,25,80,2000,0,20);hlist->Add(hsienergy_theta[d]);
    hsireco_theta[d] = new TH2F(Form("hsireco_theta_%d",d),Form("hsireco_theta_%d",d),55,25,80,2000,0,20);hlist->Add(hsireco_theta[d]);
    if(RAWCSI)
      hcal_csi[d] = new TH2F(Form("hcal_csi_%d",d),Form("hcal_csi_%d",d),512,0,4096,2000,0,20);
    else
      hcal_csi[d] = new TH2F(Form("hcal_csi_%d",d),Form("hcal_csi_%d",d),500,0,20,2000,0,20);
    hlist->Add(hcal_csi[d]);
  }
  
  int siring[NDET];
  double sitheta[NDET];
  double sienergy[NDET];
  double sibsen[NDET];
  double csien[NDET];
  double sireco[2][NDET];
  tr->SetBranchAddress("siring",&siring);
  tr->SetBranchAddress("sitheta",&sitheta);
  tr->SetBranchAddress("sienergy",&sienergy);
  tr->SetBranchAddress("sibsen",&sibsen);
  tr->SetBranchAddress("sireco",&sireco);
  tr->SetBranchAddress("csien",&csien);
  Double_t nentries = tr->GetEntries();
  cout << nentries << " entries in tree" << endl;
  if(nentries<1)
    return 4;
  if(LastEvent>0)
    nentries = LastEvent;

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
    for(unsigned short d=0;d<NDET;d++){
      siring[d] = -1;
      sitheta[d] = NAN;
      sienergy[d] = NAN;
      sireco[0][d] = NAN;
      sireco[1][d] = NAN;
      sibsen[d] = NAN;
      csien[d] = NAN;
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
    for(unsigned short d=0;d<NDET;d++){
      //cout << sienergy[d] << "\t "<< sibsen[d]<< "\t "<< csien[d] << endl;
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
	  }
	  // double range = 
	  //  calc total energy from sitheta and sienergy
	  //  plot Etot - sienergy
	}
      }
      if(sibsen[d]>0 && csien[d]>0)
	hsibsen_csi[d]->Fill(csien[d],sibsen[d]);
    }
  }
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
  if (sig == SIGINT){
    signal_received = true;
  }
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
  Nucleus *si = new Nucleus(14,14,massFile);
  Compound *dete = new Compound(si);
  
  Nucleus *prot = new Nucleus(1,0,massFile);
  Nucleus *deut = new Nucleus(1,1,massFile);

  Reconstruction *protDete = new Reconstruction(prot, dete);
  Reconstruction *deutDete = new Reconstruction(deut, dete);
  double protDete_range;
  double deutDete_range;
  double detedensity = 2.32; //g/cm^3 
  for(int d=0;d<NDET;d++){
    detectorthick[d] *= detedensity*0.1;
    cout << "calculating energy loss of " <<  prot->GetSymbol() << ", " << deut->GetSymbol() << " in detector " << d << " d = " << detectorthick[d] << " mg/cm^2" << endl;
    for(int s=0;s<16;s++){
      double alpha = 90.-detectorangle-strmap[s];
      if(DEBUG)
	cout << s << "\t" << alpha << "\t" << cos(alpha*deg2rad) << endl;
      protDete->SetTargetThickness(detectorthick[d]/cos(alpha*deg2rad));
      protDete_l2e[d][s] = protDete->EnergyLoss2Energy(50,0.1);
      deutDete->SetTargetThickness(detectorthick[d]/cos(alpha*deg2rad));
      deutDete_l2e[d][s] = deutDete->EnergyLoss2Energy(50,0.1);
    }
  }

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

