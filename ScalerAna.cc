#include "TArtEventStore.hh"
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <TEnv.h>
#include <TGraph.h> 
#include <TAxis.h> 
#include <iostream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <sys/time.h>

//command line interface library
#include "CommandLineInterface.hh"
using namespace std;

#define SCALER_ID 4
#define CLOCK1KHZ 2
int vlevel =0;
//signal handling, enables ctrl-c to exit nicely
bool signal_received = false;
void signalhandler(int sig);
double get_time();
string* readscalermap(char *filename);
int main(int argc, char** argv){
  double time_start = get_time();  
  signal(SIGINT,signalhandler);
  Int_t RunNumber = -1;
  vlevel = 0;
  int LastEvent =-1;
  char* ScalerMap = NULL;
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-r", "run number", &RunNumber);
  interface->Add("-m", "scaler channel mapping", &ScalerMap);
  interface->Add("-le", "last event to be read", &LastEvent);  
  interface->Add("-v", "verbose level", &vlevel);  

  //check flags and complain if somethings is missing
  interface->CheckFlags(argc, argv);
  if(RunNumber<0){
    cout << "invalid runnumber: " << RunNumber << endl;
    return 1;
  }
  if(ScalerMap==NULL){
    cout << "Error: no detectorsetup file given" << endl;
    return -99;
  }
  //input and output files 
  cout <<" analyzing run " << RunNumber << endl;
  char* inputfilename  = (char*)Form("data/run%04d.ridf",RunNumber);
  char* outputfilename = (char*)Form("root/sca%04d.root",RunNumber);
  TArtEventStore *estore = new TArtEventStore();
  estore->Open(inputfilename);
  TArtRawEventObject *rawevent = estore->GetRawEventObject();
  TFile *fout = new TFile(outputfilename,"RECREATE");

  int neve = 0;  // event number
  long long int values[32];
  double rate[32];
  TGraph *g[32];
  string* scalername = readscalermap(ScalerMap);
  double x[1] = {0};
  double y[1] = {0};
  for(int i=0;i<32;i++){
    if(vlevel>0)
      cout << scalername[i] << endl;
    values[i] = 0;
    rate[i] = 0;
    g[i] = new TGraph(1,x,y);
    g[i]->SetName(Form("scaler_%d",i));
    g[i]->GetXaxis()->SetTitle("time (s)");
    g[i]->GetYaxis()->SetTitle("rate (1/s)");
  }

  int scalerevents = 0;
  while(estore->GetNextEvent()){ 
    if(signal_received){
      break;
    }
    if(vlevel >1){
      cout << "-------------------------------------next event segments = "<<rawevent->GetNumScaler()<< endl;
    }
    if(neve == LastEvent)
      break;
    if(neve%100000==0){
      double time_end = get_time();
      cout << neve << " events analyzed, " << (Float_t)neve/(time_end - time_start) << " events/s\r" << flush;
    }
    for(int i=0;i<rawevent->GetNumScaler();i++){
      TArtRawScalerObject *sca = rawevent->GetScaler(i);
      if(sca->GetScalerID()!=SCALER_ID)
	cout << "unknown scaler ID: sca->GetScalerID() = " << sca->GetScalerID()<<endl;
      else
	scalerevents++;
      for(int j=0;j<sca->GetNumChannel();j++){
	if(vlevel>1&& sca->GetScaler(j)!=0)
	  cout << "sca->GetScaler("<<j<<") = " << sca->GetScaler(j)<<endl;
	rate[j] = sca->GetScaler(j) - values[j];
	values[j] = (long long int)sca->GetScaler(j);
      }
      for(int j=0;j<sca->GetNumChannel();j++){
	if(j==CLOCK1KHZ)
	  continue;
	if(vlevel>1)
	  cout << j <<"\t"<<rate[j] ;
	rate[j]/=rate[CLOCK1KHZ]/1e3;
	if(vlevel>1)
	  cout << "\t"<<rate[j] << "\t"<<values[CLOCK1KHZ]<< endl;
	g[j]->SetPoint(g[j]->GetN(),values[CLOCK1KHZ]/1e3,rate[j]);
      }
      // update the clock
      if(vlevel>0)
	  cout << CLOCK1KHZ <<"\t"<<rate[CLOCK1KHZ] ;
      rate[CLOCK1KHZ]/=rate[CLOCK1KHZ]/1e3;
	if(vlevel>0)
	  cout << "\t"<<rate[CLOCK1KHZ] << "\t"<<(long long int)values[CLOCK1KHZ]<< endl;
	g[CLOCK1KHZ]->SetPoint(g[CLOCK1KHZ]->GetN(),(long long int)values[CLOCK1KHZ]/(long long int)1e3,rate[CLOCK1KHZ]);
    }//segments
    estore->ClearData();
    neve++;
  }//events
  cout << endl;
  cout << "totals "<<scalerevents << " scalerevents,\t" << neve<< " all events" << endl;
  for(int i=0;i<32;i++){
    if(values[i]>0)
      cout << Form("ch[%2d] =\t",i) << values[i]  << "\t" << scalername[i] << endl;
    g[i]->Write("",TObject::kOverwrite);
  }

  fout->Close();
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
string* readscalermap(char *filename){
  string* scalername = new string[32];
  TEnv* sca = new TEnv(filename);
  for(int i=0;i<32;i++){
    scalername[i] = sca->GetValue(Form("Scaler.%d",i),"");
  }
  return scalername;
}
