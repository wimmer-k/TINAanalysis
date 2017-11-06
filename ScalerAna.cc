#include "TArtEventStore.hh"
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <TEnv.h>
#include <TCanvas.h> 
#include <TGraph.h> 
#include <TAxis.h> 
#include <TLegend.h> 
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
#define NDET 6
int vlevel =0;
//signal handling, enables ctrl-c to exit nicely
bool signal_received = false;
void signalhandler(int sig);
double get_time();
string* readscalermap(char *filename);
void printscalers(int runnumber, char* scalermap);
int main(int argc, char** argv){
  double time_start = get_time();  
  signal(SIGINT,signalhandler);
  Int_t RunNumber = -1;
  vlevel = 0;
  int LastEvent =-1;
  char* ScalerMap = NULL;
  bool PrintOnly = false;
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-r", "run number", &RunNumber);
  interface->Add("-p", "print only, assumes scaNNNN.root exists", &PrintOnly);
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
    cout << "Error: no scaler mapping file given" << endl;
    return -99;
  }
  if(PrintOnly){
    printscalers(RunNumber,ScalerMap);
    return 0;
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
void printscalers(int RunNumber, char* ScalerMap){
  //top left trigger live, trigger raw
  //middle left Si FS, individual and sum
  //middle right Si BS, individual and sum
  //bottom left CsI small, individual and sum
  //bottom right CsI large, individual and sum

  TFile* fsca = new TFile(Form("./root/sca%04d.root",RunNumber));
  TGraph *g[32];
  string* scalername = readscalermap(ScalerMap);
  for(int i=0;i<32;i++){
    g[i] = (TGraph*)fsca->Get(Form("scaler_%d",i));
    g[i]->SetTitle(scalername[i].c_str());
    g[i]->GetXaxis()->SetTitle("time (s)");
    g[i]->GetYaxis()->SetTitle("rate (1/s)");
  }
  TCanvas *c = new TCanvas("scalers","scalers",600,900);
  c->Divide(2,3);
  TLegend *leg[6];
  for(int i=0;i<6;i++){
    leg[i] = new TLegend(0.8,0.8,0.98,0.98);
    leg[i]->SetTextSize(0.03);
  }
  g[0]->SetLineColor(3);
  g[1]->SetLineColor(2);
  leg[0]->AddEntry(g[0],scalername[0].c_str(),"L");
  leg[0]->AddEntry(g[1],scalername[1].c_str(),"L");
  c->cd(1);
  g[1]->Draw("AL");
  g[0]->Draw("L");

  for(int j=0;j<4;j++){
    c->cd(3+j);
    g[4+j]->SetLineColor(1);
    leg[2+j]->AddEntry(g[4+j],scalername[4+j].c_str(),"L");
    g[4+j]->Draw("AL");
    for(int i=0;i<NDET;i++){
      g[8+j*NDET+i]->SetLineColor(2+i/3);
      g[8+j*NDET+i]->SetLineStyle(1+i%3);
      leg[2+j]->AddEntry(g[8+j*NDET+i],scalername[8+j*NDET+i].c_str(),"L");
      g[8+j*NDET+i]->Draw("L");
    }
    leg[2+j]->Draw();
  }

  c->SaveAs(Form("scaler/sca%04d.pdf",RunNumber));

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
