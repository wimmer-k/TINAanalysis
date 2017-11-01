#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLine.h"
#include "TEnv.h"
using namespace std;

int lbinw;
int lnpeaks;
TFile* f =NULL;
TH2F* h2d = NULL;
TH1F* h = NULL;
TEnv* params;
double llowcut = 1000;
double lhighcut = 4096;
double lfitrange = 50;
double lthresh = 0.5;
int pedestalsigmas = 5;
double searchwidth[3][32];
Double_t lmultgausbg(Double_t *x, Double_t *par);
void loadchannel(int adc, int ch);
void setfilename(char* filename){lfilename = filename;}
void setsearchsigmas();
vector<double> fitlines(int np=3, double searchsigma = 10, bool draw=true, bool noBG = false);
void cal(int adc,int ch, bool draw =true);
double fitpedestal(int adc,int ch, bool draw = true);

char* cfilename = (char*)"settings/calparams_jul26.dat";
char* pfilename = (char*)"settings/pedestals.dat";
char* lfilename = (char*)"root/run0014.root";

void run(){
  params = new TEnv(cfilename);
  lfitrange = 100;
  lthresh =0.1;
  cal(0,16,false);
  cal(0,17,false);
  cal(0,18,false);
  lfitrange = 100;
  for(int adc=1;adc<3;adc++){
    for(int ch=0;ch<32;ch++){
      lthresh =0.001;
      if(adc==1 && (ch==14 || ch==17 || ch==28 || ch==30))
	continue;
      if(adc==2 && (ch==2 || ch==10 || ch==12 || ch==14 || ch==26 || ch==28 || ch==30))
	continue;
      cal(adc,ch,false);
    }
  }
  lthresh = 0.00005;
  lfitrange = 150;
  cal(1,14,false);
  cal(2,10,false);
  cal(2,12,false);
  lthresh = 0.0006;  
  cal(2,26,false);
  cal(2,28,false);
  lfitrange = 50;
  lthresh = 0.00005;
  cal(1,28,false);
  params->SaveLevel(kEnvLocal);
}
void pedestals(){
  params = new TEnv(pfilename);
  setfilename(lfilename);
  TCanvas *c= new TCanvas("c","c",900,400);
  c->Divide(3,1);
  TLine* line[3][32];
  for(int adc=0;adc<3;adc++){
    c->cd(1+adc);
    loadchannel(adc,0);
    TH2F *h2dc = (TH2F*)h2d->Clone(Form("clone%d",adc));
    h2dc->GetYaxis()->SetRangeUser(0,500);
    h2dc->Draw("colz");
    for(int ch = 0;ch<32;ch++){
      double thresh = fitpedestal(adc,ch,false);
      if(thresh<0)
	continue;
      line[adc][ch] = new TLine(ch,thresh,ch+1,thresh);
      line[adc][ch]->SetLineWidth(2);
      line[adc][ch]->SetLineColor(1);
      line[adc][ch]->Draw();
      params->SetValue(Form("Pedestal.Module%d.Ch%d",adc,ch), thresh);

    }
  }

  params->SaveLevel(kEnvLocal);
}
double fitpedestal(int adc,int ch, bool draw){
  loadchannel(adc,ch);
  double searchsigma = 10;
  TCanvas *c;
  if(draw){
    c= new TCanvas("c","c",400,400);
    c->cd();
    h->Draw();
  }
    h->GetXaxis()->SetRangeUser(0,300);
  // peak search
  TSpectrum *sp = new TSpectrum(1);
  Int_t nfound;
  if(draw)
    nfound= sp->Search(h,searchsigma,"",lthresh);
  else
    nfound= sp->Search(h,searchsigma,"goff",lthresh);
  if(nfound!=1){
    cout << "searched 1 peak, found " << nfound << endl;
    return -1;
  }
  TLine *line;
  Float_t *xpeaks= sp->GetPositionX();
  lnpeaks = 1;
  double slope = 0;
  double bg = 0;
  double sigma;
  double content;
  lbinw = h->GetBinWidth(1);
  TF1* fu =new TF1("fu",lmultgausbg, xpeaks[0]-lfitrange,xpeaks[0]+lfitrange,5);
  sigma = 10;
  content = h->Integral(h->FindBin(xpeaks[0]-lfitrange),h->FindBin(xpeaks[0]+lfitrange));
  //cout << "parameters " << bg<<"\t"<<slope<<"\t"<<content<<"\t"<<xpeaks[i]<<"\t"<<sigma<<endl;
  fu->SetParameters(bg,slope,content,xpeaks[0],sigma);
  fu->SetParName(0,"BgConstant");
  fu->SetParName(1,"BgSlope   ");
  fu->SetParName(2,"Content   ");
  fu->SetParName(3,"Mean      ");
  fu->SetParName(4,"Sigma     ");
  h->Fit(fu,"Rnq");
  double thresh = fu->GetParameter(3) + pedestalsigmas*fu->GetParameter(4);
  if(draw){
    fu->Draw("same");
    line = new TLine(thresh,0,thresh,h->GetMaximum());
    line->SetLineColor(3);
    line->Draw();
  }
  return thresh;
}
void cal(int adc,int ch, bool draw){
  loadchannel(adc,ch);
  TCanvas *c;
  if(draw){
    c= new TCanvas("c","c",800,400);
    c->Divide(2,1);
    c->cd(1);
    gPad->SetLogy();
    h->Draw();
  }
  vector<double> ene;
  ene.push_back(0);
  ene.push_back(5.486);//dummy
  vector<double> raw;
  vector<TLine*> line;
  line.resize(ene.size());
  llowcut = 0;
  lhighcut = 4096;
  raw = fitlines(ene.size(),searchwidth[adc][ch],draw,true);
  cout << "raw.size() " << raw.size() << endl;
  if(adc==2&&ch==12)
    raw[1] = 1896;
  if(adc==2&&ch==28)
    raw[1] = 1940;
  if(draw){
    h->GetXaxis()->SetRangeUser(0,5000);
    for(unsigned short int i=0;i<ene.size();i++){
      //cout << "raw["<<i<<"] = " << raw[i] << endl;
      line[i] = new TLine(raw[i],0,raw[i],h->GetMaximum());
      line[i]->SetLineColor(3);
      line[i]->Draw();
    }
  }
  
  TGraph* g = new TGraph(ene.size(),&raw[0],&ene[0]);
  TF1* fu = new TF1("pol","pol1",0,5000);
  g->Fit(fu,"Rqn");
  cout << "adc = " << adc << ", ch = " << ch << ", gain = " << fu->GetParameter(1) << ", offset = " << fu->GetParameter(0)<< endl;
  if(draw){
    
    c->cd(2);
    g->Draw("AP*");
    fu->Draw("same");
  }
  else{
    params->SetValue(Form("Gain.Module%d.Ch%d",adc,ch), fu->GetParameter(1));
    params->SetValue(Form("Offset.Module%d.Ch%d",adc,ch), fu->GetParameter(0));
  }
}
vector<double> fitlines(int np, double searchsigma, bool draw, bool noBG){
  //set the range of the histogram to cut off the noise
  h->GetXaxis()->SetRangeUser(llowcut,lhighcut);
  if(draw)
    h->Draw(); 

  // peak search
  TSpectrum *sp = new TSpectrum(np);
  Int_t nfound;
  if(draw)
    nfound= sp->Search(h,searchsigma,"",lthresh);
  else
    nfound= sp->Search(h,searchsigma,"goff",lthresh);
  if(nfound!=np && nfound!=0){
    cout << " PROBLEM NOT "<< np<<" PEAKS found, but " << nfound << endl;
    if(draw)
      nfound= sp->Search(h,searchsigma+1,"",lthresh);
    else
      nfound= sp->Search(h,searchsigma+1,"goff",lthresh);
    if(nfound!=np && nfound!=0){
      h->Rebin(2);
      h->GetXaxis()->SetRangeUser(llowcut,lhighcut);
      if(draw)
	nfound= sp->Search(h,2,"",lthresh);
      else
	nfound= sp->Search(h,2,"goff",lthresh);
    }
  }
  Float_t *xpeaks= sp->GetPositionX();
  lnpeaks = 1;
  double slope;
  double bg;
  double sigma;
  double content;
  lbinw = h->GetBinWidth(1);
  TF1* fu[3];
  vector<double> mean;
  vector<int> order;
  for(int i=0;i<np;i++){
   order.push_back(i);
   //new function with range around the peak
   fu[i] = new TF1(Form("fu_%d",i),lmultgausbg, xpeaks[i]-lfitrange,xpeaks[i]+lfitrange,5);
   //start parameters
   slope = h->GetBinContent(h->FindBin(xpeaks[i]-lfitrange)) - h->GetBinContent(h->FindBin(xpeaks[i]+lfitrange));
   slope/=(xpeaks[i]-lfitrange-(xpeaks[i]+lfitrange));
   bg = h->GetBinContent(h->FindBin(xpeaks[i]-lfitrange)) + h->GetBinContent(h->FindBin(xpeaks[i]+lfitrange));
   bg-=slope*(xpeaks[i]-lfitrange+(xpeaks[i]+lfitrange));
   bg/=2;
      
   sigma = 10;
   content = h->Integral(h->FindBin(xpeaks[i]-lfitrange),h->FindBin(xpeaks[i]+lfitrange));
   //cout << "parameters " << bg<<"\t"<<slope<<"\t"<<content<<"\t"<<xpeaks[i]<<"\t"<<sigma<<endl;
   if(noBG){
     fu[i]->SetParameters(0,0,content,xpeaks[i],sigma);
     fu[i]->FixParameter(0,0);
     fu[i]->FixParameter(1,0);
   }
   else
     fu[i]->SetParameters(bg,slope,content,xpeaks[i],sigma);
   fu[i]->SetParName(0,"BgConstant");
   fu[i]->SetParName(1,"BgSlope   ");
   fu[i]->SetParName(2,"Content   ");
   fu[i]->SetParName(3,"Mean      ");
   fu[i]->SetParName(4,"Sigma     ");
   h->Fit(fu[i],"Rnq");
   mean.push_back(fu[i]->GetParameter(3));
  }

  /*
  for(int i=0;i<np;i++){
    cout << " before mean " << mean[i] << "  position " << order[i] << endl;
  }
  */
  int SwapCount = 0;
  for(int i=0;i<np;i++){
    for(int j=0;j<(np-1);j++){
      if( mean[j] > mean[j+1] ){
	double temp = mean[j];
	int itemp = order[j];
	mean[j] = mean[j+1];
	order[j] = order[j+1];
	mean[j+1] = temp;
	order[j+1] = itemp;
	SwapCount++;	
      }      
    }

    if(SwapCount == 0)
      break;
    else
      SwapCount = 0;
  }
  /*
  for(int i=0;i<np;i++){
    cout << " after mean " << mean[i] << "  position " << order[i] << endl;
  }
  */



  if(draw){
    for(int i=0;i<np;i++)
      fu[i]->Draw("same");
  }
  
  return mean;
  
}

void loadchannel(int adc, int ch){
  setsearchsigmas();
  if(f==NULL)
    f = new TFile(lfilename);
  h2d = (TH2F*)f->Get(Form("adc_%d",adc));
  h = (TH1F*)h2d->ProjectionY(Form("adc_%d_ch_%d",adc,ch),ch+1,ch+1);
}
void setsearchsigmas(){
  for(int adc=0;adc<3;adc++){
    for(int ch=0;ch<32;ch++){
      searchwidth[adc][ch] = 20;
    }
  }
  searchwidth[1][19] = 20;
  searchwidth[1][22] = 20;
  searchwidth[1][24] = 20;
  searchwidth[1][28] = 40;
  //searchwidth[2][12] = 40;
}
void checkparam(){
  double chan[3][32];
  double gain[3][32];
  double offs[3][32];
  TEnv *cal = new TEnv(cfilename);
  TGraph* g[2][3];
  TCanvas *ca = new TCanvas("ca","ca",600,300);
  ca->Divide(2,1);
  for(int m=0;m<3;m++){
    for(int c=0;c<32;c++){
      chan[m][c] = c;
      gain[m][c] = cal->GetValue(Form("Gain.Module%d.Ch%d",m,c),1.0);
      offs[m][c] = cal->GetValue(Form("Offset.Module%d.Ch%d",m,c),0.0);
    } 
    g[0][m] = new TGraph(32,chan[m],gain[m]);
    g[1][m] = new TGraph(32,chan[m],offs[m]);
    g[0][0]->GetYaxis()->SetRangeUser(0.002,0.004);
    g[1][0]->GetYaxis()->SetRangeUser(-1,1);
    for(int i=0;i<2;i++){
      ca->cd(1+i);
      if(m==0)
	g[i][0]->Draw("AL");
      g[i][m]->SetLineColor(1+m);
      g[i][m]->Draw("L");
    }  
  }
}


Double_t lmultgausbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  /*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss0 content
  par[3]   gauss0 mean
  par[4]   gauss0 width
  par[5]   gauss1 content
  par[6]   gauss1 mean
  par[7]   gauss1 width
  ......
  */
  Double_t result = par[0] + par[1]*x[0];

  for (Int_t p=0;p<lnpeaks;p++) {
    Double_t norm  = par[3*p+2];
    Double_t mean  = par[3*p+3];
    Double_t sigma = par[3*p+4];
    if(sigma==0)
      continue;
    arg = (x[0]-mean)/(sqrt2*sigma);
    result += lbinw/(sqrt2pi*sigma) * norm * exp(-arg*arg);
  }
  

  return result;
}
void cal_avfexp(int adc,int ch, bool draw){
  loadchannel(adc,ch);
  TCanvas *c;
  if(draw){
    c= new TCanvas("c","c",800,400);
    c->Divide(2,1);
    c->cd(1);
    h->Draw();
  }
  vector<double> ene;
  //if(ch<16)
  if(adc<2)
    ene.push_back(3183);
  else
    ene.push_back(4788);
  ene.push_back(5486);
  ene.push_back(5805);
  vector<double> raw;
  vector<TLine*> line;
  line.resize(ene.size());
  lthresh =0.5;
  if(adc==2)
    lthresh =0.3;
  llowcut = 1000;
  if(adc==0&&ch==29)
    llowcut = 1700;
  lhighcut = 4096;
  if(adc==1&&ch==25)
    lhighcut = 3500;
  raw = fitlines(ene.size(),searchwidth[adc][ch],draw);
  cout << "raw.size() " << raw.size() << endl;
  if(draw){
    h->GetXaxis()->SetRangeUser(1000,5000);
    for(unsigned short int i=0;i<ene.size();i++){
      //cout << "raw["<<i<<"] = " << raw[i] << endl;
      line[i] = new TLine(raw[i],0,raw[i],h->GetMaximum());
      line[i]->SetLineColor(3);
      line[i]->Draw();
    }
  }
  
  TGraph* g = new TGraph(ene.size(),&raw[0],&ene[0]);
  TF1* fu = new TF1("pol","pol1",0,5000);
  g->Fit(fu,"Rqn");
  cout << "adc = " << adc << ", ch = " << ch << ", gain = " << fu->GetParameter(1) << ", offset = " << fu->GetParameter(0)<< endl;
  if(draw){
    
    c->cd(2);
    g->Draw("AP*");
    fu->Draw("same");
  }
  else{
    params->SetValue(Form("Gain.Module%d.Ch%d",adc,ch), fu->GetParameter(1));
    params->SetValue(Form("Offset.Module%d.Ch%d",adc,ch), fu->GetParameter(0));
  }
}

