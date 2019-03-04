<<<<<<< HEAD
const int NDET = 4;
void csiplotcal(int run){
  
  TFile* f;
  if(run==0)
    f = new TFile("hist/csical_add.root");
  else
    f = new TFile(Form("hist/run%04d.root",run));
  TH2F* hen_csi[NDET];
  TCanvas* c = new TCanvas("c","c",800,600);
  c->Divide(2,2);
  TFile* fc = new TFile("settings/pidcuts_Feb11.root");
  TCutG* gc[NDET];
  fc->ls();
  for(unsigned short d=0;d<NDET;d++){
    cout << "getting " << Form("hsienergy_csi_%d",d) << endl;
    hen_csi[d] = (TH2F*)f->Get(Form("hsienergy_csi_%d",d));
    if(hen_csi[d]==NULL)
      continue;
    cout << "got " << Form("hsienergy_csi_%d",d) << endl;
    c->cd(1+d);
    hen_csi[d]->GetYaxis()->SetRangeUser(0,10);
    hen_csi[d]->Draw("colz");
    gc[d] = (TCutG*)fc->Get(Form("csi%d",d));
    if(gc[d]==NULL)
      continue;
    gc[d]->Draw("same");
  } 
  hen_csi[0]->GetXaxis()->SetRangeUser(0,4);
  hen_csi[1]->GetXaxis()->SetRangeUser(0,5.5);
  hen_csi[2]->GetXaxis()->SetRangeUser(0,3.5);
  hen_csi[3]->GetXaxis()->SetRangeUser(0,5.5);
}
void csiplotraw(int run){
  
  TFile* f;
  if(run==0)
    f = new TFile("hist/csirawcut_add.root");
  else
    f = new TFile(Form("hist/run%04d.root",run));
  TH2F* hen_csi[NDET];
  TCanvas* c = new TCanvas("c","c",800,600);
  c->Divide(2,2);
  TFile* fc = new TFile("settings/rawpidcuts_Feb11.root");
  TCutG* gc[NDET];
  fc->ls();
  for(unsigned short d=0;d<NDET;d++){
    cout << "getting " << Form("hsienergy_csi_%d",d) << endl;
    hen_csi[d] = (TH2F*)f->Get(Form("hsienergy_csi_%d",d));
    if(hen_csi[d]==NULL)
      continue;
    cout << "got " << Form("hsienergy_csi_%d",d) << endl;
    c->cd(1+d);
    hen_csi[d]->GetYaxis()->SetRangeUser(0,10);
    hen_csi[d]->Draw("colz");
    gc[d] = (TCutG*)fc->Get(Form("csi%d",d));
    if(gc[d]==NULL)
      continue;
    gc[d]->Draw("same");
  } 
}
void csical(){
  const int NDET = 4;
  TFile *f = new TFile("hist/csirawcut_add.root");
=======
{
  const int NDET = 6;
  TFile *f = new TFile("hist/ike0035.root");
>>>>>>> 6d3290167d355101c1169de99269e4a20f9d2022
  TH2F* hcal_csi[NDET];
  TH1F* hcal_csi_px[NDET];
  TF1* fcal[NDET];
  TCanvas *c = new TCanvas("c","c",800,600);
<<<<<<< HEAD
  c->Divide(2,2);
  int adc =0;
  TEnv* params = new TEnv("settings/calparams_feb11_run36.dat");
=======
  c->Divide(3,2);
  int adc =3;
  TEnv* params = new TEnv("settings/calparams_nov08.dat");
>>>>>>> 6d3290167d355101c1169de99269e4a20f9d2022
  for(unsigned short d=0;d<NDET;d++){
    hcal_csi[d] = (TH2F*)f->Get(Form("hcal_csiL_%d",d));
    if(hcal_csi[d]==NULL)
      continue;
    //temp rebin
    hcal_csi[d]->RebinX(10);
    hcal_csi_px[d] = (TH1F*)hcal_csi[d]->ProfileX(Form("hcal_csiL_px_%d",d));
    c->cd(1+d);
    hcal_csi_px[d]->Draw();
    fcal[d] = new TF1(Form("fcal_%d",d),"pol1",0,4000);
    hcal_csi_px[d]->Fit(fcal[d],"R");
    params->SetValue(Form("Gain.Module%d.Ch%d",adc,d+22), fcal[d]->GetParameter(1));
    params->SetValue(Form("Offset.Module%d.Ch%d",adc,d+22), fcal[d]->GetParameter(0));
    
  }
  params->SaveLevel(kEnvLocal);
}
void csires(){
  gStyle->SetOptFit(1);
  const int NDET = 4;
  TFile *f = new TFile("hist/csicalcut_add.root");
  TH1F* hres_csi[NDET];
  TF1* fres[NDET];
  TCanvas *c = new TCanvas("c","c",800,600);
  c->Divide(2,2);
  for(unsigned short d=0;d<NDET;d++){
    hres_csi[d] = (TH1F*)f->Get(Form("hres_csi_%d",d));
    if(hres_csi[d]==NULL)
      continue;
    c->cd(1+d);
    //   if(d==2)
      hres_csi[d]->Rebin(40);
      //    else
      //      hres_csi[d]->Rebin(5);
    hres_csi[d]->Draw();
    fbinw = hres_csi[d]->GetBinWidth(1);
    fnpeaks =1;
    fres[d] = new TF1(Form("fres_%d",d),multgausbg,-2,2,5);
    double pars[5] = {0,0,hres_csi[d]->Integral(),0,0.1};
    fres[d]->SetParameters(pars);
    fres[d]->FixParameter(0,0);
    fres[d]->FixParameter(1,0);
    //fres[d]->Draw("same");
    hres_csi[d]->Fit(fres[d],"R");
    fres[d]->Draw("same");
  }
}
