{
  const int NDET = 6;
  TFile *f = new TFile("hist/ike0035.root");
  TH2F* hcal_csi[NDET];
  TH1F* hcal_csi_px[NDET];
  TF1* fcal[NDET];
  TCanvas *c = new TCanvas("c","c",800,600);
  c->Divide(3,2);
  int adc =3;
  TEnv* params = new TEnv("settings/calparams_nov08.dat");
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
