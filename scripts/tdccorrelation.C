void tdccorrelation(int run){
  TFile *f = new TFile(Form("/home/daq/TiNA/root/run%04d.root",run));
  TTree *t = (TTree*)f->Get("tr");
  TH1F *h1[13];
  TH2F *h2[13][13];
  int i1[13];
  int i2[13][13];
  TFile *fout = new TFile(Form("/home/daq/TiNA/root/out%04d.root",run),"RECREATE");
  for(int c=16;c<16+12;c++){
    cout << c << "\t" << c-16 << endl;
    t->Draw(Form("tdc[%d]>>hh(4000,8e4,12e4)",c));
    h1[c-16] = (TH1F*)gROOT->FindObject("hh");
    h1[c-16]->SetName(Form("tdc%d",c));
    i1[c-16] = h1[c-16]->Integral();
    cout << i1[c-16] << endl;
    h1[c-16]->Write();
  }
  int startcsitdc =16;
  TH2F* corr = new TH2F("corr","corr",12,0,12,12,0,12);
  for(int c1=startcsitdc;c1<startcsitdc+12;c1++){
    for(int c2=startcsitdc;c2<startcsitdc+12;c2++){
      cout << "c1 = " << c1 << ", " << c1-startcsitdc << "\tc2 = "<< c2 << ", " << c2-startcsitdc<< endl;
      t->Draw(Form("tdc[%d]:tdc[%d]>>hhh(4000,8e4,12e4,4000,8e4,12e4)",c1,c2));
      h2[c1-startcsitdc][c2-startcsitdc] = (TH2F*)gROOT->FindObject("hhh");
      h2[c1-startcsitdc][c2-startcsitdc]->SetName(Form("tdc%d_%d",c1,c2));
      i2[c1-startcsitdc][c2-startcsitdc] = h2[c1-startcsitdc][c2-startcsitdc]->Integral();
      cout << i2[c1-startcsitdc][c2-startcsitdc] << endl;
      corr->Fill(c1-startcsitdc,c2-startcsitdc,i2[c1-startcsitdc][c2-startcsitdc]);
      h2[c1-startcsitdc][c2-startcsitdc]->Write();
    }
  }
  corr->Write();
  fout->Close();

}
void plot(char * file){
  TFile *f = new TFile(file);
  TH2F* corr = new TH2F("corr","corr",12,0,12,12,0,12);
  TH2F *h2[12][12];
  int i2[12][12];

  int startcsitdc =16;
  for(int c1=startcsitdc;c1<startcsitdc+12;c1++){
    for(int c2=startcsitdc;c2<startcsitdc+12;c2++){
      cout << "c1 = " << c1 << ", " << c1-startcsitdc << "\tc2 = "<< c2 << ", " << c2-startcsitdc;
      h2[c1-startcsitdc][c2-startcsitdc] = (TH2F*)f->Get(Form("tdc%d_%d",c1,c2));
      i2[c1-startcsitdc][c2-startcsitdc] = h2[c1-startcsitdc][c2-startcsitdc]->Integral();
      cout << "\t" << i2[c1-startcsitdc][c2-startcsitdc] << endl;
      corr->Fill(c1-startcsitdc,c2-startcsitdc,i2[c1-startcsitdc][c2-startcsitdc]);      
    }
  }
  corr->Draw("colztext");
}
