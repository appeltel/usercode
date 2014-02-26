{

  bool doWideBinCorrection = false;

  gStyle->SetOptStat(0);

  // cross sections from the Pythia simulations:
  double xsec[10] = {6.782e+01, 1.405e-01, 2.995e-02, 3.334e-03,
                     3.800e-04, 5.138e-05, 8.251e-06, 1.525e-06,
                     3.141e-07, 9.160e-08};
  // need to fix first entry given as 0-inf. Subtract all 
  // other xsecs to get 0-20
  double xsec20toinf = 0.0;
  double totalxsec = xsec[0];
  for( int i=1;i<10;i++) xsec20toinf += xsec[i];
  xsec[0] = xsec[0] - xsec20toinf;

  // open files and get spectra and event counts 
  TFile * f[10];
  char rpthats[10][256] = {"0to20","20to30","30to50","50to80",
                           "80to120","120to170","170to230",
                           "230to300","300to380","380toinf"};
  TH1F * spectrum[10];
  TH1F * events[10];
  double nevt[10]; 
  for( int i=0; i<10;i++)
  {
    f[i] = new TFile(Form("AnaGENJetR3_Feb25_Z2_Pt%s_numEvent640000.root",rpthats[i]));
    spectrum[i] = (TH1F *) f[i]->Get("ak3GenJetSpectrum/spectrum");
    spectrum[i]->SetName(Form("spectrum%d",i));
    spectrum[i]->Sumw2();
    events[i] = (TH1F *) f[i]->Get("ak3GenJetSpectrum/events");
    events[i]->SetName(Form("events%d",i));
    nevt[i] = events[i]->GetBinContent(1);
  }

  // scale spectra by the xsec, nevt, and bin width
  for( int i=0; i<10; i++)
  {
    for( int bin = 1; bin < spectrum[i]->GetNbinsX(); bin++)
    {
       double ptWid = spectrum[i]->GetBinWidth(bin);
       double etaWid = 2.0;
       spectrum[i]->SetBinContent(bin,
                      spectrum[i]->GetBinContent(bin) 
                      / ptWid / etaWid * xsec[i] / nevt[i] ); 
       spectrum[i]->SetBinError(bin,
                      spectrum[i]->GetBinError(bin) 
                      / ptWid / etaWid * xsec[i] / nevt[i]  ); 
    }
  }

  // add spectra
  for( int i=1; i<10; i++ ) spectrum[0]->Add(spectrum[i],1.0);

  // wide binning correction
  if( doWideBinCorrection )
  {
    TFile * fcorr = new TFile("correctionfactor.root");
    TH1F * hcorr = (TH1F*) fcorr->Get("hcorrppEtaBin-10_10");
    for( int bin = 1; bin < spectrum[0]->GetNbinsX(); bin++)
    {
      double pt = spectrum[0]->GetBinCenter(bin);
      int cbin = hcorr->GetXaxis()->FindBin(pt);
      if( hcorr->GetBinContent(cbin) > 0.0 && hcorr->GetBinContent(cbin) < 10.00)
      {
         cout << pt << "   " << hcorr->GetBinContent(cbin) << endl;
         spectrum[0]->SetBinContent(bin,
                        spectrum[0]->GetBinContent(bin)
                        / hcorr->GetBinContent(cbin) );
      }
    }
  }  
  // print results 
  for( int bin = 1; bin < spectrum[0]->GetNbinsX(); bin++)
  {
    cout << spectrum[0]->GetBinCenter(bin) << "  "
         << spectrum[0]->GetBinContent(bin) << "  "
         << spectrum[0]->GetBinError(bin) << endl;
  }

  // plot results
  TCanvas * c1 = new TCanvas("c1","c1",500,500);
  c1->SetLogy();
  c1->cd();

  TH1F * hDum = new TH1F("hDum",";p_{T} [GeV/c];d^{2}#sigma / dp_{T}d#eta",10,37,700);
  hDum->SetMaximum(1.e-3);
  hDum->SetMinimum(1.e-11);
  hDum->Draw();

  spectrum[0]->SetMarkerStyle(20);
  spectrum[0]->Draw("same");

}
