{

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

  char dirname[8][256] = {
        "ak3GenJetSpectrum_n22_n12",
        "ak3GenJetSpectrum_n12_n07",
        "ak3GenJetSpectrum_n07_n03",
        "ak3GenJetSpectrum_n03_p03",
        "ak3GenJetSpectrum_p03_p07",
        "ak3GenJetSpectrum_p07_p12",
        "ak3GenJetSpectrum_p12_p22",
        "ak3GenJetSpectrum"};

  double etaWid[8] = {1.0,0.5,0.4,0.6,0.4,0.5,1.0,2.0};

  // retrieve and normalize spectra
  TH1F * spectrum[10][10];
  TH1F * spectrumf[10][10];
  TH1F * events[10][10];
  double nevt[10][10]; 
  for( int i=0; i<10;i++)
  {
    f[i] = new TFile(Form("AnaGENJetR3_Apr2_Z2Pt%s_numEvent640000.root",rpthats[i]));
    for( int j=0; j<8; j++)
    {
      events[i][j] = (TH1F *) f[i]->Get(Form("%s/events",dirname[j]));
      events[i][j]->SetName(Form("events_%d_%d",i,j));
      nevt[i][j] = events[i][j]->GetBinContent(1);
      spectrum[i][j] = (TH1F *) f[i]->Get(Form("%s/%s",dirname[j],"spectrum"));
      spectrum[i][j]->SetName(Form("spectrum_%d_%d",i,j));
      spectrum[i][j]->Sumw2();
      spectrumf[i][j] = (TH1F *) f[i]->Get(Form("%s/%s",dirname[j],"spectrum_fine"));
      spectrumf[i][j]->SetName(Form("spectrum_fine_%d_%d",i,j));
      spectrumf[i][j]->Sumw2();
      for( int bin = 1; bin <= spectrum[i][j]->GetNbinsX(); bin++)
      {
        double ptWid = spectrum[i][j]->GetBinWidth(bin);
        spectrum[i][j]->SetBinContent(bin,
                       spectrum[i][j]->GetBinContent(bin) 
                       / ptWid / etaWid[j] * xsec[i] / nevt[i][j] ); 
        spectrum[i][j]->SetBinError(bin,
                       spectrum[i][j]->GetBinError(bin) 
                       / ptWid / etaWid[j] * xsec[i] / nevt[i][j]  ); 
      }
      for( int bin = 1; bin <= spectrumf[i][j]->GetNbinsX(); bin++)
      {
        double ptWid = spectrumf[i][j]->GetBinWidth(bin);
        spectrumf[i][j]->SetBinContent(bin,
                       spectrumf[i][j]->GetBinContent(bin)
                       / ptWid / etaWid[j] * xsec[i] / nevt[i][j] );
        spectrumf[i][j]->SetBinError(bin,
                       spectrumf[i][j]->GetBinError(bin)
                       / ptWid / etaWid[j] * xsec[i] / nevt[i][j]  );
      }
    }
  }

  char dirnameout[8][256] = {
        "ak3GenJetSpectrum_n22_n12",
        "ak3GenJetSpectrum_n12_n07",
        "ak3GenJetSpectrum_n07_n03",
        "ak3GenJetSpectrum_n03_p03",
        "ak3GenJetSpectrum_p03_p07",
        "ak3GenJetSpectrum_p07_p12",
        "ak3GenJetSpectrum_p12_p22",
        "ak3GenJetSpectrum_n10_p10"};


  // add spectra
  TH1F * spectrum_out[10];
  TH1F * spectrumf_out[10];

  TFile * f_out = new TFile("AnaGENJetR3_Apr2_Z2Combined.root","RECREATE");

  for( int j=0;j<8;j++)
  {
    TDirectory * tdir = f_out->mkdir(dirnameout[j]);
    for( int i=1;i<10;i++) spectrum[0][j]->Add(spectrum[i][j],1.0);
    for( int i=1;i<10;i++) spectrumf[0][j]->Add(spectrumf[i][j],1.0);
    spectrum_out[j] =  (TH1F *) spectrum[0][j]->Clone(Form("JetSpectrum"));
    spectrum_out[j]->SetDirectory(tdir);
    spectrumf_out[j] = (TH1F *) spectrumf[0][j]->Clone(Form("JetSpectrum_Fine"));
    spectrumf_out[j]->SetDirectory(tdir);
  }

  f_out->Write();

}
