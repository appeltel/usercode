{

  char energy[256] = "7000";

  //const int MAXDIR = 61;
  const int MAXDIR = 5;

  // cross sections from the Pythia simulations:

  // Z2 cross sections

  // for 2.76 TeV
  //double xsec[10] = {6.198e+01, 5.135e-02, 9.742e-03, 9.250e-04,
  //                   8.812e-05, 9.771e-06, 1.256e-06, 1.795e-07,
  //                   2.708e-08, 4.735e-09};

  // for 5.02 TeV
  //double xsec[10] = {6.782e+01, 1.405e-01, 2.995e-02, 3.334e-03,
  //                  3.800e-04, 5.138e-05, 8.251e-06, 1.525e-06,
  //                   3.141e-07, 9.160e-08};

  // for 7 TeV
  double xsec[10] = {7.131e+01, 2.359e-01, 5.307e-02, 6.354e-03,
                     7.837e-04, 1.152e-04, 2.015e-05, 4.094e-06,
                     9.351e-07, 3.221e-07};
  

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

  char dirname[MAXDIR][256] = {
/*        "ak3GenJetSpectrum_n22_n12",
        "ak3GenJetSpectrum_n12_n07",
        "ak3GenJetSpectrum_n07_n03",
        "ak3GenJetSpectrum_n03_p03",
        "ak3GenJetSpectrum_p03_p07",
        "ak3GenJetSpectrum_p07_p12",
        "ak3GenJetSpectrum_p12_p22",
        "ak3GenJetSpectrum_n10_p10",

        "ak4GenJetSpectrum_n22_n12",
        "ak4GenJetSpectrum_n12_n07",
        "ak4GenJetSpectrum_n07_n03",
        "ak4GenJetSpectrum_n03_p03",
        "ak4GenJetSpectrum_p03_p07",
        "ak4GenJetSpectrum_p07_p12",
        "ak4GenJetSpectrum_p12_p22",
        "ak4GenJetSpectrum_n10_p10",

        "ak5GenJetSpectrum_n22_n12",
        "ak5GenJetSpectrum_n12_n07",
        "ak5GenJetSpectrum_n07_n03",
        "ak5GenJetSpectrum_n03_p03",
        "ak5GenJetSpectrum_p03_p07",
        "ak5GenJetSpectrum_p07_p12",
        "ak5GenJetSpectrum_p12_p22",
        "ak5GenJetSpectrum_n10_p10",

        "ak7GenJetSpectrum_n22_n12",
        "ak7GenJetSpectrum_n12_n07",
        "ak7GenJetSpectrum_n07_n03",
        "ak7GenJetSpectrum_n03_p03",
        "ak7GenJetSpectrum_p03_p07",
        "ak7GenJetSpectrum_p07_p12",
        "ak7GenJetSpectrum_p12_p22",
        "ak7GenJetSpectrum_n10_p10",

        "ak3GenJetSpectrum_QCD10001_00_05",
        "ak3GenJetSpectrum_QCD10001_05_10",
        "ak3GenJetSpectrum_QCD10001_10_15",
        "ak3GenJetSpectrum_QCD10001_15_20",
        "ak3GenJetSpectrum_QCD10001_20_25",
        "ak3GenJetSpectrum_QCD10001_25_30",

        "ak4GenJetSpectrum_QCD10001_00_05",
        "ak4GenJetSpectrum_QCD10001_05_10",
        "ak4GenJetSpectrum_QCD10001_10_15",
        "ak4GenJetSpectrum_QCD10001_15_20",
        "ak4GenJetSpectrum_QCD10001_20_25",
        "ak4GenJetSpectrum_QCD10001_25_30",

        "ak5GenJetSpectrum_QCD10001_00_05",
        "ak5GenJetSpectrum_QCD10001_05_10",
        "ak5GenJetSpectrum_QCD10001_10_15",
        "ak5GenJetSpectrum_QCD10001_15_20",
        "ak5GenJetSpectrum_QCD10001_20_25",
        "ak5GenJetSpectrum_QCD10001_25_30",

        "ak7GenJetSpectrum_QCD10001_00_05",
        "ak7GenJetSpectrum_QCD10001_05_10",
        "ak7GenJetSpectrum_QCD10001_10_15",
        "ak7GenJetSpectrum_QCD10001_15_20",
        "ak7GenJetSpectrum_QCD10001_20_25",
        "ak7GenJetSpectrum_QCD10001_25_30",

        "ak7GenJetSpectrum_QCD11004_00_05",
        "ak7GenJetSpectrum_QCD11004_05_10",
        "ak7GenJetSpectrum_QCD11004_10_15",
        "ak7GenJetSpectrum_QCD11004_15_20",
        "ak7GenJetSpectrum_QCD11004_20_25"
    */
 
        "ak3GenJetSpectrum_QCD11004_00_05",
        "ak3GenJetSpectrum_QCD11004_05_10",
        "ak3GenJetSpectrum_QCD11004_10_15",
        "ak3GenJetSpectrum_QCD11004_15_20",
        "ak3GenJetSpectrum_QCD11004_20_25"  

    };

  double etaWid[MAXDIR] = { /*
                       1.0,0.5,0.4,0.6,0.4,0.5,1.0,2.0,
                       1.0,0.5,0.4,0.6,0.4,0.5,1.0,2.0,
                       1.0,0.5,0.4,0.6,0.4,0.5,1.0,2.0,
                       1.0,0.5,0.4,0.6,0.4,0.5,1.0,2.0,
                       1.0,1.0,1.0,1.0,1.0,1.0,
                       1.0,1.0,1.0,1.0,1.0,1.0,
                       1.0,1.0,1.0,1.0,1.0,1.0,
                       1.0,1.0,1.0,1.0,1.0,1.0, */
                       1.0,1.0,1.0,1.0,1.0

};

  // retrieve and normalize spectra
  TH1F * spectrum[10][MAXDIR];
  TH1F * spectrumf[10][MAXDIR];
  TH1F * events[10][MAXDIR];
  double nevt[10][MAXDIR]; 
  for( int i=0; i<10;i++)
  {
    f[i] = new TFile(Form("AnaGENJetR357_%sGeV_Apr16_Z2Pt%s_numEvent640000.root",energy,rpthats[i]));
    for( int j=0; j<MAXDIR; j++)
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

  // add spectra
  TH1F * spectrum_out[MAXDIR];
  TH1F * spectrumf_out[MAXDIR];

  TFile * f_out = new TFile(Form("AnaGENJetR357_%sGeV_Apr16_Z2Combined.root",energy),"RECREATE");

  for( int j=0;j<MAXDIR;j++)
  {
    TDirectory * tdir = f_out->mkdir(dirname[j]);
    for( int i=1;i<10;i++) spectrum[0][j]->Add(spectrum[i][j],1.0);
    for( int i=1;i<10;i++) spectrumf[0][j]->Add(spectrumf[i][j],1.0);
    spectrum_out[j] =  (TH1F *) spectrum[0][j]->Clone(Form("JetSpectrum"));
    spectrum_out[j]->SetDirectory(tdir);
    spectrumf_out[j] = (TH1F *) spectrumf[0][j]->Clone(Form("JetSpectrum_Fine"));
    spectrumf_out[j]->SetDirectory(tdir);
  }

  f_out->Write();

}
