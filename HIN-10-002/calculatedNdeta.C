void calculatedNdeta() {

  // === General Settings ===

  // Note centrality is by 2.5% bins, low bin is inclusive, high bin is exclusive
  // i.e. for 0-10% use minCent=0, maxCent=4. 


  int minCent[12] = {0,2,4,6,8,10,12,14,16,20,24,28};
  int maxCent[12] = {2,4,6,8,10,12,14,16,20,24,28,32};
  
  double ptbins[21]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,12.0,16.0,22.0};

  // absolute values of eta are used, negative eta values are included
  // for full eta range set etaMin=0 etaMax=2.4
  Double_t  etaMin = 0.0;
  Double_t  etaMax = 0.4;

double dndeta[12]= {1612,1313,1082,893.9,731.6,596.8, 482.3,383.9, 266.5, 153.2,79.5, 36.3 };
double dndetaerr[12] = {55,45, 38, 31.4, 26.3, 23.1, 18.7, 16.2, 14.2, 9.7, 6.57, 3.99};

 
double dndeta_calc[12];
double dndeta_calcsyserr[12] = { 0.004, 0.01, 0.017, 0.023, 0.031, 0.041, 0.05, 0.061, 0.08, 0.12, 0.16, 0.21};
double dndeta_calcerr[12];
 
  // Input File
  TFile *f = new TFile("LowPtSpectrum_fine_Full_b.root");
  TFile *fcorr = new TFile("validation3D_HydjetNoWedge_100k_flowSQ_goodTightdz10chi40_vertexZ10.root");

  // =========================

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);


  set_plot_style();

  char dirname[256] = "spectrumGoodTightdz10chi40";
  TH1D * spec[12];
  TH1D * eff[12];
  TH1D * fak[12];
  TH1D * sim[12];
  TH1D * rec[12];

  TGraphErrors * gSpectrum[12];
  TF1 * f1[12];
  TF1 * f1x[12];
  TF1 * f2[12];

  Double_t  meanpt[12];
  Double_t  meanpterr[12];

  
    Double_t px[100],py[100],pxe[100],pye[100];
  for( int c=0; c<12; c++)
  {
    spec[c] = getSpectrum( f, Form("%s/tracks3D",dirname), minCent[c], maxCent[c], etaMin, etaMax ); 
    eff[c] = getSpectrum( fcorr, "hitrkPixelEffAnalyzer/heff3D", minCent[c], maxCent[c], etaMin, etaMax );
    sim[c] = getSpectrum( fcorr, "hitrkPixelEffAnalyzer/hsim3D", minCent[c], maxCent[c], etaMin, etaMax );
    fak[c] = getSpectrum( fcorr, "hitrkPixelEffAnalyzer/hfak3D", minCent[c], maxCent[c], etaMin, etaMax );
    rec[c] = getSpectrum( fcorr, "hitrkPixelEffAnalyzer/hrec3D", minCent[c], maxCent[c], etaMin, etaMax );

    // determine correction factors


    TH1F * nevts = (TH1F *) f->Get(Form("%s/nevts",dirname));
    double Nevt = nevts->Integral( nevts->GetXaxis()->FindBin(minCent[c]+0.001), 
                                   nevts->GetXaxis()->FindBin(maxCent[c]-0.001) );

    double maxy = 10000.;
    double miny = 0.001;
 
    //cout << "Centrality " << (int)(minCent[c]*2.5) << "-" << (int)(maxCent[c]*2.5)
    //     << ": \n\n";

    //cout << "Nevt = " << Nevt << endl << endl;

    for( int i=0; i<=18;i++)
    {
      double ptmin = ptbins[i]; double ptmax = ptbins[i+1];

      double iptmin = spec[c]->FindBin(ptmin+1e-3);
      double iptmax = spec[c]->FindBin(ptmax-1e-3);
      double icptmin = eff[c]->FindBin(ptmin+1e-3);
      double icptmax = eff[c]->FindBin(ptmax-1e-3);
      double pt = 0.;
      for( int k=iptmin;k<=iptmax;k++) pt += spec[c]->GetBinCenter(k) * spec[c]->GetBinContent(k);
      pt /= spec[c]->Integral(iptmin,iptmax);

      double yielderr;
      double yield = spec[c]->IntegralAndError(iptmin,iptmax,yielderr);
      yield /= (ptmax-ptmin) * Nevt * 2 * (etaMax-etaMin); 
      yielderr /= (ptmax-ptmin) * Nevt * 2 * (etaMax-etaMin); 

      double efmin = eff[c]->GetBinContent(icptmin)/sim[c]->GetBinContent(icptmin);
      double efmax = eff[c]->GetBinContent(icptmax)/sim[c]->GetBinContent(icptmax);
      double famin = fak[c]->GetBinContent(icptmin)/rec[c]->GetBinContent(icptmin);
      double famax = fak[c]->GetBinContent(icptmax)/rec[c]->GetBinContent(icptmax);

      double ef = (pt-ptmin)*(efmax-efmin)/(ptmax-ptmin)+efmin;
      double fa = (pt-ptmin)*(famax-famin)/(ptmax-ptmin)+famin;
      yield *= (1-fa)/ef;
      
      px[i] = pt;
      py[i] = yield;
      if( py[i] > maxy ) maxy = py[i];
      if( py[i] < miny ) miny = py[i];
      pxe[i] = 0.;
      pye[i] = yielderr ;
      //cout << px[i] << "  " << py[i] << "  " << pye[i] << endl;

    }  

    gSpectrum[c] = new TGraphErrors( 18, px, py, pxe, pye );
    gSpectrum[c]->SetMarkerStyle(20);
//    f1[c] = new TF1(Form("f1_%d",c),PowerFuncScaled,0,6,3);
//    f1[c]->SetParNames("n","p0","dN/dy");
//    f1[c]->SetParameters(28,5, 2000./ (double)(c+1) );
//    f1[c]->SetParameters(28,5,dndeta[c]);
//    f1[c]->FixParameter(2,dndeta[c]);
//f1[c] = new TF1(Form("f1_%d",c),"6.28318531 * x * [2]*pow((1+x/[1]),-[0])",0.3,5.1);
f1[c] = new TF1(Form("f1_%d",c),"6.28318531 * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.3,3.1);
//f1[c]->SetParameters(10.8327,2.221289,1201);
f1[c]->SetParameters(10.8327,2.221289,dndeta[c]);
    f1[c]->SetLineColor(kBlue);
    f1[c]->SetLineWidth(1);
    f1[c]->SetLineStyle(1);
    //cout << endl;
    gSpectrum[c]->Fit(Form("f1_%d",c),"0RNSQ");
meanpt[c] = 2. * f1[c]->GetParameter(1) /
               ( f1[c]->GetParameter(0) - 3. );

f1x[c] = new TF1(Form("f1_%d",c),"6.28318531 * x * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.0,3.1);
    f1x[c]->SetParameters(f1[c]->GetParameter(0),f1[c]->GetParameter(1), f1[c]->GetParameter(2) );

    double numerator = f1x[c]->Integral(0.0,0.3);
    double denominator = f1[c]->Integral(0.0,0.3);
    for( int p = 0; p<gSpectrum[c]->GetN(); p++)
    {
       Double_t  ppx, ppy;
       gSpectrum[c]->GetPoint(p,ppx,ppy);
       if ( 0.3 < ppx && ppx < 0.6 ) numerator += ppx * ppy * 0.1;
       if ( 0.6 < ppx && ppx < 1.2 ) numerator += ppx * ppy * 0.2;
       if ( 1.2 < ppx && ppx < 2.0 ) numerator += ppx * ppy * 0.4;
       if ( 2.0 < ppx && ppx < 4.0 ) numerator += ppx * ppy * 0.5;
       if ( 4.0 < ppx && ppx < 6.0 ) numerator += ppx * ppy * 1.0;
       if ( 0.3 < ppx && ppx < 0.6 ) denominator += ppy * 0.1;
       if ( 0.6 < ppx && ppx < 1.2 ) denominator += ppy * 0.2;
       if ( 1.2 < ppx && ppx < 2.0 ) denominator += ppy * 0.4;
       if ( 2.0 < ppx && ppx < 4.0 ) denominator += ppy * 0.5;
       if ( 4.0 < ppx && ppx < 6.0 ) denominator += ppy * 1.0;
       cout << ppx << "  " << numerator/denominator  <<  endl;
    }

    cout << endl;

   // meanpt[c] = f1[c]->Mean(0,6);
    meanpt[c] = numerator / denominator;
    
  //  dndeta_calc[c] = f1[c]->Integral(0,6);
      dndeta_calc[c] = denominator;
   // dndeta_calcerr[c] = f1[c]->IntegralError(0,6); 

    dndeta_calcerr[c] = 0.06 * dndeta_calc[c];
    //cout << endl << "Chi2/ndof: " << f1[c]->GetChisquare() / (double)(f1[c]->GetNDF() ) << endl;

    //cout << "dN/deta - Calculated " << dndeta_calc[c] << "  HIN-10-001: " << dndeta[c] << endl;

 //   TF1 * me = new TF1(Form("me_%d",c),PowerFuncME,0,6,3);
 //   double intmeerror =  me->IntegralError(0,6,r->GetParams(), r->GetCovarianceMatrix()->GetMatrixArray() );
 //   double intme = me->Integral(0,6,r->GetParams() );
 //   double intspec = f1[c]->Integral(0,6);
 //   double intspecerr = f1[c]->IntegralError(0,6);
 //   meanpterr[c] = intme / intspec * sqrt( intmeerror * intmeerror / intme / intme
  //                                        + intspecerr * intspecerr / intspec / intspec );
    meanpterr[c] = 0.;
 //     cout   << "  Mean pT: " << meanpt[c] <<  "  Error: " << meanpterr[c] << endl << endl << endl;
    
  }

//  TLegend *legPt = new TLegend(0.45,0.7,0.89,0.89);
//  legPt->SetFillColor(0); legPt->SetBorderSize(0);
//  legPt->AddEntry(spec96,"96 eff","lp");
//  legPt->AddEntry(spec99,"99 eff","lp");
//  legPt->AddEntry(spec102,"102 eff","lp");
  TH1D* hDumPt = new TH1D("hDumPt",";p_{T} [GeV/c];N",40,0.,20.); hDumPt->SetMaximum(4000);
  hDumPt->GetXaxis()->CenterTitle(); hDumPt->GetYaxis()->SetTitleOffset(1.1);
  hDumPt->SetMinimum(0.00001);
  TCanvas *c8 = new TCanvas("c8","Spectrum",600,600);
  hDumPtEff=(TH1F*) hDumPt->Clone("hDumPtEff");
  hDumPtEff->SetTitle(Form("Spectrum %4.1f < |#eta| < %4.1f, %4.0f -%4.0f Pct Centrality", 
                          etaMin,etaMax,
                           minCent[0]*2.5,maxCent[0]*2.5));
  c8->cd(); 
 // c8->SetLogy();
//  hDumPtEff->SetMaximum(maxy*1.2);
//  hDumPtEff->SetMinimum(miny>0.?miny*0.5:0.0001);
  hDumPtEff->Draw(); 

  gSpectrum[0]->Draw("p");
  f1[0]->SetRange(0.0,6.0);
  f1[0]->Draw("same");
 
  TH1D* hDumCen = new TH1D("hDumCen",";Centrality (%);<p_{T}> (GeV/c)",40,0.,80.); hDumCen->SetMaximum(1.0);
  hDumCen->GetXaxis()->CenterTitle(); hDumCen->GetYaxis()->SetTitleOffset(1.1);

  TCanvas *c1 = new TCanvas("c1","Mean Pt",600,600);
 
  c1->cd();
  hDumCen->Draw();

  Double_t mx[12], mxe[12], mye[12];
  for(int i=0;i<12;i++) 
  { 
    mx[i] = (double)(maxCent[i]+minCent[i]) * 1.25;
    mxe[i] = 0; mye[i] = meanpterr[i];
  }  

  TGraphErrors * gmean = new TGraphErrors(12, mx, meanpt, mxe, mye);
  gmean->SetMarkerStyle(20);
  gmean->SetMarkerColor(kRed);
  gmean->Draw("p");

  cout << "meanpt values \n\n";
  for(int c=0; c<12; c++) cout <<c << "  " << meanpt[c] << endl;
  cout << "\n\n";

  cout << "dN/deta values \n\n";
  cout << "cent HIN-10-001  error    HIN-10-002  error  ratio   error\n";
  for( int c =0; c<12; c++)
  {
    cout << c << "  ";
    cout << dndeta[c] << "  ";
    cout << dndetaerr[c] << "  ";
    cout << dndeta_calc[c] << "  ";
    cout << sqrt( dndeta_calcerr[c]*dndeta_calcerr[c] + 0.0036* dndeta_calc[c]*dndeta_calc[c])  << endl;
  }
  TH1D* hDumCen2 = new TH1D("hDumCen2",";Centrality (%); dN/d#eta ratio",40,0.,80.); 
  hDumCen2->GetXaxis()->CenterTitle(); hDumCen2->GetYaxis()->SetTitleOffset(1.1);
  hDumCen2->SetMaximum(1.3);
  hDumCen2->SetMinimum(0.7);
  TCanvas *c3 = new TCanvas("c3","dNdeta ratio",600,600);

  c3->cd();
  hDumCen2->Draw();

  Double_t rx[12], rxe[12], rye[12], ry[12];
  for(int i=0;i<12;i++)
  {
    rx[i] = (double)(maxCent[i]+minCent[i]) * 1.25;
    rxe[i] = 0; 
    ry[i] = dndeta_calc[i] / dndeta[i];
    double nsyserr = dndeta_calc[i] * dndeta_calcsyserr[i];
    double nerr = sqrt( dndeta_calcerr[i]*dndeta_calcerr[i] + nsyserr*nsyserr) ;
    double derr = dndetaerr[i];
    rye[i] = ry[i] * sqrt( nerr*nerr/dndeta_calc[i]/dndeta_calc[i] + derr*derr/ dndeta[i]/dndeta[i] ); 
  }

  TGraphErrors * gdndeta = new TGraphErrors(12, rx, ry, rxe, rye);
  gdndeta->SetMarkerStyle(20);
  gdndeta->SetMarkerColor(kRed);
  gdndeta->Draw("p");

  TLatex *  tex = new TLatex(10,1.2,"Calculated / HIN-10-001");
  tex->SetTextSize(0.035);
  tex->Draw();
  TLatex *  tex = new TLatex(10,1.15,"99% Trigger Efficiency");
  tex->SetTextSize(0.035);
  tex->Draw();
  TLatex *  tex = new TLatex(10,1.1,"Tsallis Fit");
  tex->SetTextSize(0.035);
  tex->Draw();



}

TH1D * getSpectrum(TFile * f, char * type, int minCent, int maxCent, double etaMin, double etaMax )
{
  TH3F *h = (TH3F*) f->Get(Form("%s",type));
  h->GetZaxis()->SetRangeUser( minCent+0.001, maxCent-0.001);
  TH2F * ch = (TH2F*) h->Project3D("yxe");

  TH1D* hPt;
  TH1D* hPtN;
  Int_t etabin10m=ch->GetXaxis()->FindBin(etaMin+0.01);
  Int_t etabin10p=ch->GetXaxis()->FindBin(etaMax-0.01);
  Int_t Netabin10p=ch->GetXaxis()->FindBin(-etaMin-0.01);
  Int_t Netabin10m=ch->GetXaxis()->FindBin(-etaMax+0.01);
  hPt = (TH1D*) ch->ProjectionY(Form("hPt_%d",rand()),etabin10m,etabin10p,"e");
  hPtN = (TH1D*) ch->ProjectionY("hPtN",Netabin10m,Netabin10p,"e");
  hPt->Add( hPtN, 1.0);

  delete h;
  delete ch;
  delete hPtN;

  return hPt;
}



void set_plot_style() {

  // nicer colz plots

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}


Double_t PowerFunc(Double_t * x, Double_t* par) {
  Double_t f;
  Double_t n = par[0];
  Double_t p0 = par[1];
  Double_t dNdy = par[2];
  f = ((n - 1.)*(n - 2.))/(2.*TMath::Pi()*p0*p0);

  f = f*dNdy;

  f = f*pow((1. + fabs(x[0])/p0),-n);
  return f;
}


Double_t PowerFuncScaled(Double_t * x, Double_t* par) {
  Double_t f;
  Double_t n = par[0];
  Double_t p0 = par[1];
  Double_t dNdy = par[2];
  f = ((n - 1.)*(n - 2.))/(2.*TMath::Pi()*p0*p0);

  f = f*dNdy;

  f = f*pow((1. + fabs(x[0])/p0),-n);
  f*=2.*TMath::Pi()*x[0];
  return f;
}

Double_t PowerFuncME(Double_t * x, Double_t* par) {
  Double_t f;
  Double_t n = par[0];
  Double_t p0 = par[1];
  Double_t dNdy = par[2];
  f = ((n - 1.)*(n - 2.))/(2.*TMath::Pi()*p0*p0);

  f = f*dNdy;

  f = f*pow((1. + fabs(x[0])/p0),-n);
  f*=2.*TMath::Pi()*x[0]*x[0];
  return f;
}


