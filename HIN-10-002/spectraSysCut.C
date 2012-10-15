void spectraSysCut() {

  // === General Settings ===

  // Note centrality is by 2.5% bins, low bin is inclusive, high bin is exclusive
  // i.e. for 0-10% use minCent=0, maxCent=4. 


  int minCent[12] = {0,2,4,6,8,10,12,14,16,20,24,28};
  int maxCent[12] = {2,4,6,8,10,12,14,16,20,24,28,32};
  
  // absolute values of eta are used, negative eta values are included
  // for full eta range set etaMin=0 etaMax=2.4
  Double_t  etaMin = 0.0;
  Double_t  etaMax = 0.4;

double dndeta[12]= {1612,1313,1082,893.9,731.6,596.8, 482.3,383.9, 266.5, 153.2,79.5, 36.3 };

  double ptbins[18]={0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,12.0};

  
  // Input File
  TFile *f = new TFile("LowPtSpectrum_fine_Full_d.root");
  TFile *fcorr[4] ;
  fcorr[0]= new TFile("validation3D_HydjetNoWedge_100k_flowSQ_goodTightdz10chi40_vertexZ10.root");
  fcorr[1]= new TFile("validation3D_HydjetNoWedge_100k_flowSQ_goodTightdz14chi80_vertexZ10.root");
  fcorr[2]= new TFile("validation3D_HydjetNoWedge_100k_flowSQ_vertexZ10.root");
  fcorr[3]= new TFile("validation3D_HydjetNoWedge_100k_flowSQ_goodTightdz6chi20_vertexZ10.root");
  // =========================

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);


  set_plot_style();

  char dirname[4][256] = {"spectrumGoodTightdz10chi40",
                          "spectrumGoodTightdz14chi80",
                          "spectrumGoodz14chi80",
                          "spectrumGoodTightdz6chi20"};
  TH1D * spec[10][12];
  TH1D * eff[10][12];
  TH1D * fak[10][12];
  TH1D * sim[10][12];
  TH1D * rec[10][12];

  TGraphErrors * gSpectrum[10][12];
  TF1 * f1[10][12];

  Double_t  meanpt[10][12];
  Double_t  meanpterr[10][12];

  
    Double_t px[100],py[100],pxe[100],pye[100];
  for( int c=0; c<12; c++)
  {
    for( int j=0; j<4; j++)
    {
      spec[j][c] = getSpectrum( f, Form("%s/tracks3D",dirname[j]), minCent[c], maxCent[c], etaMin, etaMax ); 
      eff[j][c] = getSpectrum( fcorr[j], "hitrkPixelEffAnalyzer/heff3D", minCent[c], maxCent[c], etaMin, etaMax );
      sim[j][c] = getSpectrum( fcorr[j], "hitrkPixelEffAnalyzer/hsim3D", minCent[c], maxCent[c], etaMin, etaMax );
      fak[j][c] = getSpectrum( fcorr[j], "hitrkPixelEffAnalyzer/hfak3D", minCent[c], maxCent[c], etaMin, etaMax );
      rec[j][c] = getSpectrum( fcorr[j], "hitrkPixelEffAnalyzer/hrec3D", minCent[c], maxCent[c], etaMin, etaMax );

      // determine correction factors


      TH1F * nevts = (TH1F *) f->Get(Form("%s/nevts",dirname[j]));
      double Nevt = nevts->Integral( nevts->GetXaxis()->FindBin(minCent[c]+0.001), 
                                     nevts->GetXaxis()->FindBin(maxCent[c]-0.001) );


      Double_t maxy = 0.;
      Double_t miny = 100000.;

    for( int i=0; i<=16;i++)
    {
      double ptmin = ptbins[i]; double ptmax = ptbins[i+1];

      double iptmin = spec[j][c]->FindBin(ptmin+1e-3);
      double iptmax = spec[j][c]->FindBin(ptmax-1e-3);
      double icptmin = eff[j][c]->FindBin(ptmin+1e-3);
      double icptmax = eff[j][c]->FindBin(ptmax-1e-3);
      double pt = 0.;
      for( int k=iptmin;k<=iptmax;k++) pt += spec[j][c]->GetBinCenter(k) * spec[j][c]->GetBinContent(k);
      pt /= spec[j][c]->Integral(iptmin,iptmax);

      double yielderr;
      double yield = spec[j][c]->IntegralAndError(iptmin,iptmax,yielderr);
      yield /= (ptmax-ptmin) * Nevt * 2 * (etaMax-etaMin);
      yielderr /= (ptmax-ptmin) * Nevt * 2 * (etaMax-etaMin);

      double efmin = eff[j][c]->GetBinContent(icptmin)/sim[j][c]->GetBinContent(icptmin);
      double efmax = eff[j][c]->GetBinContent(icptmax)/sim[j][c]->GetBinContent(icptmax);
      double famin = fak[j][c]->GetBinContent(icptmin)/rec[j][c]->GetBinContent(icptmin);
      double famax = fak[j][c]->GetBinContent(icptmax)/rec[j][c]->GetBinContent(icptmax);

      double ef = (pt-ptmin)*(efmax-efmin)/(ptmax-ptmin)+efmin;
      double fa = (pt-ptmin)*(famax-famin)/(ptmax-ptmin)+famin;
      yield *= (1-fa)/ef;

      px[i] = pt;
      py[i] = yield;
      if( py[i] > maxy ) maxy = py[i];
      if( py[i] < miny ) miny = py[i];
      pxe[i] = 0.;
      pye[i] = yielderr;
    }

      gSpectrum[j][c] = new TGraphErrors( 16, px, py, pxe, pye );
      gSpectrum[j][c]->SetMarkerStyle(20);
      f1[j][c] = new TF1(Form("f1_%d_%d",c,j),"6.28318531 * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.3,3.1);
     f1[j][c]->SetParameters(10.8327,2.221289,dndeta[c]);
//      f1[j][c]->FixParameter(2,dndeta[c]);
      f1[j][c]->SetLineColor(kBlue);
      f1[j][c]->SetLineWidth(1);
      f1[j][c]->SetLineStyle(1);
      gSpectrum[j][c]->Fit(Form("f1_%d_%d",c,j),"0RQ");
  //    TFitResultPtr r;
      TF1 * fitptmx = new TF1("fitptmx","6.28318531 * x * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.0,3.1);
      fitptmx->SetParameters(f1[j][c]->GetParameter(0),f1[j][c]->GetParameter(1), f1[j][c]->GetParameter(2) );
    double numerator = fitptmx->Integral(0.0,0.3);
    double denominator = f1[j][c]->Integral(0.0,0.3);
    for( int p = 0; p<gSpectrum[j][c]->GetN(); p++)
    {
       Double_t  ppx, ppy;
       gSpectrum[j][c]->GetPoint(p,ppx,ppy);
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
    }
    meanpt[j][c] = numerator / denominator;


      meanpterr[j][c] = 0.;
    }
    
  }

 
  TH1D* hDumCen = new TH1D("hDumCen",";Centrality (%);<p_{T}> Ratio",40,0.,80.); 
  hDumCen->GetXaxis()->CenterTitle(); hDumCen->GetYaxis()->SetTitleOffset(1.1);
  hDumCen->SetMaximum(1.25); hDumCen->SetMinimum(0.75);

  TCanvas *c1 = new TCanvas("c1","Mean Pt",600,600);
 
  c1->cd();
  hDumCen->Draw();

  TGraphErrors * gmean[4];
  
  int colors[4] = {0,kRed, kBlack, kBlue};

  for( int j=1; j<4; j++)
  {
    Double_t mx[12], mxe[12], my[12], mye[12];
    for(int i=0;i<12;i++) 
    { 
      mx[i] = (double)(maxCent[i]+minCent[i]) * 1.25;
      mxe[i] = 0; 
      my[i] = meanpt[j][i] / meanpt[0][i];
      mye[i] = my[i] * sqrt( meanpterr[j][i]*meanpterr[j][i]/meanpt[j][i]/meanpt[j][i] 
                           + meanpterr[0][i]*meanpterr[0][i]/meanpt[0][i]/meanpt[0][i] );
    
    }  

    gmean[j] = new TGraphErrors(12, mx, my, mxe, mye);
    gmean[j]->SetMarkerStyle(20);
    gmean[j]->SetMarkerColor(colors[j]);
    gmean[j]->SetLineColor(colors[j]);
    gmean[j]->Draw("p");
  
  }
  TLegend *legPt = new TLegend(0.45,0.7,0.89,0.89);
  legPt->SetFillColor(0); legPt->SetBorderSize(0);
  legPt->AddEntry(gmean[1],"GoodTight dz14 chi80","lp");
  legPt->AddEntry(gmean[2],"GoodMerged","lp");
  legPt->AddEntry(gmean[3],"GoodTight dz6 chi20","lp");

  legPt->Draw();

  TLatex *  tex = new TLatex(20,0.8,Form( "%4.1f < |#eta| < %4.1f",etaMin, etaMax ));
    tex->SetTextSize(0.035);
    tex->Draw();

  TH1D* hDumPt = new TH1D("hDumPt",";p_{T} [GeV/c];ratio 1/N_{evt} d^{2}N_{ch}/d#etadp_{T}",40,0.,6.); 
  hDumPt->SetMaximum(1.5);
  hDumPt->GetXaxis()->CenterTitle(); hDumPt->GetYaxis()->SetTitleOffset(1.1);

  hDumPt->SetMaximum(1.2); hDumPt->SetMinimum(0.8);

  TCanvas *c8 = new TCanvas("c8","Spectrum",1000,500);
  c8->Divide(3,1,0,0);

  TGraphErrors * gRatio[4][12]; 

  for( int q=0;q<3; q++)
  {
    c8->cd(q+1);
    int i;
    if (q==0) i=0; 
    if (q==1) i=5; 
    if (q==2) i=11; 
    hDumPt->Draw();
    for( int j=1; j<4; j++)
    { 
      Double_t rx[100],ry[100],rxe[100],rye[100];
      Double_t dx, nx, ny, dy, dye, nye; 
      for( int p=0;p<gSpectrum[0][i]->GetN();p++)
      {
        gSpectrum[0][i]->GetPoint(p,dx,dy);
        dye = gSpectrum[0][i]->GetErrorY(p);
        gSpectrum[j][i]->GetPoint(p,nx,ny);
        nye = gSpectrum[j][i]->GetErrorY(p);
        rx[p] = nx; rxe[p]=0;
        ry[p] = ny / dy;
        rye[p] = ry[p] * sqrt ( nye*nye/ny/ny + dye*dye/dy/dy);        
      }
      gRatio[j][i] = new TGraphErrors(gSpectrum[0][i]->GetN(), rx, ry, rxe, rye);
      gRatio[j][i]->SetMarkerStyle(20);
      gRatio[j][i]->SetMarkerColor(colors[j]);
      gRatio[j][i]->SetLineColor(colors[j]);
      gRatio[j][i]->Draw("p"); 
    tex = new TLatex(1,0.87,Form("%d-%d Cent.",(int)(minCent[i]*2.5),(int)(maxCent[i]*2.5)));
    tex->SetTextSize(0.045);
    tex->Draw();

    if( i==0) legPt->Draw();
  
    tex = new TLatex(1,0.85,Form("%4.1f < |#eta| < %4.1f",etaMin,etaMax));
    tex->SetTextSize(0.045);
    if( i==0) tex->Draw();
    }
  } 
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

Double_t ExpFuncScaled(Double_t * x, Double_t * par) {
  Double_t f;
  Double_t T = par[0];
  Double_t m = par[1];
  Double_t dNdy = par[2];

  Double_t mt = TMath::Sqrt( m*m + x[0]*x[0] );
  f = 1./(2.*TMath::Pi()*T*(T+m));

  f = f*dNdy;

  f = f*TMath::Exp(-(mt-m)/T);
   f*=2.*TMath::Pi()*x[0];

  return f;

}

Double_t ExpFuncME(Double_t * x, Double_t * par) {
  Double_t f;
  Double_t T = par[0];
  Double_t m = par[1];
  Double_t dNdy = par[2];

  Double_t mt = TMath::Sqrt( m*m + x[0]*x[0] );
  f = 1./(2.*TMath::Pi()*T*(T+m));

  f = f*dNdy;

  f = f*TMath::Exp(-(mt-m)/T);
   f*=2.*TMath::Pi()*x[0]*x[0];

  return f;

}

