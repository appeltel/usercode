void meanPtSysFit7() {

  // === General Settings ===

  // Note centrality is by 2.5% bins, low bin is inclusive, high bin is exclusive
  // i.e. for 0-10% use minCent=0, maxCent=4. 


  int minCent[12] = {0,2,4,6,8,10,12,14,16,20,24,28};
  int maxCent[12] = {2,4,6,8,10,12,14,16,20,24,28,32};
  
  // absolute values of eta are used, negative eta values are included
  // for full eta range set etaMin=0 etaMax=2.4


  Double_t  etaMin = 0.0; 
  Double_t  etaMax = 0.4;
  char cfn[12][256];
  for( int i=0;i<12;i++) strncpy(cfn[i],Form("./finaldata/Spectra_00_04_%d.txt",i),256); 

  double spectraSysErr;
  if( etaMax < 0.9 ) spectraSysErr = 0.06;
  else if( etaMax < 1.3 ) spectraSysErr = 0.09;
  else if( etaMax < 2.5 ) spectraSysErr = 0.14;

double dndeta[12]= {1612,1313,1082,893.9,731.6,596.8, 482.3,383.9, 266.5, 153.2,79.5, 36.3 };

//  double ptbins[18]={0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,12.0};
  double ptbins[19]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,12.0};

  
  // Input File
  TFile *f = new TFile("LowPtSpectrum_fine_Full_d.root");
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
  TGraphErrors * gSpectrumSys[12];
  TF1 * f1[8][12];

  Double_t  meanpt[8][12];
  Double_t  meanpterr[8][12];


  double goodnessMin[10];
  double goodnessMax[10];
  for( int i=0;i<10;i++)
  {
    goodnessMin[i]=1e20;  
    goodnessMax[i]=0.;
  }  
    Double_t px[100],py[100],pxe[100],pye[100], pyse[100];
  for( int c=0; c<12; c++)
  {

    ifstream infile;
    infile.open(cfn[c]);
    if (infile.is_open())
    {
      cout << "reading " << cfn[c] << endl;
      std::string line;
      int np=0;
      while(std::getline(infile,line,'\n'))
      {
        double v1,v2,v3,v4,v5,v6,v7;
        int res = sscanf(line.c_str(),"%lg %lg %lg",&v1,&v2,&v3);
        if (res!=3) continue;
        if (v1 < 0.2 ) continue;
      //  cout <<"v1="<<v1<<" v2="<<v2<<" v3="<<v3 << endl;
        px[np]=v1;
        pxe[np]=0.0;
        py[np]=v2;
        pye[np]=v3;

        pyse[np]=0.;


        np++;
      }
      infile.close();
    }


    gSpectrum[c] = new TGraphErrors( 16, px, py, pxe, pye );
    gSpectrumSys[c] = new TGraphErrors( 16, px, py, pxe, pyse );
    gSpectrum[c]->SetMarkerStyle(20);
    for( int j=0; j<5; j++)
    {
      if( j==0) f1[j][c] = new TF1(Form("f1_%d_%d",c,j),"6.28318531 * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.3,3.1);
      if( j==1) f1[j][c] = new TF1(Form("f1_%d_%d",c,j),"6.28318531 * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.2,3.1);
      if( j==2) f1[j][c] = new TF1(Form("f1_%d_%d",c,j),"6.28318531 * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.4,3.1);
      if ( j<3) f1[j][c]->SetParameters(10.8327,2.221289,dndeta[c]);
      if ( j==2) f1[j][c]->SetParameters(10.8327,2.221289,dndeta[c]/3.);
      if (j==3) f1[j][c] = new TF1(Form("f1_%d_%d",c,j),ExpFuncScaled,0.3,0.7,3);
      if (j==3) f1[j][c]->SetParameters(0.2,0.,dndeta[c]);
      if (j==4 ) f1[j][c] = new TF1(Form("f1_%d_%d",c,j),"6.28318531 * x * [2]*pow((1+x/[1]),-[0])",0.3,5.1);
      if (j==4 )  f1[j][c]->SetParameters(10.8327,2.221289,1201);

      f1[j][c]->SetLineColor(kBlue);
      f1[j][c]->SetLineWidth(1);
      f1[j][c]->SetLineStyle(1);

      gSpectrum[c]->Fit(Form("f1_%d_%d",c,j),"0RQ");

      meanpterr[j][c] = 0.;
         
      TF1 * f1x;
      if(j==3) f1x = new TF1(Form("f1x_%d_%d",c,j),ExpFuncME,0.0,3.1,3);
      if(j==4) f1x = new TF1(Form("f1x_%d_%d",c,j),"6.28318531 *x * x * [2]*pow((1+x/[1]),-[0])",0.0,3.1);
      if(j<3) f1x = new TF1(Form("f1x_%d_%d",c,j),
            "6.28318531 * x * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.0,3.1);
      f1x->SetParameters(f1[j][c]->GetParameter(0),f1[j][c]->GetParameter(1), f1[j][c]->GetParameter(2) );
      double numerator = f1x->Integral(0.0,0.3);
      if( j==1) numerator = f1x->Integral(0.0,0.2);
      if( j==2) numerator = f1x->Integral(0.0,0.4);
      double denominator = f1[j][c]->Integral(0.0,0.3);
      if( j==1) denominator = f1[j][c]->Integral(0.0,0.2);
      if( j==2) denominator = f1[j][c]->Integral(0.0,0.4);
      for( int p = 0; p<gSpectrum[c]->GetN(); p++)
      {
         Double_t  ppx, ppy;
         gSpectrum[c]->GetPoint(p,ppx,ppy);
         if ( 0.2 < ppx && ppx < 0.3 && j==1 ) numerator += ppx * ppy * 0.1;
         if ( 0.3 < ppx && ppx < 0.6 && j!=2 ) numerator += ppx * ppy * 0.1;
         if ( 0.4 < ppx && ppx < 0.6 && j==2 ) numerator += ppx * ppy * 0.1;
         if ( 0.6 < ppx && ppx < 1.2 ) numerator += ppx * ppy * 0.2;
         if ( 1.2 < ppx && ppx < 2.0 ) numerator += ppx * ppy * 0.4;
         if ( 2.0 < ppx && ppx < 4.0 ) numerator += ppx * ppy * 0.5;
         if ( 4.0 < ppx && ppx < 6.0 ) numerator += ppx * ppy * 1.0;
         if ( 0.2 < ppx && ppx < 0.3 && j==1 ) denominator += ppy * 0.1;
         if ( 0.3 < ppx && ppx < 0.6 && j!=2 ) denominator += ppy * 0.1;
         if ( 0.4 < ppx && ppx < 0.6 && j==2 ) denominator += ppy * 0.1;
         if ( 0.6 < ppx && ppx < 1.2 ) denominator += ppy * 0.2;
         if ( 1.2 < ppx && ppx < 2.0 ) denominator += ppy * 0.4;
         if ( 2.0 < ppx && ppx < 4.0 ) denominator += ppy * 0.5;
         if ( 4.0 < ppx && ppx < 6.0 ) denominator += ppy * 1.0;
      }
      meanpt[j][c] = numerator / denominator;
      if ( j==0 ) cout << c << "  " << meanpt[0][c] << endl;

    }
    
  }

 
  TH1D* hDumCen = new TH1D("hDumCen",";Centrality (%);<p_{T}> Ratio",40,0.,80.); hDumCen->SetMaximum(2.0);
  hDumCen->GetXaxis()->CenterTitle(); hDumCen->GetYaxis()->SetTitleOffset(1.1);

  TCanvas *c1 = new TCanvas("c1","Mean Pt",600,600);
 
  c1->cd();
  hDumCen->Draw();

  TGraphErrors * gmean[8];
  
  int colors[7] = {kViolet,kRed, kBlack, kBlue, kGreen+1, kOrange+1, kGray};

  for( int j=0; j<5; j++)
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
  legPt->AddEntry(gmean[0],"Fit A pT > 0.3","lp");
  legPt->AddEntry(gmean[1],"Fit A pT > 0.2","lp");
  legPt->AddEntry(gmean[2],"Fit A pT > 0.4","lp");
  legPt->AddEntry(gmean[3],"Fit B pT > 0.3","lp");
  legPt->AddEntry(gmean[4],"Fit C pT > 0.3","lp");
  //legPt->AddEntry(gmean[5],"Fit C pT > 0.3 fixed dN/d#eta","lp");

  legPt->Draw();

  TLatex *  tex = new TLatex(0.2,0.3,Form("%4.1f < |#eta| < %4.1f",etaMin,etaMax));
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

