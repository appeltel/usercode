void spectraCompareHIN10005() {

  // === General Settings ===

  // Note centrality is by 2.5% bins, low bin is inclusive, high bin is exclusive
  // i.e. for 0-10% use minCent=0, maxCent=4. 


  int minCent[6] = {0,2,4,12,20,28};
  int maxCent[6] = {2,4,12,20,28,36};
  
  // HIN-11-005
  //int minCent[12] = {0,2,4,6,8,10,12,14,16,18,20,24};
  //int maxCent[12] = {2,4,6,8,10,12,14,16,18,20,24,28};

  int showCent[3] = {0,1,2};  
  
  // absolute values of eta are used, negative eta values are included
  // for full eta range set etaMin=0 etaMax=2.4
  Double_t  etaMin = 0.0;
  Double_t  etaMax = 1.0;

  double ptbins[20]={0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,12.0,16.0,22.0};


double dndeta[12]= {1612,1313,1082,893.9,731.6,596.8, 482.3,383.9, 266.5, 153.2,79.5, 36.3 };


  
  // Input File
//  TFile *f = new TFile("LowPtSpectrum_fine_Full_b.root");
  TFile *f = new TFile("LowPtSpectrum_fine_Full_d.root");
  TFile *fcorr;
  fcorr= new TFile("validation3D_HydjetNoWedge_100k_flowSQ_goodTightdz10chi40_vertexZ10.root");
  //fcorr= new TFile("validation3D_HydjetNoWedge_100k_flowSQ_vertexZ10.root");
//  TFile *fsteve = new TFile("EPSpectra_gt10dz40chi_13_V2.root");
  //  TFile *fsteve = new TFile("EPSpectra_gt10dz40chi_13_V2_HIN11005.root");
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
  TGraphErrors * gSpectrumSteve[12];
  TF1 * f1[12];
  TSpline3 * spline[12];
  
  Double_t px[100],py[100],pxe[100],pye[100];
  for( int c=0; c<6; c++)
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


      Double_t maxy = 0.;
      Double_t miny = 100000.;
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
      pye[i] = sqrt(yielderr*yielderr + py[i]*py[i]*0.06*0.06);
      if( c==5) cout << pt << "  " << yield << "  " << yielderr << "  "  << pye[i] << endl;
    }

      gSpectrum[c] = new TGraphErrors( 19, px, py, pxe, pye );
      gSpectrum[c]->SetMarkerStyle(20);
      f1[c] = new TF1(Form("f1_%d",c),PowerFuncScaled,0.8,5,3);
      f1[c]->SetParNames("n","p0","dN/dy");
      f1[c]->SetParameters(28,5,dndeta[c]);
      //f1[c]->FixParameter(2,dndeta[c]);
      f1[c]->SetLineColor(kBlue);
      f1[c]->SetLineWidth(1);
      f1[c]->SetLineStyle(1);

      spline[c] = new TSpline3(Form("grs_%d",c),gSpectrum[c]);
      spline[c]->SetLineColor(kBlue);
      spline[c]->SetLineWidth(1);
      spline[c]->SetLineStyle(1);

      int nnp=0;
      ifstream infile;
      infile.open(Form("spectrumHIN10005_%d_%d.txt",(int)(minCent[c]*2.5),(int)(maxCent[c]*2.5)));
      if (infile.is_open())
      {
      std::string line;
      while(std::getline(infile,line,'\n'))
      {
        double v1,v2,v3,v4,v5,v6,;
        int res = sscanf(line.c_str(),"%lg %lg %lg %lg %lg %lg",&v1,&v2,&v3,&v4,&v5,&v6);
        if (res!=6) continue;
      //  cout <<"v1="<<v1<<" v2="<<v2<<" v3="<<v3 << endl;
        px[nnp]=v1;
        pxe[nnp]=0.0;
        py[nnp]=v4;
        pye[nnp]=v5;

        if (c==1) py[nnp] /= 1e3;
        if (c==2) py[nnp] /= 1e6;
        if (c==3) py[nnp] /= 1e9;
        if (c==4) py[nnp] /= 1e12;
        if (c==5) py[nnp] /= 1e15;
        if (c==1) pye[nnp] /= 1e3;
        if (c==2) pye[nnp] /= 1e6;
        if (c==3) pye[nnp] /= 1e9;
        if (c==4) pye[nnp] /= 1e12;
        if (c==5) pye[nnp] /= 1e15;

        py[nnp]=py[nnp]*2*3.14159265*px[nnp];
        pye[nnp]=pye[nnp]*2*3.14159265*px[nnp];
        nnp++;
      }
      infile.close();
      }


    gSpectrumSteve[c] = new TGraphErrors(nnp, px,py,pxe, pye);
    gSpectrumSteve[c]->SetMarkerStyle(24);
    gSpectrumSteve[c]->SetMarkerColor(kBlue);
    gSpectrumSteve[c]->SetLineColor(kBlue);
    gSpectrum[c]->Fit(Form("f1_%d",c),"0RNQS");
    cout << "Centrality " << (int)(minCent[c]*2.5) << "-" <<  (int)(maxCent[c]*2.5) << endl << endl;
    cout << "pT    yield     error" << endl;
    for(int d = 0; d<gSpectrumSteve[c]->GetN(); d++) 
    {  
      Double_t dx, dy, dye;
      gSpectrumSteve[c]->GetPoint(d,dx,dy);
      dye = gSpectrumSteve[c]->GetErrorY(d);
      cout << dx << " " << dy << " " << dye << endl;
    }

    cout << endl << endl;   
  }

 

  
  

  TLegend *legPt = new TLegend(0.45,0.7,0.89,0.89);
  legPt->SetFillColor(0); legPt->SetBorderSize(0);
  legPt->AddEntry(gSpectrum[0],"HIN-10-002","lp");
  legPt->AddEntry(gSpectrumSteve[0],"HIN-10-005","lp");
  //legPt->AddEntry(spline[0],"Cubic Spline Fit","lp");
  legPt->AddEntry(f1[0],"Power Law Fit","lp");


  TH1D* hDumPt = new TH1D("hDumPt",";p_{T} [GeV/c];ratio 1/N_{evt} d^{2}N_{ch}/d#etadp_{T}",40,0.8,5.1); 
  hDumPt->GetXaxis()->CenterTitle(); hDumPt->GetYaxis()->SetTitleOffset(1.1);

  TH1D * hDumPtRat = hDumPt->Clone("hDumPtRat");
  hDumPtRat->SetMaximum(1.5); hDumPtRat->SetMinimum(0.5);
  hDumPt->SetMaximum(400.); hDumPt->SetMinimum(0.01);

  TCanvas *c8 = new TCanvas("c8","Spectrum",1000,800);
  c8->Divide(3,2,0,0);
  TGraphErrors * gRatio[4][12]; 

  for( int q=0;q<3; q++)
  {
    c8->cd(q+4);
    //c8->GetPad(q+4)->SetLogx();
    int i=showCent[q];
    hDumPtRat->Draw();
    for( int j=0; j<2; j++)
    { 
      Double_t rx[100],ry[100],rxe[100],rye[100];
      Double_t dx, nx, ny, dy, dye, nye;
      int np;
      if (j==0) np = gSpectrum[i]->GetN();
      if (j==1) np = gSpectrumSteve[i]->GetN();  
      for( int p=0;p<np;p++)
      {
        if( j==0) gSpectrum[i]->GetPoint(p,dx,dy);
        if( j==1) gSpectrumSteve[i]->GetPoint(p,dx,dy);
        if(j==0)dye = gSpectrum[i]->GetErrorY(p);
        if(j==1)dye = gSpectrumSteve[i]->GetErrorY(p);
        rx[p] = dx; rxe[p]=0;
        
        ry[p] = dy / f1[i]->Eval(dx);
        //ry[p] = dy / spline[i]->Eval(dx);
        rye[p] = dy * 0.06 / f1[i]->Eval(dx);
      }

      gRatio[j][i] = new TGraphErrors(np, rx, ry, rxe, rye);
      gRatio[j][i] = new TGraphErrors(np, rx, ry, rxe, rye);
      if(j==0)
      { 
        gRatio[j][i]->SetMarkerStyle(20);
        gRatio[j][i]->SetMarkerColor(kBlack);
        gRatio[j][i]->SetLineColor(kBlack);
      }  
      if(j==1)
      { 
        gRatio[j][i]->SetMarkerStyle(24);
        gRatio[j][i]->SetMarkerColor(kBlue);
        gRatio[j][i]->SetLineColor(kBlue);
      }
      gRatio[j][i]->Draw("p"); 


    }
    tex = new TLatex(1.3,0.7,Form("%d-%d Cent.",(int)(minCent[i]*2.5),(int)(maxCent[i]*2.5)));
    tex->SetTextSize(0.045);
    tex->Draw();
  
    tex = new TLatex(1.3,0.8,Form("%4.1f < |#eta| < %4.1f",etaMin,etaMax));
    tex->SetTextSize(0.045);
    if( i==3 || i==0) tex->Draw();

    tex = new TLatex(1.3,0.6,"97% Trigger Efficiency" );
    tex->SetTextSize(0.045);
    if( i==3 || i==0) tex->Draw();

    c8->cd(q+1);
    //c8->GetPad(q+1)->SetLogx();
    c8->GetPad(q+1)->SetLogy();
    hDumPt->Draw();

    if ( i==3 || i==0) legPt->Draw();
    f1[i]->Draw("same");
    //spline[i]->Draw("same");
    gSpectrum[i]->Draw("p");
    gSpectrumSteve[i]->Draw("p");
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

