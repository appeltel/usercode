void analyzeVNfake() {

  // === General Settings ===
  //
  // harmonic is v_n{psi_n} for n < 10. Above 10 are special values:
  //
  // 20 - v_2{Generator RP}
  // 42 - v_4{psi_2}
  // 62 - v_6{psi_2}
  // 420 - v_4{Generator RP}
  // 620 - v_6{Generator RP}
  int harmonic = 420;

  // Note centrality is by 2.5% bins, low bin is inclusive, high bin is exclusive
  // i.e. for 0-10% use minCent=0, maxCent=4. 
  int minCent = 0;
  int maxCent = 4;
  
  // absolute values of eta are used, negative eta values are included
  // for full eta range set etaMin=0 etaMax=2.4
  Double_t  etaMin = 0.0;
  Double_t  etaMax = 2.4;

  // maximum pT to be shown in some plots
  // integration range is hard coded as 0.3 < pT < 3.0
  Double_t ptMax = 2.0;
  
  // Input File
  TFile *f = new TFile("validation3D_HYDJETNoWedge_33k_flow40_goodTightdz18chi120_vertexZ10_2ep.root");

  // Additional files for comparing cuts, and switch to turn the comparison on or off
  bool doCuts = true;
  TFile * fCut1;
  TFile * fCut2;
  if (doCuts) fCut1 = new TFile("validation3D_AMPTNoWedge_33k_flow40_goodTightdz14chi80_vertexZ10_2ep.root");
  if (doCuts) fCut2 = new TFile("validation3D_AMPTNoWedge_33k_flow40_goodTightdz10chi40_vertexZ10_2ep.root");
  int dzcuts[3] = { 18, 14, 10 };
  int chi2cuts[3] = { 120, 80, 40 };

  // Do a linear fit on the 4-panel centrality comparison
  bool doLinear = false;

  // =========================

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);

  int planeharm = harmonic;
  if( harmonic > 10 ) planeharm = 2;

  set_plot_style();


  // get spectra

  // TH1D * hSimPt = getSpectrum( f, "sim", minCent, maxCent, etaMin, etaMax);
  // TH1D * hEffPt = getSpectrum( f, "reff", minCent, maxCent, etaMin, etaMax);
  // TH1D * hRecPt = getSpectrum( f, "rec", minCent, maxCent, etaMin, etaMax);
  // TH1D * hFakPt = getSpectrum( f, "fak", minCent, maxCent, etaMin, etaMax);
   
  // Calculate Event Plane Resolution

  TH2F * epres[5];
  TH1D * cepres[5];
  Double_t eventplaneres[5];
  Double_t eventplanereserr[5];
  for( int i=0;i<5;i++)
  {
    epres[i] = (TH2F*) f->Get(Form("hitrkPixelEffAnalyzer/epres_%d",i));
    Int_t cmax = epres[0]->GetXaxis()->FindBin(maxCent-0.01); 
    Int_t cmin = epres[0]->GetXaxis()->FindBin(minCent+0.01); 
    cepres[i] = epres[i]->ProjectionY(Form("cepres_%d",i),cmin,cmax,"e");
    eventplaneres[i] = cepres[i]->GetMean()>0. ? sqrt(cepres[i]->GetMean()) : 0.0001; 
    cout << "Event Plane Res " << i+2 << "  " << eventplaneres[i] << endl; 
    eventplanereserr[i] = 0.5 * cepres[i]->GetRMS() / sqrt(cepres[i]->GetEntries()) / eventplaneres[i] ;
  }

  double er = eventplaneres[planeharm-2];
  double ers = eventplanereserr[planeharm-2];

  // Make v_n(pT) TGraphs (no resolution correction for the vno ones)

  TGraphErrors * gSim_vno = get_vno( f, "sim", minCent, maxCent, etaMin, etaMax, harmonic );
  TGraphErrors * gEff_vno = get_vno( f, "reff", minCent, maxCent, etaMin, etaMax, harmonic );
  TGraphErrors * gFak_vno = get_vno( f, "fak", minCent, maxCent, etaMin, etaMax, harmonic );
  TGraphErrors * gRec_vno = get_vno( f, "rec", minCent, maxCent, etaMin, etaMax, harmonic );

  TGraphErrors * gSim_vn = get_vn( f, "sim", minCent, maxCent, etaMin, etaMax, harmonic, er, ers );
  TGraphErrors * gEff_vn = get_vn( f, "reff", minCent, maxCent, etaMin, etaMax, harmonic, er, ers );
  TGraphErrors * gFak_vn = get_vn( f, "fak", minCent, maxCent, etaMin, etaMax, harmonic, er, ers );
  TGraphErrors * gRec_vn = get_vn( f, "rec", minCent, maxCent, etaMin, etaMax, harmonic, er, ers );


  // Plot EP resolution for all harmonics

  TH1D* hDumER = new TH1D("hDumER",Form(";Centrality;R",harmonic),40,0.0,100.); hDumER->SetMaximum(1.0);
  TGraphErrors * gepres[5];
  for( int h=0;h<5;h++) 
  {
    Double_t ypt[10];
    Double_t xpt[10] = {5.,15.,25.,35.,45.,55.,65.,75.,85.,95.};
    Double_t xpterr[10] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Double_t ypterr[10];
    for( int i =0; i<10;i++)
    {
      Int_t cmax = epres[h]->GetXaxis()->FindBin((double)(i+1) * 4 - 0.001 );
      Int_t cmin = epres[h]->GetXaxis()->FindBin((double)(i) * 4 + 0.001 );
      TH1D * cep = epres[h]->ProjectionY("cep",cmin,cmax,"e");
      ypt[i] = cep->GetMean()>0. ? sqrt(cep->GetMean()) : - sqrt( - cep->GetMean());
      ypterr[i] = 0.5 * cep->GetRMS() / sqrt(cep->GetEntries()) /  fabs( ypt[i] );
    }
    gepres[h] = new TGraphErrors(10,xpt,ypt,xpterr,ypterr);
  }

  TCanvas *c5 = new TCanvas("c5","Resolution",500,500);
  c5->cd(); hDumER->Draw();

  int rescolor[5] = {kBlack, kRed, kBlue, kGreen+1, kViolet};
  for( int i=0;i<5;i++)
  {
    gepres[i]->SetMarkerStyle(20);
    gepres[i]->SetMarkerSize(1.1);
    gepres[i]->SetMarkerColor(rescolor[i]);
    gepres[i]->SetLineColor(rescolor[i]);
    gepres[i]->Draw("p");
  }
  
  TLegend *legER = new TLegend(0.45,0.7,0.89,0.89);
  legER->SetFillColor(0); legER->SetBorderSize(0);
  for( int i=0;i<5;i++)
  {
    legER->AddEntry(gepres[i],Form("#Psi_{%d}",i+2),"lp");
  } 
  legER->Draw();


  // Now plot the corrected v_n signals

  TH1D* hDumPt = new TH1D("hDumPt",Form(";p_{T} [GeV/c];v_{%d}",harmonic),40,0.1,ptMax); hDumPt->SetMaximum(1.0);
  hDumPt->GetXaxis()->CenterTitle(); hDumPt->GetYaxis()->SetTitleOffset(1.1);
  TH1D* hDumPt2 = (TH1D*) hDumPt->Clone("hDumPt2"); hDumPt2->SetMaximum(0.1); hDumPt2->SetMinimum(-0.02); 
  TH1D* hDumPt3 = (TH1D*) hDumPt->Clone("hDumPt3"); hDumPt3->SetMaximum(0.00049); 
  TH1D* hDumPt4 = new TH1D("hDumPt4",Form(";p_{T} [GeV/c];v_{%d}",harmonic),40,0.1,1.0); hDumPt->SetMaximum(1.0);
  hDumPt4->GetXaxis()->CenterTitle(); hDumPt4->GetYaxis()->SetTitleOffset(1.1);


  gSim_vn->SetMarkerStyle(20);
  gSim_vn->SetMarkerSize(0.8);
  gSim_vn->SetMarkerColor(kBlue);
  gSim_vn->SetLineColor(kBlue);
  
  gEff_vn->SetMarkerStyle(20);
  gEff_vn->SetMarkerSize(0.6);
  gEff_vn->SetMarkerColor(kGreen);
  gEff_vn->SetLineColor(kGreen);
  
  gRec_vn->SetMarkerStyle(20);
  gRec_vn->SetMarkerSize(1.0);
  gRec_vn->SetMarkerColor(kBlack);
  gRec_vn->SetLineColor(kBlack);
  
  gFak_vn->SetMarkerStyle(20);
  gFak_vn->SetMarkerSize(0.9);
  gFak_vn->SetMarkerColor(kRed);
  gFak_vn->SetLineColor(kRed);
  

  TLegend *legEta = new TLegend(0.35,0.3,0.65,0.5);
  legEta->SetFillColor(0); legEta->SetBorderSize(0);

  TLegend *legPt = new TLegend(0.45,0.7,0.89,0.89);
  legPt->SetFillColor(0); legPt->SetBorderSize(0);
  legPt->AddEntry(gSim_vn,"Generated Particles","lp");
  legPt->AddEntry(gRec_vn,"Reconstructed Tracks","lp");
  legPt->AddEntry(gFak_vn,"Fake Tracks","lp");
  legPt->AddEntry(gEff_vn,"Matched Tracks","lp");

  TCanvas *c8 = new TCanvas("c8","Flow",600,600);
  hDumPtEff=(TH1F*) hDumPt2->Clone("hDumPtEff");
  hDumPtEff->SetTitle(Form("v%d coefficient %4.1f < |#eta| < %4.1f, %4.0f -%4.0f Pct Centrality",
                          harmonic, 
                          etaMin,etaMax,
                           minCent*2.5,maxCent*2.5));
  c8->cd(); hDumPtEff->Draw(); 

  gFak_vn->Draw("p");
  gRec_vn->Draw("p");
  gSim_vn->Draw("p");  
  gEff_vn->Draw("p");

  legPt->Draw();


  // Draw ratio plot of reco to matched, and fit

  TCanvas *c2 = new TCanvas("c2","Ratio Reco/Matched",600,600);
  c2->cd();
  hDumPtRat=(TH1F*) hDumPt2->Clone("hDumPtRat");
  hDumPtRat->SetTitle(Form("v%d Reco / v%d Matched %4.1f < |#eta| < %4.1f, %4.0f -%4.0f Pct Centrality;p_{T} [GeV/c];v_{%d} ratio",
                      harmonic, harmonic,
                      etaMin,etaMax,
                      minCent*2.5,maxCent*2.5,harmonic));
  hDumPtRat->SetMaximum(1.5);
  hDumPtRat->SetMinimum(0.5);
  hDumPtRat->Draw();

  Double_t xc2[100], xec2[100], yc2[100], yec2[100];

  for( int p=0; p<gEff_vno->GetN(); p++)
  {
    Double_t ay, by, aye, bye;
    gRec_vno->GetPoint(p,xc2[p],ay);
    gEff_vno->GetPoint(p,xc2[p],by);
    aye = gRec_vno->GetErrorY(p);
    bye = gEff_vno->GetErrorY(p);
    xec2[p] = 0.0;
    yc2[p] = ay / by;
    yec2[p] = fabs(yc2[p]) * sqrt( aye*aye/ay/ay + bye*bye/by/by );
  }

  TGraphErrors * gRat = new TGraphErrors( gEff_vno->GetN(), xc2, yc2, xec2, yec2 );

  //for( int i=0;i<gRatFit->GetN();i++) gRatFit->SetPointError(i,0.0,0.0);
  gRat->SetMarkerStyle(20);
  gRat->SetMarkerSize(1.1);
  gRat->Draw("p");

  TF1 *f2 = new TF1("f2","[0]/x + [1]/(x*x) + 1.0",0.3,3.0);
  f2->SetLineColor(kBlue);
  f2->SetLineWidth(1);
  f2->SetLineStyle(1);
  gRat->Fit("f2","0RROB");
  f2->Draw("same"); 


  // Now plot ratio of fake v_n(p_T) to integral v_n

  TCanvas *c3 = new TCanvas("c3","Fake to Integral Ratio",600,600);
  c3->cd();
  hDumPtDif=(TH1F*) hDumPt4->Clone("hDumPtDif");
  hDumPtDif->SetTitle(Form("v%d Fake / v%d Integral %4.1f < |#eta| < %4.1f, %4.0f -%4.0f Pct Centrality;p_{T} [GeV/c];v_{%d} ratio",
                      harmonic, harmonic,
                      etaMin,etaMax,
                      minCent*2.5,maxCent*2.5,harmonic));
  hDumPtDif->SetMaximum(2.5);
  hDumPtDif->SetMinimum(0.0);
  hDumPtDif->Draw();

  double vnInt;

  TGraphErrors * gRat2 = get_vnFakRat( f, minCent, maxCent, etaMin, etaMax, harmonic, vnInt );


  gRat2->SetMarkerStyle(20);
  gRat2->SetMarkerSize(1.1);

  TLatex  * tex;
    tex = new TLatex(0.2,0.3,Form("<v_{n}^{EP,reco}> = %f",vnInt));
    tex->SetTextSize(0.035);
    tex->Draw();

  TF1 *f3 = new TF1("f3","pol0",0.2,1.0);
  f3->SetLineColor(kBlue);
  f3->SetLineWidth(1);
  f3->SetLineStyle(1);
  gRat2->Fit("f3","0R");
//  f3->Draw("same"); 
  
  double midfit = f3->GetParameter(0);
  double upfit = midfit + f3->GetParError(0); 
  double downfit = midfit - f3->GetParError(0); 
  TLine * lmean = new TLine(0.1,midfit,1.,midfit);
  lmean->SetLineColor(kBlue);
  lmean->Draw();
  TLine * lup = new TLine(0.1,upfit,1.,upfit);
  lup->SetLineColor(kGreen);
  lup->Draw();
  TLine * ldown = new TLine(0.1,downfit,1.,downfit);
  ldown->SetLineColor(kGreen);
  ldown->Draw();


  gRat2->Draw("p");


  // ==== 9-panel comparison plot ====

  TCanvas *c4 = new TCanvas("c4","Fake to Integral Ratio Cent and Eta Comparison",900,900);
  c4->Divide( 3,3,0.,0.);
  TH1D * hDumComp = hDumPtDif->Clone("hDumComp");
  hDumComp->SetTitle("");

  TGraphErrors * gRatComp[3][3];
  TLine * lupComp[3][3], * lmeanComp[3][3], * ldownComp[3][3];
  TF1 * fcomp[3][3];
  double vnDum;
  for( int c=0; c<3; c++)
  {
    for( int e=0;e<3;e++)
    {
      c4->cd(c*3 + e + 1);
      hDumComp->Draw();
      gRatComp[c][e] = get_vnFakRat( f, c*4, (c+1)*4, 0.8*(double)e, 0.8*(double)(e+1), harmonic, vnDum );
      fcomp[c][e] = new TF1(Form("fcomp%d%d",c,e),"pol0",0.2,1.0);
      gRatComp[c][e]->Fit(Form("fcomp%d%d",c,e),"0R");
      double mid = fcomp[c][e]->GetParameter(0);
      double up = mid + fcomp[c][e]->GetParError(0);
      double down = mid - fcomp[c][e]->GetParError(0);
      lmeanComp[c][e] = new TLine(0.1,mid,1.,mid);
      lmeanComp[c][e]->SetLineColor(kBlue);
      lmeanComp[c][e]->SetLineWidth(2);
      lupComp[c][e] = new TLine(0.1,up,1.,up);
      lupComp[c][e]->SetLineColor(kGreen);
      lupComp[c][e]->SetLineWidth(2);
      ldownComp[c][e] = new TLine(0.1,down,1.,down);
      ldownComp[c][e]->SetLineColor(kGreen);
      ldownComp[c][e]->SetLineWidth(2);
      lupComp[c][e]->Draw();
      lmeanComp[c][e]->Draw();
      ldownComp[c][e]->Draw();

      gRatComp[c][e]->SetMarkerStyle(20);
      gRatComp[c][e]->SetMarkerSize(1.2);
      gRatComp[c][e]->Draw("p");

      tex = new TLatex(0.2,0.4,Form("%d - %d Centrality",10*c, 10*(c+1) ));
      tex->SetTextSize(0.055);
      tex->Draw();
      tex = new TLatex(0.2,0.2,Form("%4.1f < |#eta| < %4.1f",0.8*(double)e, 0.8*(double)(e+1) ));
      tex->SetTextSize(0.055);
      tex->Draw();
      tex = new TLatex(0.6,2.0,Form("SF = %4.3f #pm %4.3f",mid, up-mid));
      tex->SetTextSize(0.055);
      tex->Draw();



    }
  }


  // ==== 3-panel cut comparison plot ====

  if ( doCuts)
  {
    TCanvas *c5 = new TCanvas("c5","Fake to Integral Ratio Cut Comparison",900,400);
    c5->Divide( 3,0,0.,0.);

    TGraphErrors * gRatCut[3];
    TLine * lupCut[3], * lmeanCut[3], * ldownCut[3];
    TF1 * fcut[3];
    TFile * ftmp[3];
    ftmp[0] = f; ftmp[1] = fCut1; ftmp[2] = fCut2;
    double vnCut[3];
    for( int c=0; c<3; c++)
    {
      c5->cd(c + 1);
      hDumComp->Draw();
      gRatCut[c] = get_vnFakRat( ftmp[c], minCent, maxCent, etaMin, etaMax, harmonic, vnCut[c] );
      fcut[c] = new TF1(Form("fcut%d",c),"pol0",0.2,1.0);
      gRatCut[c]->Fit(Form("fcut%d",c),"0R");
      double mid = fcut[c]->GetParameter(0);
      double up = mid + fcut[c]->GetParError(0);
      double down = mid - fcut[c]->GetParError(0);
      lmeanCut[c] = new TLine(0.1,mid,1.,mid);
      lmeanCut[c]->SetLineColor(kBlue);
      lmeanCut[c]->SetLineWidth(2);
      lupCut[c] = new TLine(0.1,up,1.,up);
      lupCut[c]->SetLineColor(kGreen);
      lupCut[c]->SetLineWidth(2);
      ldownCut[c] = new TLine(0.1,down,1.,down);
      ldownCut[c]->SetLineColor(kGreen);
      ldownCut[c]->SetLineWidth(2);
      lupCut[c]->Draw();
      lmeanCut[c]->Draw();
      ldownCut[c]->Draw();

      gRatCut[c]->SetMarkerStyle(20);
      gRatCut[c]->SetMarkerSize(1.2);
      gRatCut[c]->Draw("p");
      
      if( c==0)
      {
        tex = new TLatex(0.2,0.4,Form("%3.1f - %3.1f Centrality",minCent*2.5, maxCent*2.5 ));
        tex->SetTextSize(0.045);
        tex->Draw();
        tex = new TLatex(0.2,0.2,Form("%4.1f < |#eta| < %4.1f",etaMin, etaMax ));
        tex->SetTextSize(0.045);
        tex->Draw();
      }
      tex = new TLatex(0.65,0.4,Form("dz < %d",dzcuts[c] ));
      tex->SetTextSize(0.045);
      tex->Draw();
      tex = new TLatex(0.65,0.2,Form("chi2 < %d",chi2cuts[c] ));
      tex->SetTextSize(0.045);
      tex->Draw();
      tex = new TLatex(0.6,2.2,Form("<v_{n}^{EP,reco}> = %6.5f",vnCut[c]));
      tex->SetTextSize(0.045);
      tex->Draw();
      tex = new TLatex(0.6,2.0,Form("SF = %4.3f #pm %4.3f",mid, up-mid));
      tex->SetTextSize(0.045);
      tex->Draw();



    }
  }

  // ==== 4-panel centrality comparison plot ====

    TCanvas *c6 = new TCanvas("c6","Fake to Integral Ratio Centrality Comparison",1200,400);
    c6->Divide( 4,0,0.,0.);

    TGraphErrors * gRatCen[4];
    TLine * lupCen[3], * lmeanCen[4], * ldownCen[4];
    TF1 * fCen[4];
    TF1 * fLin[4];
    double vnCen[4];
    int minC[4] = {0,4,8,12};
    int maxC[4] = {4,8,12,20};
    for( int c=0; c<4; c++)
    {
      c6->cd(c + 1);
      hDumComp->Draw();
      gRatCen[c] = get_vnFakRat( f, minC[c], maxC[c], 0.0, 2.4, harmonic, vnCen[c] );
      fCen[c] = new TF1(Form("fCen%d",c),"pol0",0.2,1.0);
      fLin[c] = new TF1(Form("fLin%d",c),"pol1",0.0,1.0);
      gRatCen[c]->Fit(Form("fCen%d",c),"0R");
      gRatCen[c]->Fit(Form("fLin%d",c),"0R");
      double mid = fCen[c]->GetParameter(0);
      double up = mid + fCen[c]->GetParError(0);
      double down = mid - fCen[c]->GetParError(0);
      lmeanCen[c] = new TLine(0.1,mid,1.,mid);
      lmeanCen[c]->SetLineColor(kBlue);
      lmeanCen[c]->SetLineWidth(2);
      lupCen[c] = new TLine(0.1,up,1.,up);
      lupCen[c]->SetLineColor(kGreen);
      lupCen[c]->SetLineWidth(2);
      ldownCen[c] = new TLine(0.1,down,1.,down);
      ldownCen[c]->SetLineColor(kGreen);
      ldownCen[c]->SetLineWidth(2);
      lupCen[c]->Draw();
      lmeanCen[c]->Draw();
      ldownCen[c]->Draw();

      fLin[c]->SetLineWidth(2);
      fLin[c]->SetLineColor(kRed);
      if( doLinear ) fLin[c]->Draw("same");

      gRatCen[c]->SetMarkerStyle(20);
      gRatCen[c]->SetMarkerSize(1.2);
      gRatCen[c]->Draw("p");
      tex = new TLatex(0.15,0.4,Form("%3.1f - %3.1f Cent.",minC[c]*2.5, maxC[c]*2.5 ));
      tex->SetTextSize(0.045);
      tex->Draw();
      
      if( c==0)
      {
        tex = new TLatex(0.15,0.2,"0.0 < |#eta| < 2.4");
        tex->SetTextSize(0.045);
        tex->Draw();
      }
      tex = new TLatex(0.6,2.2,Form("<v_{n}^{EP,reco}> = %6.5f",vnCen[c]));
      tex->SetTextSize(0.045);
      tex->Draw();
      tex = new TLatex(0.6,2.0,Form("SF = %4.3f #pm %4.3f",mid, up-mid));
      tex->SetTextSize(0.045);
      tex->Draw();
      tex = new TLatex(0.55,0.4,Form("pol0 #Chi^{2}/ndof = %5.3f",
                       fCen[c]->GetChisquare()/fCen[c]->GetNDF()));
      tex->SetTextSize(0.045);
      tex->Draw();
      tex = new TLatex(0.55,0.2,Form("pol1 #Chi^{2}/ndof = %5.3f",
                       fLin[c]->GetChisquare()/fLin[c]->GetNDF()));
      tex->SetTextSize(0.045);
      if ( doLinear) tex->Draw();

    }
 
}


TH1D * getSpectrum(TFile * f, char * type, int minCent, int maxCent, double etaMin, double etaMax )
{
  TH3F *h = (TH3F*) f->Get(Form("hitrkPixelEffAnalyzer/h%s3D",type));
  h->GetZaxis()->SetRangeUser( minCent, maxCent);
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

TGraphErrors * get_vno(TFile * f, char * type, int minCent, int maxCent, double etaMin, double etaMax, int h )
{
  TH3F *h_vn[40];
  for( int i=minCent;i<maxCent;i++)
  {
      if ( h < 10 )
        h_vn[i] = (TH3F *) f->Get(Form("hitrkPixelEffAnalyzer/ch%s_v%dRP_%d",type,h,i));
      if ( h == 20 )
        h_vn[i] = (TH3F *) f->Get(Form("hitrkPixelEffAnalyzer/ch%s_v2_%d",type,i));
      if ( h == 42 )
        h_vn[i] = (TH3F *) f->Get(Form("hitrkPixelEffAnalyzer/ch%s_v42e_%d",type,i));
      if ( h == 62 )
        h_vn[i] = (TH3F *) f->Get(Form("hitrkPixelEffAnalyzer/ch%s_v62e_%d",type,i));
      if ( h == 420 )
        h_vn[i] = (TH3F *) f->Get(Form("hitrkPixelEffAnalyzer/ch%s_v42_%d",type,i));
      if ( h == 620 )
        h_vn[i] = (TH3F *) f->Get(Form("hitrkPixelEffAnalyzer/ch%s_v62_%d",type,i));
  }

  TH3F * ch_vn = (TH3F*) h_vn[minCent]->Clone("cch_vn");
  
  for( int i=minCent+1;i<maxCent;i++)
  {
    ch_vn->Add(h_vn[i],1.0);
    delete h_vn[i];
  }

  Int_t etabin10m=ch_vn->GetXaxis()->FindBin(etaMin+0.01);
  Int_t etabin10p=ch_vn->GetXaxis()->FindBin(etaMax-0.01);
  Int_t Netabin10p=ch_vn->GetXaxis()->FindBin(-etaMin-0.01);
  Int_t Netabin10m=ch_vn->GetXaxis()->FindBin(-etaMax+0.01);

  TH2D* hPt_vn;
  TH2D* hPtN_vn;

  ch_vn->GetXaxis()->SetRangeUser(etaMin+0.01,etaMax-0.01);
  hPt_vn = (TH2D*) ch_vn->Project3D("yze");
  ch_vn->GetXaxis()->SetRangeUser(-etaMax+0.01,-etaMin-0.01);
  hPtN_vn = (TH2D*) ch_vn->Project3D("yze");
  hPt_vn->Add( hPtN_vn, 1.0);

  Double_t x[100],xe[100],y[100],ye[100];

  int numBin = hPt_vn->GetNbinsY();
  for( int pTbin = 0; pTbin < numBin; pTbin++)
  {
    TH1D * vn = hPt_vn->ProjectionX("vs",pTbin+1,pTbin+1,"e");
    TH1D * vnx = hPt_vn->ProjectionY("vsx",1,1,"e");
    x[pTbin] = vnx->GetBinCenter(pTbin+1); xe[pTbin] = 0;
    y[pTbin] = vn->GetMean(); ye[pTbin] = vn->GetRMS() / sqrt(vn->GetEntries());
    delete vn;
    delete vnx;
  }

  delete ch_vn;
  delete hPtN_vn; 
//  delete hPt_vn;

  return new TGraphErrors( numBin, x, y, xe, ye ); 
}


TGraphErrors * get_vn(TFile * f, char * type, int minCent, int maxCent, 
                      double etaMin, double etaMax, int h, double er, double ers )
{
  TH3F *h_vn[40];

  if ( h == 20 || h == 420 || h == 620 )
    return get_vno( f, type, minCent, maxCent, etaMin, etaMax, h ); 

  for( int i=minCent;i<maxCent;i++)
  {
      if ( h < 10 )
        h_vn[i] = (TH3F *) f->Get(Form("hitrkPixelEffAnalyzer/ch%s_v%dRP_%d",type,h,i));
      if ( h == 42 )
        h_vn[i] = (TH3F *) f->Get(Form("hitrkPixelEffAnalyzer/ch%s_v42e_%d",type,i));
      if ( h == 62 )
        h_vn[i] = (TH3F *) f->Get(Form("hitrkPixelEffAnalyzer/ch%s_v62e_%d",type,i));
  }

  TH3F * ch_vn = (TH3F*) h_vn[minCent]->Clone("cch_vn");
  
  for( int i=minCent+1;i<maxCent;i++)
  {
    ch_vn->Add(h_vn[i],1.0);
    delete h_vn[i];
  }

  Int_t etabin10m=ch_vn->GetXaxis()->FindBin(etaMin+0.01);
  Int_t etabin10p=ch_vn->GetXaxis()->FindBin(etaMax-0.01);
  Int_t Netabin10p=ch_vn->GetXaxis()->FindBin(-etaMin-0.01);
  Int_t Netabin10m=ch_vn->GetXaxis()->FindBin(-etaMax+0.01);

  TH2D* hPt_vn;
  TH2D* hPtN_vn;

  ch_vn->GetXaxis()->SetRangeUser(etaMin+0.01,etaMax-0.01);
  hPt_vn = (TH2D*) ch_vn->Project3D("yze");
  ch_vn->GetXaxis()->SetRangeUser(-etaMax+0.01,-etaMin-0.01);
  hPtN_vn = (TH2D*) ch_vn->Project3D("yze");
  hPt_vn->Add( hPtN_vn, 1.0);

  Double_t x[100],xe[100],y[100],ye[100];

  int numBin = hPt_vn->GetNbinsY();
  for( int pTbin = 0; pTbin < numBin; pTbin++)
  {
    TH1D * vn = hPt_vn->ProjectionX("vs",pTbin+1,pTbin+1,"e");
    TH1D * vnx = hPt_vn->ProjectionY("vsx",1,1,"e");

    x[pTbin] = vnx->GetBinCenter(pTbin+1); xe[pTbin] = 0;
    double yo = vn->GetMean(); double yoe = vn->GetRMS()/sqrt(vn->GetEntries());
    y[pTbin] = yo / er ;
    ye[pTbin] = fabs(y[pTbin]) * sqrt( yoe*yoe/yo/yo + ers*ers/er/er );
    delete vn;
    delete vnx;
  }

  delete ch_vn;
  delete hPtN_vn; 
//  delete hPt_vn;

  return new TGraphErrors( numBin, x, y, xe, ye ); 
}


TGraphErrors * get_vnFakRat(TFile * f, int minCent, int maxCent, 
                            double etaMin, double etaMax, int h, double & vnInt )
{

  TH1D * hRecPt = getSpectrum( f, "rec", minCent, maxCent, etaMin, etaMax); 
  TH1D * hEffPt = getSpectrum( f, "eff", minCent, maxCent, etaMin, etaMax); 
  TH1D * hFakPt = getSpectrum( f, "fak", minCent, maxCent, etaMin, etaMax); 
  TH1D * hSimPt = getSpectrum( f, "sim", minCent, maxCent, etaMin, etaMax); 
  TGraphErrors * gRec_vno = get_vno( f, "rec", minCent, maxCent, etaMin, etaMax, h );
  TGraphErrors * gFak_vno = get_vno( f, "fak", minCent, maxCent, etaMin, etaMax, h );

  // Integration range is hard coded here as 0.3<pT<3.0
  int lim_ptL = hRecPt->GetXaxis()->FindBin(0.305);
  int lim_ptH = hRecPt->GetXaxis()->FindBin(2.995);
  double corrected[100];
  double corrint = 0.0;
  vnInt = 0.0;
  for( int i=lim_ptL; i<=lim_ptH; i++)
  {
    double xDum, yDum;
    corrected[i] = hRecPt->GetBinContent(i);
    double eff = hEffPt->GetBinContent(i)/hSimPt->GetBinContent(i);
    double fak = hFakPt->GetBinContent(i)/hRecPt->GetBinContent(i);
    corrected[i] = corrected[i] * (1 - fak ) / eff;
    corrint += corrected[i];
    gRec_vno->GetPoint(i-1,xDum,yDum);
    vnInt += yDum * corrected[i];
  }
  vnInt /= corrint;

  Double_t xc3[100], xec3[100], yc3[100], yec3[100];

  for( int p=0; p<gFak_vno->GetN(); p++)
  {
    Double_t ay, aye,;
    gFak_vno->GetPoint(p,xc3[p],ay);
    aye = gFak_vno->GetErrorY(p);
    xec3[p] = 0.0;
    yc3[p] = ay / vnInt;
    yec3[p] = fabs(yc3[p]) * sqrt( aye*aye/ay/ay );
  }
 
  int n = gFak_vno->GetN();

  delete hFakPt;
  delete hEffPt;
  delete hRecPt;
  delete hSimPt;
  delete gRec_vno;
  delete gFak_vno;

  return new TGraphErrors( n, xc3, yc3, xec3, yec3 );
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

