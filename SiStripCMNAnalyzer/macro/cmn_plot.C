{

  // ------- CONFIGURATION -----------------

  // Name of gallery to plot
  //
  // galA - calculated baseline 116-140 for 25th percentile and median
  // galB - calculated baseline 2-100 for 25th percentile and median
  // galD - calculated baseline >160 for 25th percentile and median
  // galE - calculated baseline >160 for median but 116-140 for 25th percentile
 
  char gallery[10] = "galE";

  // File to analyze

  TFile *f1 = TFile::Open("cmn_1sigB.root");

  // Plot signal-only rawDigis (for hacked Digitizer events only)

  bool gensignal = false;

  // ----------------------------------------


  gStyle->SetCanvasColor(kWhite);     // background is no longer mouse-dropping white

  gStyle->SetTitleOffset(1.0,"x");
  gStyle->SetTitleOffset(1.25,"y");

  char galTitle[255] = "APV25 plot #";

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);

  TCanvas * c[4];
  TH1F * qmedians[4];
  TH1F * imedians[4];
  TH1F * quarters[4];
  TH1F * gens[4];
  TLegend * leg[4];

  char detName[100];  
  
  for( int i = 0; i<4; i++)
  {
    switch (i)
    { 
      case 0:
        sprintf(detName,"TEC");
        break;
      case 1:
        sprintf(detName,"TIB");
        break;
      case 2:
        sprintf(detName,"TOB");
        break;
      case 3:
        sprintf(detName,"TID");
        break;
    } 
    c[i] = new TCanvas(Form("c%s",detName),Form("c%s",detName),800.,600.);
    c[i]->SetFillColor(kWhite);
    c[i]->SetFrameFillColor(kWhite);
    c[i]->SetBorderMode(0);
    c[i]->SetFrameBorderMode(0);
    c[i]->cd();
    qmedians[i] = (TH1F*) f1->Get(Form("cmn/%sMedians",detName));
    qmedians[i]->Rebin(5);
    imedians[i] = (TH1F*) f1->Get(Form("cmn/%sIterMeds",detName));
    imedians[i]->Rebin(5);
    quarters[i] = (TH1F*) f1->Get(Form("cmn/%s25ths",detName));
    quarters[i]->Rebin(5);
    gens[i] = (TH1F*) f1->Get(Form("cmn/%sGens",detName));
    gens[i]->Rebin(5);
    gens[i]->SetLineColor(kGreen);
    gens[i]->SetLineWidth(2);
    gens[i]->SetTitle(Form("ADC Common Mode Noise Offset for %s",detName));
    gens[i]->GetXaxis()->SetTitle("Offset");
    gens[i]->SetFillColor(kGreen);
    gens[i]->SetFillStyle(3020);
    gens[i]->Draw();
    quarters[i]->SetLineColor(kRed);
    quarters[i]->SetLineWidth(2);
    imedians[i]->SetLineColor(kViolet);
    imedians[i]->SetLineWidth(2);
    qmedians[i]->SetLineWidth(2);
    quarters[i]->Draw("same");
    qmedians[i]->Draw("same");
    imedians[i]->Draw("same");
    c[i]->SetLogy();
    leg[i] = new TLegend(0.68,0.65,0.98,0.85);
    leg[i]->AddEntry(qmedians[i],"Median Algorithm");
    leg[i]->AddEntry(quarters[i],"25th Percentile Algorithm");
    leg[i]->AddEntry(imedians[i],"Iterated Median Algorithm");
    leg[i]->AddEntry(gens[i],"Generated CMN");
    leg[i]->SetFillColor(0);
    leg[i]->Draw();
  }

  TCanvas * cComp = new TCanvas("cComp","cComp",800.,600.);
  cComp->SetFillColor(kWhite);
  cComp->SetFrameFillColor(kWhite);
  cComp->SetBorderMode(0);
  cComp->SetFrameBorderMode(0);
  cComp->cd();

  

  TH2F * comparison = (TH2F*) f1->Get("cmn/medvs25");
  comparison->RebinX(2);
  comparison->RebinY(5);
  comparison->SetTitle("Comparison of algorithms for all subdetectors");
  comparison->GetXaxis()->SetTitle("Median Algorithm Offset");
  comparison->GetYaxis()->SetTitle("Percentile (25) Algorithm Offset");
  comparison->Draw("colz");
  cComp->SetLogz();

  TCanvas * cWidth = new TCanvas("cWidth","cWidth",800.,600.);
  cWidth->SetFillColor(kWhite);
  cWidth->SetFrameFillColor(kWhite);
  cWidth->SetBorderMode(0);
  cWidth->SetFrameBorderMode(0);
  cWidth->cd();
  cWidth->SetLogy();


  TH1F * medCW = (TH1F*) f1->Get("cmn/medianClusWidth"); 
  medCW->SetLineColor(kBlue);
  medCW->GetXaxis()->SetTitle("Width in Strips");
  medCW->GetXaxis()->SetLimits(0.,150.);
  medCW->SetTitle("Cluster Width Comparison for Three CMN Algorithms");
  //medCW->SetTitleFillColor(0);
  medCW->SetLineWidth(2);
  medCW->Draw();
  TH1F * iterMedCW = (TH1F*) f1->Get("cmn/iterMedClusWidth"); 
  iterMedCW->SetLineColor(kViolet);
  iterMedCW->SetLineWidth(2);
  iterMedCW->Draw("same");
  TH1F * per25CW = (TH1F*) f1->Get("cmn/per25ClusWidth"); 
  per25CW->SetLineColor(kRed);
  per25CW->SetLineWidth(2);
  per25CW->Draw("same");
  legCW = new TLegend(0.68,0.65,0.98,0.85);
  legCW->AddEntry(medCW,"Median Algorithm");
  legCW->AddEntry(per25CW,"25th Percentile Algorithm");
  legCW->AddEntry(iterMedCW,"Iterated Median Algorithm");
  legCW->SetFillColor(0);
  legCW->Draw();

 
  TCanvas * bad[10];
  TGraph * bg[10];
  TH1F * bgf[10];
  TGraph * bsig[10];
  TLine * median[10];
  TLine * percent25[10];
  TLine * imedian[10];
  
  TGraph * clusM[10];
  TGraph * clusIM[10];
  TGraph * clus25[10];
  TLine *  clusterLineM[1280];
  TLine *  clusterLineIM[1280];
  TLine *  clusterLine25[1280];
  TLine *  clusterLineMF[1280];
  TLine *  clusterLineIMF[1280];
  TLine *  clusterLine25F[1280];
  TLine *  clusterLineME[1280];
  TLine *  clusterLineIME[1280];
  TLine *  clusterLine25E[1280];

  for( int j=0; j<10; j++)
  {
    bad[j] = new TCanvas(Form("bad%d",j),Form("bad%d",j),1200.,500.);
    bad[j]->Divide(2,1,0.005,0.005);
    bad[j]->SetFillColor(kWhite);
    bad[j]->SetFrameFillColor(kWhite);
    bad[j]->SetBorderMode(0);
    bad[j]->SetFrameBorderMode(0);
    bad[j]->cd(1);
    bg[j] = (TGraph*) f1->Get(Form("cmn/%sCount%d",gallery,j));
    if ( gensignal)
      bsig[j] = (TGraph*) f1->Get(Form("cmn/%sSignalCount%d",gallery,j));
    bg[j]->SetTitle(Form("%s%d",galTitle,j));
    Double_t dummy;
    Double_t med;
    Double_t per25;
    Double_t imed;
    bg[j]->GetPoint(129,dummy,per25);
    bg[j]->GetPoint(128,dummy,med);
    bg[j]->GetPoint(130,dummy,imed);
    for( int i=139;i>127;i--) bg[j]->RemovePoint(i);
    bg[j]->GetXaxis()->SetLimits(0.,128.);
    bg[j]->GetYaxis()->SetRangeUser(0.,900.);
    bg[j]->GetYaxis()->SetTitle("ADC Count");
    bg[j]->GetXaxis()->SetTitle("Si Strip Number");
    median[j] = new TLine(0,med,128,med);
    imedian[j] = new TLine(0,imed,128,imed);
    percent25[j] = new TLine(0,per25,128,per25);
    median[j]->SetLineColor(kBlue);
    imedian[j]->SetLineColor(kViolet);
    percent25[j]->SetLineColor(kRed);
    bg[j]->SetMarkerStyle(7);
    if (gensignal) 
    {
      bsig[j]->SetMarkerStyle(7);
      bsig[j]->SetLineColor(kGreen);
      bsig[j]->SetMarkerColor(kGreen);
    }
    bg[j]->Draw("apl");
    if (gensignal)
      bsig[j]->Draw("pl");
    median[j]->Draw("same");
    imedian[j]->Draw("same");
    percent25[j]->Draw("same");
    clusM[j] = (TGraph*) f1->Get(Form("cmn/%sMedianClusters%d",gallery,j));
    clusIM[j] = (TGraph*) f1->Get(Form("cmn/%sIterMedClusters%d",gallery,j));
    clus25[j] = (TGraph*) f1->Get(Form("cmn/%sPer25Clusters%d",gallery,j));
    for( int k = 0; k < 128; k++ )
    {
      Double_t first;
      Double_t last;
      clusM[j]->GetPoint(k,first,last);
      if ( first != 0. && last != 0. )
      {
        clusterLineM[j*128+k] = new TLine(first, 800, last, 800);
        clusterLineM[j*128+k]->SetLineColor(kBlue);
        clusterLineM[j*128+k]->SetLineWidth(2);
        clusterLineM[j*128+k]->Draw("same");
        clusterLineMF[j*128+k] = new TLine(first, 775, first, 825);
        clusterLineMF[j*128+k]->SetLineColor(kBlue);
        clusterLineMF[j*128+k]->SetLineWidth(2);
        clusterLineMF[j*128+k]->Draw("same");
        clusterLineME[j*128+k] = new TLine(last, 775, last, 825);
        clusterLineME[j*128+k]->SetLineColor(kBlue);
        clusterLineME[j*128+k]->SetLineWidth(2);
        clusterLineME[j*128+k]->Draw("same");
      }
    }
    for( int k = 0; k < 128; k++ )
    {
      Double_t first;
      Double_t last;
      clusIM[j]->GetPoint(k,first,last);
      if ( first != 0. && last != 0. )
      {
        clusterLineIM[j*128+k] = new TLine(first, 750, last, 750);
        clusterLineIM[j*128+k]->SetLineColor(kViolet);
        clusterLineIM[j*128+k]->SetLineWidth(2);
        clusterLineIM[j*128+k]->Draw("same");
        clusterLineIMF[j*128+k] = new TLine(first, 725, first, 775);
        clusterLineIMF[j*128+k]->SetLineColor(kViolet);
        clusterLineIMF[j*128+k]->SetLineWidth(2);
        clusterLineIMF[j*128+k]->Draw("same");
        clusterLineIME[j*128+k] = new TLine(last, 725, last, 775);
        clusterLineIME[j*128+k]->SetLineColor(kViolet);
        clusterLineIME[j*128+k]->SetLineWidth(2);
        clusterLineIME[j*128+k]->Draw("same");
      }
    }
    for( int k = 0; k < 128; k++ )
    {
      Double_t first;
      Double_t last;
      clus25[j]->GetPoint(k,first,last);
      if ( first != 0. && last != 0. )
      {
        clusterLine25[j*128+k] = new TLine(first, 700, last, 700);
        clusterLine25[j*128+k]->SetLineColor(kRed);
        clusterLine25[j*128+k]->SetLineWidth(2);
        clusterLine25[j*128+k]->Draw("same");
        clusterLine25F[j*128+k] = new TLine(first, 675, first, 725);
        clusterLine25F[j*128+k]->SetLineColor(kRed);
        clusterLine25F[j*128+k]->SetLineWidth(2);
        clusterLine25F[j*128+k]->Draw("same");
        clusterLine25E[j*128+k] = new TLine(last, 675, last, 725);
        clusterLine25E[j*128+k]->SetLineColor(kRed);
        clusterLine25E[j*128+k]->SetLineWidth(2);
        clusterLine25E[j*128+k]->Draw("same");
      }
    }

    bad[j]->cd(2);
    bgf[j] = new TH1F(Form("bgf%d",j),"",100,0.,500.);
    for( int k = 0; k<128; k++ )
    {
      Double_t x;
      Double_t y;
      bg[j]->GetPoint(k,x,y);
      bgf[j]->Fill( y );
    }
    bgf[j]->GetXaxis()->SetTitle("ADC Count");
    bgf[j]->GetYaxis()->SetTitle("frequency");
    bgf[j]->Draw();
   // bad[j]->GetPad(2)->SetLogy();
    Double_t rx1, rx2, ry1, ry2;
    bad[j]->GetRange( rx1, ry1, rx2, ry2 );
    bad[j]->GetPad(1)->SetPad(rx1, ry1, rx2 * 0.6, ry2);
    bad[j]->GetPad(2)->SetPad(rx2 * 0.6, ry1, rx2, ry2);
    bad[j]->GetPad(1)->SetFrameBorderMode(0); 
    bad[j]->GetPad(2)->SetFrameBorderMode(0); 
  }
}
