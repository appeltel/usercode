// Macro for plotting and comparing the output of several
// RpPbTrackAnalyzer outputs.
//
// Eric Appelt, Jan 2014, based on Rylan Conway's eariler macro 
// from Jan 2013.

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include <cstdarg>

using namespace std;

void drawTrackPlots( double ptMin, double ptMax, 
                     double etaMin, double etaMax,  
                     string refFile, string refDir, string refName,
                     string anaFile, string anaDir, string anaName,
                     int numExtraFiles, ... )
{

  /////////////
  // setup plotting style

  double msize = 0.85;
  int mstyle[10] = {20,20,20,20,20,20,20,20,20,20};
  int mcol[10] = {kBlack, kRed, kBlue, kGreen, kViolet, kOrange, kGray, kCyan, kYellow, kPink};
  gStyle->SetOptStat(0);

  /////////////
  // open data files

  TFile * f[10];
  char * dir[10];
  char * dataName[10];
 
  f[0] = new TFile(refFile.c_str());
  f[1] = new TFile(anaFile.c_str());
  dir[0] = refDir.c_str();
  dir[1] = anaDir.c_str();
  dataName[0] = refName.c_str(); 
  dataName[1] = anaName.c_str(); 

  if( numExtraFiles > 8 ) 
    cout << "Warning: 10 datasets maximum! Additional files will be ignored." << endl;

  va_list extraFiles;
  va_start ( extraFiles, numExtraFiles*3 );
  for( int x = 0; x < numExtraFiles*3; x += 3)
  {
    if( x >= 24) break;
    string fstr = va_arg( extraFiles, string );
    f[x/3+2] = new TFile(fstr.c_str()); 
    string dstr = va_arg(extraFiles, string );    
    dir[x/3+2] = dstr.c_str();
    string astr = va_arg(extraFiles, string );
    dataName[x/3+2] = astr.c_str();
  }
  va_end ( extraFiles );  
  
  int nDatasets = numExtraFiles+2;  

  //////////////
  // load track performance histograms
  // and related data

  cout << "nDatasets = " << nDatasets << endl; 
             
  TH1D * nEvt[10], * nTrk[10], * nVtx[10];

  TH1D * vtxX[10];
  TH1D * vtxY[10];
  TH1D * vtxZ[10];
  TH2D * vtxNtrk2D[10];
  TH2D * trkEtaPhi[10];
  
  TH3F * spectrum[10];
  TH3F * trkNhit3D[10];
  TH3F * trkPtErr3D[10];
  TH3F * trkPhi3D[10];
  TH3F * trkChi23D[10];
  TH3F * trkDxyErr3D[10];
  TH3F * trkDzErr3D[10];
   
  TH1D * trkPt[10];
  TH1D * trkEta[10];
  TH1D * trkNhit[10];
  TH1D * trkPtErr[10];
  TH1D * trkPhi[10];
  TH1D * trkChi2[10];
  TH1D * trkDxyErr[10];
  TH1D * trkDzErr[10];

  TH1D * ratPt[10];
  TH1D * ratEta[10];
  TH1D * ratNhit[10];
  TH1D * ratPtErr[10];
  TH1D * ratPhi[10];
  TH1D * ratChi2[10];
  TH1D * ratDxyErr[10];
  TH1D * ratDzErr[10];

  for( int i=0; i<nDatasets; i++ )
  {
    cout << "Grabbing histos i= " << i << endl;
    nEvt[i] = (TH1D*) f[i]->Get(Form("%s/events",dir[i]));
    nTrk[i]  = (TH1D*) f[i]->Get(Form("%s/evtNtrk",dir[i]));
    nVtx[i]  = (TH1D*) f[i]->Get(Form("%s/evtNvtx",dir[i])); 

    vtxX[i] = (TH1D*) f[i]->Get(Form("%s/vtxX",dir[i]));
    vtxY[i] = (TH1D*) f[i]->Get(Form("%s/vtxY",dir[i]));
    vtxZ[i] = (TH1D*) f[i]->Get(Form("%s/vtxZ",dir[i]));
    vtxNtrk2D[i] = (TH2D*) f[i]->Get(Form("%s/vtxNtrk2D",dir[i]));

    spectrum[i] = (TH3F*) f[i]->Get(Form("%s/trkSpectrum",dir[i]));
    trkNhit3D[i] = (TH3F*) f[i]->Get(Form("%s/trkNhit3D",dir[i]));
    trkPtErr3D[i] = (TH3F*) f[i]->Get(Form("%s/trkPterr3D",dir[i]));
    trkPhi3D[i] = (TH3F*) f[i]->Get(Form("%s/trkPhi3D",dir[i]));
    trkChi23D[i] = (TH3F*) f[i]->Get(Form("%s/trkChi23D",dir[i]));
    trkDxyErr3D[i] = (TH3F*) f[i]->Get(Form("%s/trkDxyErr3D",dir[i]));
    trkDzErr3D[i] = (TH3F*) f[i]->Get(Form("%s/trkDzErr3D",dir[i]));

    nEvt[i]->Sumw2();
    nTrk[i]->Sumw2();
    nVtx[i]->Sumw2();
    vtxX[i]->Sumw2();
    vtxY[i]->Sumw2();
    vtxZ[i]->Sumw2();
    vtxNtrk2D[i]->Sumw2();
    spectrum[i]->Sumw2();
    trkNhit3D[i]->Sumw2();
    trkPtErr3D[i]->Sumw2();
    trkPhi3D[i]->Sumw2();
    trkChi23D[i]->Sumw2();
    trkDxyErr3D[i]->Sumw2();
    trkDzErr3D[i]->Sumw2();
     
    // run projections

    trkPt[i] = spectrum[i]->ProjectionY(Form("trkPt_%d",i),
               spectrum[i]->GetXaxis()->FindBin(etaMin+0.001),
               spectrum[i]->GetXaxis()->FindBin(etaMax-0.001),
               0,-1,"e");
    trkEta[i] = spectrum[i]->ProjectionX(Form("trkEta_%d",i),
               spectrum[i]->GetYaxis()->FindBin(ptMin+0.001),
               spectrum[i]->GetYaxis()->FindBin(ptMax-0.001),
               0,-1,"e");

    trkNhit[i] = trkNhit3D[i]->ProjectionZ(Form("trkNhit_%d",i),
                 trkNhit3D[i]->GetXaxis()->FindBin(etaMin+0.001),
                 trkNhit3D[i]->GetXaxis()->FindBin(etaMax-0.001),
                 trkNhit3D[i]->GetYaxis()->FindBin(ptMin+0.001),
                 trkNhit3D[i]->GetYaxis()->FindBin(ptMax-0.001), "e");
    trkPtErr[i] = trkPtErr3D[i]->ProjectionZ(Form("trkPtErr_%d",i),
                 trkPtErr3D[i]->GetXaxis()->FindBin(etaMin+0.001),
                 trkPtErr3D[i]->GetXaxis()->FindBin(etaMax-0.001),
                 trkPtErr3D[i]->GetYaxis()->FindBin(ptMin+0.001),
                 trkPtErr3D[i]->GetYaxis()->FindBin(ptMax-0.001), "e");
    trkPhi[i] = trkPhi3D[i]->ProjectionZ(Form("trkPhi_%d",i),
                 trkPhi3D[i]->GetXaxis()->FindBin(etaMin+0.001),
                 trkPhi3D[i]->GetXaxis()->FindBin(etaMax-0.001),
                 trkPhi3D[i]->GetYaxis()->FindBin(ptMin+0.001),
                 trkPhi3D[i]->GetYaxis()->FindBin(ptMax-0.001), "e");
    trkChi2[i] = trkChi23D[i]->ProjectionZ(Form("trkChi2_%d",i),
                 trkChi23D[i]->GetXaxis()->FindBin(etaMin+0.001),
                 trkChi23D[i]->GetXaxis()->FindBin(etaMax-0.001),
                 trkChi23D[i]->GetYaxis()->FindBin(ptMin+0.001),
                 trkChi23D[i]->GetYaxis()->FindBin(ptMax-0.001), "e");
    trkDxyErr[i] = trkDxyErr3D[i]->ProjectionZ(Form("trkDxyErr_%d",i),
                 trkDxyErr3D[i]->GetXaxis()->FindBin(etaMin+0.001),
                 trkDxyErr3D[i]->GetXaxis()->FindBin(etaMax-0.001),
                 trkDxyErr3D[i]->GetYaxis()->FindBin(ptMin+0.001),
                 trkDxyErr3D[i]->GetYaxis()->FindBin(ptMax-0.001), "e");
    trkDzErr[i] = trkDzErr3D[i]->ProjectionZ(Form("trkDzErr_%d",i),
                 trkDzErr3D[i]->GetXaxis()->FindBin(etaMin+0.001),
                 trkDzErr3D[i]->GetXaxis()->FindBin(etaMax-0.001),
                 trkDzErr3D[i]->GetYaxis()->FindBin(ptMin+0.001),
                 trkDzErr3D[i]->GetYaxis()->FindBin(ptMax-0.001), "e");

    trkPhi3D[i]->GetYaxis()->SetRangeUser(etaMin+0.001, etaMin-0.001);
    trkEtaPhi[i] =  (TH2D*) trkPhi3D[i]->Project3D("xze");
    trkEtaPhi[i]->SetName(Form("trkEtaPhi_%d",i));
  
    // normalize histograms

    nTrk[i]->Scale(1./nTrk[i]->Integral());
    nVtx[i]->Scale(1./nVtx[i]->Integral());
    vtxX[i]->Scale(1./vtxX[i]->Integral());
    vtxY[i]->Scale(1./vtxY[i]->Integral());
    vtxZ[i]->Scale(1./vtxZ[i]->Integral());
    trkPt[i]->Scale(1./trkPt[i]->Integral());
    trkEta[i]->Scale(1./trkEta[i]->Integral());
    trkNhit[i]->Scale(1./trkNhit[i]->Integral());
    trkPhi[i]->Scale(1./trkPhi[i]->Integral());
    trkPtErr[i]->Scale(1./trkPtErr[i]->Integral());
    trkChi2[i]->Scale(1./trkChi2[i]->Integral());
    trkDxyErr[i]->Scale(1./trkDxyErr[i]->Integral());
    trkDzErr[i]->Scale(1./trkDzErr[i]->Integral());

  }

  //////////////////
  // create ratio hists

  for( int i=1; i<nDatasets; i++)
  { 
    ratPt[i] = (TH1D*) trkPt[i]->Clone(Form("ratPt_%d",i));
    ratPt[i]->Divide(trkPt[0]);
    ratEta[i] = (TH1D*) trkEta[i]->Clone(Form("ratEta_%d",i));
    ratEta[i]->Divide(trkEta[0]);
    ratNhit[i] = (TH1D*) trkNhit[i]->Clone(Form("ratNhit_%d",i));
    ratNhit[i]->Divide(trkNhit[0]);
    ratPhi[i] = (TH1D*) trkPhi[i]->Clone(Form("ratPhi_%d",i));
    ratPhi[i]->Divide(trkPhi[0]);
    ratPtErr[i] = (TH1D*) trkPtErr[i]->Clone(Form("ratPtErr_%d",i));
    ratPtErr[i]->Divide(trkPtErr[0]);
    ratChi2[i] = (TH1D*) trkChi2[i]->Clone(Form("ratChi2_%d",i));
    ratChi2[i]->Divide(trkChi2[0]);
    ratDxyErr[i] = (TH1D*) trkDxyErr[i]->Clone(Form("ratDxyErr_%d",i));
    ratDxyErr[i]->Divide(trkDxyErr[0]);
    ratDzErr[i] = (TH1D*) trkDzErr[i]->Clone(Form("ratDzErr_%d",i));
    ratDzErr[i]->Divide(trkDzErr[0]);
  }

  ///////////////////
  // Vertex properties plot

  TCanvas* c0 = new TCanvas("c0","Vertex Information", 850,750);
  c0->Divide(2,2);

  c0->cd(1); 
  c0->GetPad(1)->SetLogy();
  nVtx[0]->SetAxisRange(1.0e-5,1.3,"Y");
  for( int i=0;i<nDatasets;i++)
  {
    nVtx[i]->SetXTitle("Number of vertices per event");
    nVtx[i]->SetMarkerSize(msize);
    nVtx[i]->SetLineWidth(msize);
    nVtx[i]->SetMarkerStyle(mstyle[i]);
    nVtx[i]->SetMarkerColor(mcol[i]);
    if( i==0) nVtx[i]->Draw("PE");
    else nVtx[i]->Draw("PSAME");
  }
  TLegend* leg = new TLegend(0.5, 0.7, 0.8, 0.85);
  for( int i=0;i<nDatasets;i++)
    leg->AddEntry(nVtx[i],dataName[i],"p");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw("SAME");

  c0->cd(2);
  for( int i=0;i<nDatasets;i++)
  {
    vtxX[i]->SetXTitle("V_{x}");
    vtxX[i]->SetMarkerSize(msize);
    vtxX[i]->SetLineWidth(msize);
    vtxX[i]->SetMarkerStyle(mstyle[i]);
    vtxX[i]->SetMarkerColor(mcol[i]);
    if( i==0) vtxX[i]->Draw("HIST");
    else vtxX[i]->Draw("HISTSAME");
  }
  gPad->Update();

  c0->cd(3);
  for( int i=0;i<nDatasets;i++)
  {
    vtxY[i]->SetXTitle("V_{y}");
    vtxY[i]->SetMarkerSize(msize);
    vtxY[i]->SetLineWidth(msize);
    vtxY[i]->SetMarkerStyle(mstyle[i]);
    vtxY[i]->SetMarkerColor(mcol[i]);
    if( i==0) vtxY[i]->Draw("HIST");
    else vtxY[i]->Draw("HISTSAME");
  }
  gPad->Update();

  c0->cd(4);
  for( int i=0;i<nDatasets;i++)
  {
    vtxZ[i]->SetXTitle("V_{z}");
    vtxZ[i]->SetMarkerSize(msize);
    vtxZ[i]->SetLineWidth(msize);
    vtxZ[i]->SetMarkerStyle(mstyle[i]);
    vtxZ[i]->SetMarkerColor(mcol[i]);
    if( i==0) vtxZ[i]->Draw("PE");
    else vtxZ[i]->Draw("PESAME");
  }
  gPad->Update();


  ////////////////
  // Tracking pt, eta, phi distributions

  TCanvas* c1 = new TCanvas("c1","Track Variables 1", 1100,750);
  c1->Divide(3,2);

  c1->cd(1);
  c1->GetPad(1)->SetLogx();
  c1->GetPad(1)->SetLogy();
  trkPt[0]->GetXaxis()->SetRangeUser(0.1,100);
  trkPt[0]->SetTitle("Track p_{T} Distribution");
  for( int i=0;i<nDatasets;i++)
  {
    trkPt[i]->SetMarkerSize(msize);
    trkPt[i]->SetLineWidth(msize);
    trkPt[i]->SetMarkerStyle(mstyle[i]);
    trkPt[i]->SetMarkerColor(mcol[i]);
    if( i==0) trkPt[i]->Draw("PE");
    else trkPt[i]->Draw("PESAME");
  }
  TLegend* leg1 = new TLegend(0.5, 0.7, 0.8, 0.85);
  for( int i=0;i<nDatasets;i++)
    leg1->AddEntry(trkPt[i],dataName[i],"p");
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->Draw("SAME");

  c1->cd(2);
  trkPhi[0]->SetTitle("Track #phi Distribution");
  trkPhi[0]->GetXaxis()->SetTitle("#phi");
  for( int i=0;i<nDatasets;i++)
  {
    trkPhi[i]->SetMarkerSize(msize);
    trkPhi[i]->SetLineWidth(msize);
    trkPhi[i]->SetMarkerStyle(mstyle[i]);
    trkPhi[i]->SetMarkerColor(mcol[i]);
    if( i==0) trkPhi[i]->Draw("PE");
    else trkPhi[i]->Draw("PESAME");
  }

  c1->cd(3);
  trkEta[0]->SetTitle("Track #eta Distribution");
  trkEta[0]->GetXaxis()->SetTitle("#eta");
  for( int i=0;i<nDatasets;i++)
  {
    trkEta[i]->SetMarkerSize(msize);
    trkEta[i]->SetLineWidth(msize);
    trkEta[i]->SetMarkerStyle(mstyle[i]);
    trkEta[i]->SetMarkerColor(mcol[i]);
    if( i==0) trkEta[i]->Draw("PE");
    else trkEta[i]->Draw("PESAME");
  }

  c1->cd(4);
  ratPt[1]->SetAxisRange(0.7,1.3,"Y");
  c1->GetPad(4)->SetLogx();
  ratPt[1]->GetXaxis()->SetRangeUser(0.1,100);
  ratPt[1]->SetTitle("Track p_{T} Ratio");
  for( int i=1;i<nDatasets;i++)
  {
    ratPt[i]->SetMarkerSize(msize);
    ratPt[i]->SetLineWidth(msize);
    ratPt[i]->SetMarkerStyle(mstyle[i]);
    ratPt[i]->SetMarkerColor(mcol[i]);
    if( i==1) ratPt[i]->Draw("HIST P");
    else ratPt[i]->Draw("HIST P SAME");
  }
  TLegend* legR = new TLegend(0.14, 0.7, 0.75, 0.9);
  for( int i=1;i<nDatasets;i++)
  {
    legR->AddEntry(ratPt[i],Form("%s / %s",dataName[i],dataName[0]), "p");
  }
  legR->SetBorderSize(0);
  legR->SetFillColor(0);
  legR->Draw("PESAME");

  c1->cd(5);
  ratPhi[1]->SetAxisRange(0.7,1.3,"Y");
  ratPhi[1]->SetTitle("Track #phi Ratio");
  for( int i=1;i<nDatasets;i++)
  {
    ratPhi[i]->SetMarkerSize(msize);
    ratPhi[i]->SetLineWidth(msize);
    ratPhi[i]->SetMarkerStyle(mstyle[i]);
    ratPhi[i]->SetMarkerColor(mcol[i]);
    if( i==1) ratPhi[i]->Draw("HIST P");
    else ratPhi[i]->Draw("HIST P SAME");
  }

  c1->cd(6);
  ratEta[1]->SetAxisRange(0.7,1.3,"Y");
  ratEta[1]->SetTitle("Track #eta  Ratio");
  for( int i=1;i<nDatasets;i++)
  {
    ratEta[i]->SetMarkerSize(msize);
    ratEta[i]->SetLineWidth(msize);
    ratEta[i]->SetMarkerStyle(mstyle[i]);
    ratEta[i]->SetMarkerColor(mcol[i]);
    if( i==1) ratEta[i]->Draw("HIST P");
    else ratEta[i]->Draw("HIST P SAME");
  }

  ////////////////
  // Tracking pterr, dxyerr, dzerr distributions

  TCanvas* c2 = new TCanvas("c2","Track Variables 2", 1100,750);
  c2->Divide(3,2);

  c2->cd(1);
  trkPtErr[0]->SetTitle("Track p_{T} error Distribution");
  for( int i=0;i<nDatasets;i++)
  {
    trkPtErr[i]->SetMarkerSize(msize);
    trkPtErr[i]->SetLineWidth(msize);
    trkPtErr[i]->SetMarkerStyle(mstyle[i]);
    trkPtErr[i]->SetMarkerColor(mcol[i]);
    if( i==0) trkPtErr[i]->Draw("PE");
    else trkPtErr[i]->Draw("PESAME");
  }
  TLegend* leg2 = new TLegend(0.5, 0.7, 0.8, 0.85);
  for( int i=0;i<nDatasets;i++)
    leg2->AddEntry(trkPtErr[i],dataName[i],"p");
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->Draw("SAME");

  c2->cd(2);
  trkDxyErr[0]->SetTitle("Transverse DCA Signficance");
  for( int i=0;i<nDatasets;i++)
  {
    trkDxyErr[i]->SetMarkerSize(msize);
    trkDxyErr[i]->SetLineWidth(msize);
    trkDxyErr[i]->SetMarkerStyle(mstyle[i]);
    trkDxyErr[i]->SetMarkerColor(mcol[i]);
    if( i==0) trkDxyErr[i]->Draw("PE");
    else trkDxyErr[i]->Draw("PESAME");
  }

  c2->cd(3);
  trkDzErr[0]->SetTitle("Longitudinal DCA Signficance");
  for( int i=0;i<nDatasets;i++)
  {
    trkDzErr[i]->SetMarkerSize(msize);
    trkDzErr[i]->SetLineWidth(msize);
    trkDzErr[i]->SetMarkerStyle(mstyle[i]);
    trkDzErr[i]->SetMarkerColor(mcol[i]);
    if( i==0) trkDzErr[i]->Draw("PE");
    else trkDzErr[i]->Draw("PESAME");
  }

  c2->cd(4);
  ratPtErr[1]->SetTitle("Track p_{T} error Ratio");
  ratPtErr[1]->SetAxisRange(0.7,1.3,"Y");
  for( int i=1;i<nDatasets;i++)
  {
    ratPtErr[i]->SetMarkerSize(msize);
    ratPtErr[i]->SetLineWidth(msize);
    ratPtErr[i]->SetMarkerStyle(mstyle[i]);
    ratPtErr[i]->SetMarkerColor(mcol[i]);
    if( i==0) ratPtErr[i]->Draw("PE");
    else ratPtErr[i]->Draw("PESAME");
  }
  TLegend* legR2 = new TLegend(0.14, 0.7, 0.8, 0.85);
  for( int i=1;i<nDatasets;i++)
    legR2->AddEntry(ratPtErr[i],Form("%s / %s",dataName[i],dataName[0]), "p");
  legR2->SetBorderSize(0);
  legR2->SetFillColor(0);
  legR2->Draw("SAME");

  c2->cd(5);
  ratDxyErr[1]->SetTitle("Transverse DCA Ratio");
  ratDxyErr[1]->SetAxisRange(0.7,1.3,"Y");
  for( int i=1;i<nDatasets;i++)
  {
    ratDxyErr[i]->SetMarkerSize(msize);
    ratDxyErr[i]->SetLineWidth(msize);
    ratDxyErr[i]->SetMarkerStyle(mstyle[i]);
    ratDxyErr[i]->SetMarkerColor(mcol[i]);
    if( i==0) ratDxyErr[i]->Draw("PE");
    else ratDxyErr[i]->Draw("PESAME");
  }

  c2->cd(6);
  ratDzErr[1]->SetTitle("Longitudinal DCA Ratio");
  ratDzErr[1]->SetAxisRange(0.7,1.3,"Y");
  for( int i=1;i<nDatasets;i++)
  {
    ratDzErr[i]->SetMarkerSize(msize);
    ratDzErr[i]->SetLineWidth(msize);
    ratDzErr[i]->SetMarkerStyle(mstyle[i]);
    ratDzErr[i]->SetMarkerColor(mcol[i]);
    if( i==0) ratDzErr[i]->Draw("PE");
    else ratDzErr[i]->Draw("PESAME");
  }

  ////////////////
  // Tracking Nhit, Chi2 distributions
  // and eta-phi maps (only 1st two datasets)

  TCanvas* c3 = new TCanvas("c3","Track Variables 3", 1100,750);
  c3->Divide(3,2);

  c3->cd(1);
  trkChi2[0]->SetTitle("Track Chi2/ndof Distribution");
  for( int i=0;i<nDatasets;i++)
  {
    trkChi2[i]->SetMarkerSize(msize);
    trkChi2[i]->SetLineWidth(msize);
    trkChi2[i]->SetMarkerStyle(mstyle[i]);
    trkChi2[i]->SetMarkerColor(mcol[i]);
    if( i==0) trkChi2[i]->Draw("PE");
    else trkChi2[i]->Draw("PESAME");
  }
  TLegend* leg3 = new TLegend(0.5, 0.7, 0.8, 0.85);
  for( int i=0;i<nDatasets;i++)
    leg3->AddEntry(trkChi2[i],dataName[i],"p");
  leg3->SetBorderSize(0);
  leg3->SetFillColor(0);
  leg3->Draw("SAME");

  c3->cd(2);
  trkNhit[0]->SetTitle("Track Nhit Distribution");
  for( int i=0;i<nDatasets;i++)
  {
    trkNhit[i]->SetMarkerSize(msize);
    trkNhit[i]->SetLineWidth(msize);
    trkNhit[i]->SetMarkerStyle(mstyle[i]);
    trkNhit[i]->SetMarkerColor(mcol[i]);
    if( i==0) trkNhit[i]->Draw("PE");
    else trkNhit[i]->Draw("PESAME");
  }

  c3->cd(4);
  ratChi2[1]->SetTitle("Track Chi2/ndof Ratio");
  ratChi2[1]->SetAxisRange(0.7,1.3,"Y");
  for( int i=1;i<nDatasets;i++)
  {
    ratChi2[i]->SetMarkerSize(msize);
    ratChi2[i]->SetLineWidth(msize);
    ratChi2[i]->SetMarkerStyle(mstyle[i]);
    ratChi2[i]->SetMarkerColor(mcol[i]);
    if( i==0) ratChi2[i]->Draw("PE");
    else ratChi2[i]->Draw("PESAME");
  }
  TLegend* legR3 = new TLegend(0.14, 0.7, 0.8, 0.85);
  for( int i=1;i<nDatasets;i++)
    legR3->AddEntry(ratChi2[i],Form("%s / %s",dataName[i],dataName[0]), "p");
  legR3->SetBorderSize(0);
  legR3->SetFillColor(0);
  legR3->Draw("SAME");

  c3->cd(5);
  ratNhit[1]->SetTitle("Track Nhit Ratio");
  ratNhit[1]->SetAxisRange(0.7,1.3,"Y");
  for( int i=1;i<nDatasets;i++)
  {
    ratNhit[i]->SetMarkerSize(msize);
    ratNhit[i]->SetLineWidth(msize);
    ratNhit[i]->SetMarkerStyle(mstyle[i]);
    ratNhit[i]->SetMarkerColor(mcol[i]);
    if( i==0) ratNhit[i]->Draw("PE");
    else ratNhit[i]->Draw("PESAME");
  }

  c3->cd(3);
  trkEtaPhi[0]->Draw("colz");
  trkEtaPhi[0]->SetXTitle(Form("%s    #eta",dataName[0]));
  trkEtaPhi[0]->SetTitle(Form("Track #eta - #phi distribution (%s)",dataName[0])); 
  c3->cd(6);
  trkEtaPhi[1]->Draw("colz");
  trkEtaPhi[1]->SetXTitle(Form("%s    #eta",dataName[1]));
  trkEtaPhi[1]->SetTitle(Form("Track #eta - #phi distribution (%s)",dataName[1])); 


}


