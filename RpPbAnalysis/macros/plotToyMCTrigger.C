#include "trackingPlotUtilities.C"

plotToyMCTrigger()
{

  double ptbins[38] = {        0.0, 0.1, 0.2, 0.3, 0.4, 
        0.5, 0.6, 0.7, 0.8, 0.9, 
        1.1, 1.2, 1.4, 1.6, 1.8, 
        2.0, 2.2, 2.4, 3.2, 4.0, 
        4.8, 5.6, 6.4, 7.2, 9.6, 
        12.0, 14.4, 19.2, 24.0, 28.8, 
        35.2, 41.6, 48.0, 60.8, 73.6, 
        86.4, 103.6, 130.0 };

  double eff[38] = { 1.000, 0.027, 0.573, 0.693, 0.790,
     0.856, 0.862, 0.901, 0.902, 0.903,
     0.903, 0.902, 0.901, 0.900, 0.899,
     0.899, 0.897, 0.892, 0.888, 0.888,
     0.891, 0.889, 0.890, 0.894, 0.895,
     0.898, 0.893, 0.896, 0.890, 0.878,
     0.875, 0.863, 0.858, 0.836, 0.824,
     0.797, 0.772 }; 
    
   
    

  TFile * f = new TFile("toyMCTrigger_PythiaZ2star_ptHat15_100M.root");
  //TFile * f = new TFile("toyMCTrigger_PythiaZ2star_ptHat15.root");
  //TFile * f = new TFile("pythiaZ2_pthat15_100.root");

  // array goes as 0=MB, 1=Trk12, 2=Trk20, 3=Trk30
  
  TH1F *hGen;
  TH1F *hTrk[4];
  TH1F *hTrk12;

  hGen = (TH1F *) f->Get("trkAna/genSpectrum"); 

  hTrk[0] = (TH1F*) f->Get("trkAna/mbSpectrumSel");
  hTrk[1] = (TH1F*) f->Get("trkAna/trk12SpectrumSel");
  hTrk[2] = (TH1F*) f->Get("trkAna/trk20SpectrumSel");
  hTrk[3] = (TH1F*) f->Get("trkAna/trk30SpectrumSel");

  hTrk12 = (TH1F*) f->Get("trkAna/trk12Spectrum");

  TH1F * hgenevt;
  TH1F * hnevt[4];
 
  hgenevt = (TH1F *) f->Get("trkAna/genLeads");
  hnevt[0] = (TH1F *) f->Get("trkAna/mbLeads");
  hnevt[1] = (TH1F *) f->Get("trkAna/trk12Leads");
  hnevt[2] = (TH1F *) f->Get("trkAna/trk20Leads");
  hnevt[3] = (TH1F *) f->Get("trkAna/trk30Leads");
  
  // Calculate Nevt for each sample

  double genevt;
  double nevt[4];
  double nevtTrk12;
    
  // Adding the nevent histogram
    
  genevt = hgenevt->Integral();
  

  nevt[0] = hnevt[0]->Integral( hnevt[0]->GetXaxis()->FindBin(0.01),
                               hnevt[0]->GetXaxis()->FindBin(999.9) );
    
  nevt[1] = nevt[0] * hnevt[1]->Integral( hnevt[1]->GetXaxis()->FindBin(12.1),
                                         hnevt[1]->GetXaxis()->FindBin(19.1) )
            / hnevt[0]->Integral( hnevt[0]->GetXaxis()->FindBin(12.1),
                                  hnevt[0]->GetXaxis()->FindBin(19.1) );
  
  nevt[2] = nevt[1] * hnevt[2]->Integral( hnevt[2]->GetXaxis()->FindBin(19.3),
                                         hnevt[2]->GetXaxis()->FindBin(28.7) )
            / hnevt[1]->Integral( hnevt[1]->GetXaxis()->FindBin(19.3),
                                  hnevt[1]->GetXaxis()->FindBin(28.7) );
  nevt[3] = nevt[2] * hnevt[3]->Integral( hnevt[3]->GetXaxis()->FindBin(28.9),
                                          hnevt[3]->GetXaxis()->FindBin(999.9) )
            / hnevt[2]->Integral( hnevt[2]->GetXaxis()->FindBin(28.9),
                                  hnevt[2]->GetXaxis()->FindBin(999.9) );

  nevtTrk12 = nevt[0] * hnevt[1]->Integral( hnevt[1]->GetXaxis()->FindBin(12.1),
                                         hnevt[1]->GetXaxis()->FindBin(999.9) )
            / hnevt[0]->Integral( hnevt[0]->GetXaxis()->FindBin(12.1),
                                  hnevt[0]->GetXaxis()->FindBin(999.9) );

  cout << "Nevt Gen :   " << genevt << endl;
  cout << "Nevt MB :    " << nevt[0] << endl;
  cout << "Nevt Trk12 : " << nevt[1] << endl;
  cout << "Nevt Trk20 : " << nevt[2] << endl;
  cout << "Nevt Trk30 : " << nevt[3] << endl;
  cout << endl;

  double error;
  cout << "MB 12-20: " << hnevt[0]->IntegralAndError( hnevt[0]->GetXaxis()->FindBin(12.1),
                                  hnevt[0]->GetXaxis()->FindBin(19.1), error )  << " +/- ";
  cout << error << endl;   
  cout << "Trk12 20-30: " << hnevt[1]->IntegralAndError( hnevt[1]->GetXaxis()->FindBin(19.3),
                                  hnevt[1]->GetXaxis()->FindBin(28.7), error  ) << " +/- ";
  cout << error << endl; 
  cout << "Trk20 30-inf: " << hnevt[1]->IntegralAndError( hnevt[1]->GetXaxis()->FindBin(28.9),
                                  hnevt[1]->GetXaxis()->FindBin(999.9), error  ) << " +/- "; 
  cout << error << endl; 
 
/////////////////// Applying efficiency corrections first ///////////////////
    
  for( int i=0;i<4;i++)
  {
    for(int x=1; x <= hTrk[i]->GetNbinsX(); x++)
    {
      double xpos =  hTrk[i]->GetXaxis()->GetBinCenter(x);
      double xeff = 1.0;
      for(int j=0;j<37;j++)
      { 
        if( xpos > ptbins[j] && xpos <= ptbins[j+1] ) xeff = eff[j];
      }
      hTrk[i]->SetBinContent(x, hTrk[i]->GetBinContent(x)/xeff);
      hTrk[i]->SetBinError(x, hTrk[i]->GetBinError(x)/xeff);
    }    

  }

  for(int x=1; x <= hTrk12->GetNbinsX(); x++)
  {
    double xpos =  hTrk12->GetXaxis()->GetBinCenter(x);
    double xeff = 1.0;
    for(int j=0;j<37;j++)
    {
      if( xpos > ptbins[j] && xpos <= ptbins[j+1] ) xeff = eff[j];
    }
    hTrk12->SetBinContent(x, hTrk12->GetBinContent(x)/xeff);
    hTrk12->SetBinError(x, hTrk12->GetBinError(x)/xeff);
  }

/////////////////////////////////////////////////////////////////////////////
    

  normalizeBinAndEtaWidth( hGen, 2.0 );
  hGen->Scale(1./genevt);
  normalizeBinAndEtaWidth( hTrk12, 2.0 );
  hTrk12->Scale(1./nevtTrk12);

  for( int i=0; i<4; i++)
  {
    normalizeBinAndEtaWidth( hTrk[i], 2.0 );
    hTrk[i]->Scale(1./nevt[i]);
  }
     
  // add together spectra
     
  for( int i=1;i<4;i++)
  {
     hTrk[0]->Add(hTrk[i],1.);
  }
    
  ////////////////////////////////////////////////////////////
    
  hTrk[0]->SetMarkerStyle(20);
  hGen->SetMarkerStyle(21);
  hGen->SetMarkerColor(kBlue);
  hTrk12->SetMarkerStyle(24);
  hTrk12->SetMarkerColor(kRed);
   TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   //c1->SetLogy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.13);
   c1->SetRightMargin(0.06);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.16);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);

   //c1->Divide(1,2,0,0);
   c1->cd();
   c1->SetLogy();
   c1->SetLogx();
   c1->SetLogx();
  // c1->SetLogy();
  // c1->SetLogx();

   TH1F * hDum = new TH1F("hhdum","",10,0.4,110.);
   hDum->SetMinimum(1e-10);
   hDum->SetMaximum(100);
   hDum->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   hDum->GetXaxis()->SetTitleSize(0.035);
   hDum->GetXaxis()->SetTitleOffset(1.5);
   hDum->GetXaxis()->SetLabelSize(0.035);
   hDum->GetYaxis()->SetTitle("1/N_{evt} d^{2}N_{ch}/d_{T}d#eta");
   hDum->GetYaxis()->SetTitleSize(0.035);
   hDum->GetYaxis()->SetLabelSize(0.035);
   hDum->GetYaxis()->SetTitleOffset(1.5);

  hDum->Draw();
  hGen->Draw("same");
  hTrk12->Draw("same");
  hTrk[0]->Draw("same");

  TLatex * tex = new TLatex(0.6,1.0e-5,"|#eta| < 1.0");
  tex->SetTextSize(0.036);
  tex->Draw();
  TLatex * tex = new TLatex(0.6,1.0e-4,"PYTHIA Z2* p_{T}-hat = 15 GeV/c");
  tex->SetTextSize(0.036);
  tex->Draw();


  TLegend * leg = new TLegend(0.22,0.18,0.8,0.4);
  leg->AddEntry(hGen,"Generated Particles","lp");
  leg->AddEntry(hTrk[0],"Corrected Tracks","lp");
  leg->AddEntry(hTrk12,"Corrected Tracks (Trk12 only)","lp");
 // leg->AddEntry(hTrkPtRev[0],"Corrected Tracks Pbp","lp");
  leg->SetFillColor(kWhite);
  leg->Draw();

   TCanvas *c2 = new TCanvas("c2a", "c2",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

   c2->Range(0,0,1,1);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   //c2->SetLogy();
   c2->SetTickx(1);
   c2->SetTicky(1);
   c2->SetLeftMargin(0.13);
   c2->SetRightMargin(0.06);
   c2->SetTopMargin(0.05);
   c2->SetBottomMargin(0.16);
   c2->SetFrameFillStyle(0);
   c2->SetFrameBorderMode(0);

   c2->cd();
   c2->SetLogx();

   TH1F * hDum2 = new TH1F("hhdum2","",10,0.4,110.);
   hDum2->SetMinimum(0.5);
   hDum2->SetMaximum(1.5);
   hDum2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   hDum2->GetXaxis()->SetTitleSize(0.035);
   hDum2->GetXaxis()->SetTitleOffset(1.5);
   hDum2->GetXaxis()->SetLabelSize(0.035);
   hDum2->GetYaxis()->SetTitle("ratio reco/gen");
   hDum2->GetYaxis()->SetTitleSize(0.035);
   hDum2->GetYaxis()->SetLabelSize(0.035);
   hDum2->GetYaxis()->SetTitleOffset(1.5);

   hDum2->Draw();

   TH1F * hRat = hTrk[0]->Clone("hRat");
   hRat->Divide(hGen);
   hRat->Draw("same");
   TH1F * hRat12 = hTrk12->Clone("hRat12");
   hRat12->Divide(hGen);
   hRat12->Draw("same");

}
