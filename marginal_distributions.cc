//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//The script asks for a sigm res. value, obtained from the 1D analysis. However by default it uses the value only as initial par. value for the 2D hist fit, and produces the lambda plot by fitting with the sigma res. found with the 2D hist. (i.e. Sig0). If you want to use the sigma res. from the 1D analyses, you need to change that line in the fit.
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "TVirtualFitter.h"
#include "TFile.h"
#include "TLegend.h"
#include <TMarker.h>
#include "TApplication.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TColor.h"
#include <chrono> 

using namespace std;
using namespace std::chrono;

// 2D Gauss fit function
Double_t gaus2(Double_t *x, Double_t *par) {

   Double_t bkg_x = 1.+par[8]*x[0] + par[9]*x[0]*x[0] + par[10]*x[0]*x[0]*x[0];
   Double_t bkg_y = 1.+par[5]*x[1] + par[6]*x[1]*x[1]+par[7]*x[1]*x[1]*x[1];
   Double_t sig = par[3]; // Angular resolution
   Double_t amp = par[4];
   Double_t r1 = ((x[0]-par[1])/sig);
   Double_t r2 = ((x[1]-par[2])/sig);
   return par[0]*bkg_x*bkg_y*(1. - (amp*TMath::Pi()*0.26*0.26)/(2.*TMath::Pi()*sig*sig)*TMath::Exp(-0.5*(r1*r1+r2*r2)));
   //return par[0]*bkg_x*bkg_y*(1. - amp/(2.*TMath::Pi()*sig*sig)*TMath::Exp(-0.5*(r1*r1+r2*r2)));
}

// 2D Gauss fit function at [0,0]
Double_t gaus20(Double_t *x, Double_t *par) {

   Double_t bkg_x = 1.+par[6]*x[0] + par[7]*x[0]*x[0] + par[8]*x[0]*x[0]*x[0];
   Double_t bkg_y = 1.+par[3]*x[1] + par[4]*x[1]*x[1]+par[5]*x[1]*x[1]*x[1];
   //Double_t SizeMS = TMath::Pi()*0.26*0.26; //angular size of moon or sun
   Double_t sig = par[1]; // Angular resolution 
   //Double_t amp = SizeMS/(2.*TMath::Pi()*sig*sig);
   Double_t amp = par[2];
   //amp *= par[2];
   Double_t r1 = ((x[0])/sig);
   Double_t r2 = ((x[1])/sig);
   return par[0]*bkg_x*bkg_y*(1. - (amp*TMath::Pi()*0.26*0.26)/(2.*TMath::Pi()*sig*sig)*TMath::Exp(-0.5*(r1*r1+r2*r2)));
   //return par[0]*bkg_x*bkg_y*(1. - amp/(2.*TMath::Pi()*sig*sig)*TMath::Exp(-0.5*(r1*r1+r2*r2)));
}

// 2D Background fit function
Double_t bkg(Double_t *x, Double_t *par) {

  Double_t bkg_x = 1.+par[4]*x[0] + par[5]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0];
  Double_t bkg_y = 1.+par[1]*x[1] + par[2]*x[1]*x[1]+par[3]*x[1]*x[1]*x[1];
  return par[0]*bkg_x*bkg_y;
  //return par[0]+par[1]*x[1]+par[2]*x[1]*x[1]+par[3]*x[1]*x[1]*x[1]+par[4]*x[0]+par[5]*x[0]*x[0]+par[6]*x[0]*x[0]*x[0];
}



int main(int argc, char** argv) {
  
  if (argc < 3) {
    cout << "3 arguments are required: <input_file> <beta_cut> <max_deg> "<<endl;
    exit(1);
  }

  
  auto start = high_resolution_clock::now();


  //TApplication app("app",NULL,NULL); 
  TFile *file=new TFile("data_sun.root", "RECREATE");
  
  string filename, daz, dalt, fitinf0, sundist, MC_sun_dist, lik, len, time, run;
  double chi;
  int nHitFit;
 

 double beta_cut = std::stof(argv[2]);

 float size = 0.05; // label and title size for various axes


 double MaxGrad = std::stof(argv[3]);
 double Bin = 0.3;
 int Nbin2 = (2.*MaxGrad + Bin/2. ) / Bin;

// ---------- define 1D projection of event distribution to obtain BG Altitude fit

 TH1D *hAlt = new TH1D("hAlt","",Nbin2,-MaxGrad,MaxGrad);
 //hAlt->SetStats(kFALSE);
 hAlt->SetXTitle("Altitude [deg]");
 hAlt->SetYTitle("Events");
 hAlt->GetXaxis()->SetTitleSize(size);
 hAlt->GetXaxis()->SetTitleOffset(0.9);  
 hAlt->GetXaxis()->SetTitleFont(62);  
 hAlt->GetXaxis()->SetLabelFont(62);  
 hAlt->GetXaxis()->SetLabelSize(size);
 hAlt->GetYaxis()->SetTitleSize(size);
 hAlt->GetYaxis()->SetTitleOffset(0.9);  
 hAlt->GetYaxis()->SetTitleFont(62);  
 hAlt->GetYaxis()->SetLabelFont(62);  
 hAlt->GetYaxis()->SetLabelSize(size);

// ---------- define 1D projection of event distribution to obtain BG Azimuth fit

 TH1D *hAzi = new TH1D("hAzi","",Nbin2,-MaxGrad,MaxGrad);
 //hAzi->SetStats(kFALSE);
 hAzi->SetXTitle("Azimuth [deg]");
 hAzi->SetYTitle("Events");
 hAzi->GetXaxis()->SetTitleSize(size);
 hAzi->GetXaxis()->SetTitleOffset(0.9);  
 hAzi->GetXaxis()->SetTitleFont(62);  
 hAzi->GetXaxis()->SetLabelFont(62);  
 hAzi->GetXaxis()->SetLabelSize(size);
 hAzi->GetYaxis()->SetTitleSize(size);
 hAzi->GetYaxis()->SetTitleOffset(0.9);  
 hAzi->GetYaxis()->SetTitleFont(62);  
 hAzi->GetYaxis()->SetLabelFont(62);  
 hAzi->GetYaxis()->SetLabelSize(size);
 
 // ---------- loop over n_beta  --------------------------

 cout << "Analyzing: "<<argv[1]<<endl;
 cout << "Beta0 cut = "<<argv[2]<<endl;
 cout << "Maxdeg = "<<argv[3]<<endl;

 int n_h1 =0;
 double p_nomoon;
 double chi2Fit0;
 // ---------- loop over ntuple muon events --------------------------
 double dtr = TMath::DegToRad();
 ifstream fstream;
   
   fstream.open(argv[1]);
   

 {
      while ( ! fstream.eof() )
    {
      fstream >> daz >> dalt >> fitinf0 >> sundist >> lik >> len >> time ;//>> run;
     double DeltaAzi = std::stof(daz);
     double DeltaAlt = std::stof(dalt);
     double beta0 = std::stof(fitinf0)*TMath::RadToDeg();
     double sun_dist = std::stof(sundist);
     double likelihood = std::stof(lik);
     double track_length = std::stof(len);
     double time_since_epoch = std::stoi(time);
     //int run_id = std::stoi(run);
     if ( beta0 < beta_cut ) // <------------------------------------------------------ CUT ----------------
       
       { 
	 if (DeltaAzi < MaxGrad &&  DeltaAzi > -MaxGrad && DeltaAlt < MaxGrad && DeltaAlt > -MaxGrad ) {
	     hAzi->Fill(DeltaAzi,1./Nbin2);

	 }
	 if (DeltaAzi < MaxGrad &&  DeltaAzi > -MaxGrad && DeltaAlt < MaxGrad && DeltaAlt > -MaxGrad ) {
	     hAlt->Fill(DeltaAlt,1./Nbin2);
	 }
	 
       }

    }
 }// end event loop
 fstream.close();
  
  // --------------- event loop finished ------------ make fits !! 

 cout<<"N. of events in x distribution = "<<hAzi->GetEntries()<<endl;
 cout<<"N. of events in y distribution = "<<hAlt->GetEntries()<<endl;


   // fit projection of 2D histogram to correct for changing event rate
 TF1 *pol0 = new TF1("pol0", "pol0",-MaxGrad,MaxGrad);
 TF1 *pol2 = new TF1("pol2", "pol2",-MaxGrad,MaxGrad);
 TF1 *pol3 = new TF1("pol3", "pol3",-MaxGrad,MaxGrad);

 cout<<"Fit minchi2 Azimuth: "<<endl;
 hAzi->Fit("pol0","N");
 double a0 = pol0->GetParameter(0);
 double x1 = 0;
 double x2 = 0;
 double x3 = 0;
 //double x1 = pol3->GetParameter(1)/a0;
 //double x2 = pol3->GetParameter(2)/a0;
 //double x3 = pol3->GetParameter(3)/a0;
 
 double chi2azi = pol0->GetChisquare(); //<------------ no need of multiplying by 2 as the minimum fcn value is already the pearson chi2 value!! 
 


 hAlt->Fit("pol2","N");
 a0 = pol2->GetParameter(0);
 double y1 = pol2->GetParameter(1)/a0;
 double y2 = pol2->GetParameter(2)/a0;
 double y3 = 0; 
//double y3 = pol3->GetParameter(3)/a0;
 double chi2alt = pol2->GetChisquare();
  
 cout << "--------- 1D FITS --------  "<<endl;
 cout << "--------------------------  "<<endl;
 cout << "chi2azi_ = " << chi2azi<<endl;
 cout << "chi2alt = " << chi2alt<<endl;
 cout << "p(azi) = " << TMath::Prob(chi2azi,Nbin2-(4))<<endl;
 cout << "p(alt) = " << TMath::Prob(chi2alt,Nbin2-(4))<<endl;



// define 2d fit function, initialize parameters

 //gROOT->SetStyle("Plain");
 gStyle->SetOptFit(1111);
 gStyle->SetOptStat("nemr");


 //--------------PLOT MARGINAL DISTRIBUTIONS-------------------
 TCanvas *c0 = new TCanvas("c0", "Azimuth", 200, 10, 600, 400);
 c0->cd();
 hAzi->Draw("H");
 hAzi->Draw("E, P, SAME");
 pol0->Draw("SAME");
 c0->SaveAs("/sps/km3net/users/fbenfe/Moon_shadow/data_azimuth_distr.pdf");
 

 TCanvas *c1 = new TCanvas("c1", "Altitude", 200, 10, 600, 400);
 c1->cd();
 hAlt->Draw("E, P, SAME");
 hAlt->Draw("H");
 pol2->Draw("SAME");
 
 /*
 gPad->Update();
 auto stat = dynamic_cast<TPaveStats*>(hAlt->FindObject("stats"));
 if (stat) {
   //cout << " X1NDC: " << stat->GetX1NDC() << " X2NDC: " << stat->GetX2NDC() << endl;
   stat->SetX1NDC(0.1); stat->SetX2NDC(0.3);
   stat->Draw();
 } else {
   cout << "No stats box found!\n";
 }
 */
 c1->SaveAs("/sps/km3net/users/fbenfe/Moon_shadow/data_altitude_distr.pdf");

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<minutes>(stop - start);
    cout << "-------------------------------------"<<endl;
    cout << "Time elapsed: "
	 << duration.count() << " minutes" << endl;
    cout << "-------------------------------------"<<endl;
    
    //app.Run(kTRUE); 
    cout<<"finished"<<endl;
    return 0;
}

