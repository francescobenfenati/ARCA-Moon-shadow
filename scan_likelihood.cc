#include "TVirtualFitter.h"
#include "TFile.h"
#include "TLegend.h"
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
 
using namespace std;

// 2D Gauss fit function at [0,0]
Double_t gaus20(Double_t *x, Double_t *par) {

   Double_t bgx = 1.+par[6]*x[0] + par[7]*x[0]*x[0] + par[8]*x[0]*x[0]*x[0];
   Double_t bgy = 1.+par[3]*x[1] + par[4]*x[1]*x[1]+par[5]*x[1]*x[1]*x[1];
   //Double_t SizeMS = TMath::Pi()*0.26*0.26; //angular size of moon or sun
   Double_t sig = par[1]; // Angular resolution 
   //Double_t amp = SizeMS/(2.*TMath::Pi()*sig*sig);
   Double_t amp = par[2];
   //amp *= par[2];
   Double_t r1 = ((x[0])/sig);
   Double_t r2 = ((x[1])/sig);
   return par[0]*bgx*bgy*(1. - (amp*TMath::Pi()*0.26*0.26)/(2.*TMath::Pi()*sig*sig)*TMath::Exp(-0.5*(r1*r1+r2*r2)));
}

// 2D Background fit function
Double_t bkg(Double_t *x, Double_t *par) {

  Double_t bgx = 1.+par[4]*x[0] + par[5]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0];
  Double_t bgy = 1.+par[1]*x[1] + par[2]*x[1]*x[1]+par[3]*x[1]*x[1]*x[1];
  return par[0]*bgx*bgy;
  //return par[0]+par[1]*x[1]+par[2]*x[1]*x[1]+par[3]*x[1]*x[1]*x[1]+par[4]*x[0]+par[5]*x[0]*x[0]+par[6]*x[0]*x[0]*x[0];
}



int scan_likelihood() {
 

  string daz, dalt, fitinf0, sundist, MC_sun_dist, lik, len;
 int nHitFit;
 
 double sig = 0.5;
 
 double like_min =  0.;
 double like_max = 200.;
 double like_cut = like_min;
 double like_step = 10;

 int nh2;
 
 float size = 0.05; // label and title size for various axes

// ---------- define 2D histogram of event number distribution

 double MaxGrad = 6.;
 double Bin = 0.1;
 int Nbin2 = (2.*MaxGrad + Bin/2. ) / Bin;
 TH2D *h2_noshad = new TH2D("h2_noshad","",Nbin2,-MaxGrad,MaxGrad,Nbin2,-MaxGrad,MaxGrad);
 h2_noshad->SetBinErrorOption(TH2::kPoisson);
 h2_noshad->SetStats(kFALSE);
 h2_noshad->SetXTitle("Azimuth [deg]");
 h2_noshad->SetYTitle("Altitude [deg]");
 h2_noshad->GetXaxis()->SetTitleSize(size);
 h2_noshad->GetXaxis()->SetTitleOffset(0.9);  
 h2_noshad->GetXaxis()->SetTitleFont(62);  
 h2_noshad->GetXaxis()->SetLabelFont(62);  
 h2_noshad->GetXaxis()->SetLabelSize(size);
 h2_noshad->GetYaxis()->SetTitleSize(size);
 h2_noshad->GetYaxis()->SetTitleOffset(0.9);  
 h2_noshad->GetYaxis()->SetTitleFont(62);  
 h2_noshad->GetYaxis()->SetLabelFont(62);  
 h2_noshad->GetYaxis()->SetLabelSize(size);

 TH2D *h2_shad = new TH2D("h2_shad","",Nbin2,-MaxGrad,MaxGrad,Nbin2,-MaxGrad,MaxGrad);
 h2_shad->SetBinErrorOption(TH2::kPoisson);
 h2_shad->SetStats(kFALSE);
 h2_shad->SetXTitle("Azimuth [deg]");
 h2_shad->SetYTitle("Altitude [deg]");
 h2_shad->GetXaxis()->SetTitleSize(size);
 h2_shad->GetXaxis()->SetTitleOffset(0.9);
 h2_shad->GetXaxis()->SetTitleFont(62);
 h2_shad->GetXaxis()->SetLabelFont(62);
 h2_shad->GetXaxis()->SetLabelSize(size);
 h2_shad->GetYaxis()->SetTitleSize(size);
 h2_shad->GetYaxis()->SetTitleOffset(0.9);
 h2_shad->GetYaxis()->SetTitleFont(62);
 h2_shad->GetYaxis()->SetLabelFont(62);
 h2_shad->GetYaxis()->SetLabelSize(size);


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
 
 // ---------- loop over n_like  --------------------------

   int number = 0;
   double p_nomoon;

   // ---------- loop over ntuple muon events --------------------------
   double dtr = TMath::DegToRad();
   ifstream fstream;
 
   while(like_cut<like_max+like_step) {

 fstream.open("csv/arca8/MC/arca8rbr_mu_pre+anoise+combined_shadow.txt");
 {
      while ( ! fstream.eof() )
    {
      fstream >> daz >> dalt >> fitinf0 >> sundist >> MC_sun_dist >> lik >> len;
     double DeltaAzi = std::stof(daz);
     double DeltaAlt = std::stof(dalt);
     double beta0 = std::stof(fitinf0)*TMath::RadToDeg();
     double sun_dist = std::stof(sundist);
     double mc_sun_dist = std::stof(MC_sun_dist);
     double likelihood = std::stof(lik);
     double track_length = std::stof(len);
    
     if ( likelihood > like_cut ) // <------------------------------------------------------ CUT ----------------
       { 
	 h2_noshad->Fill(DeltaAzi,DeltaAlt);
	 if (mc_sun_dist > 0.26) { // <-------------------------------PREPARE S1 SAMPLE----------
	   if (DeltaAzi < MaxGrad &&  DeltaAzi > -MaxGrad && DeltaAlt < MaxGrad && DeltaAlt > -MaxGrad) {
	     hAzi->Fill(DeltaAzi,1./Nbin2);
	   }
	   if ( DeltaAlt < MaxGrad && DeltaAlt > -MaxGrad && DeltaAzi < MaxGrad &&  DeltaAzi > -MaxGrad) {
	     hAlt->Fill(DeltaAlt,1./Nbin2);
	   }
	   h2_shad->Fill(DeltaAzi,DeltaAlt); 
	 }
       }

    } // end event loop
 fstream.close();
 }
   

  // --------------- FITTING ------------
 cout<<"like_cut :  "<<like_cut<<endl;

 //cout<<"N. of events in S0 = "<<h2_noshad->GetEntries()<<endl;
 //cout<<"N. of events in S1 = "<<h2_shad->GetEntries()<<endl;
 //cout<<"N. of events in x distribution = "<<hAzi->GetEntries()<<endl;
 //cout<<"N. of events in y distribution = "<<hAlt->GetEntries()<<endl;


 TF1 *pol0 = new TF1("pol0","pol0",-MaxGrad,MaxGrad);
 TF1 *pol2 = new TF1("pol2", "pol2",-MaxGrad,MaxGrad);
 TF1 *pol3 = new TF1("pol3", "pol3",-MaxGrad,MaxGrad);

 hAzi->Fit("pol0","QN");
 double a0 = pol0->GetParameter(0);
 double x1 = 0;
 double x2 = 0;
 double x3 = 0;
 
 
 TVirtualFitter *fitterazi = TVirtualFitter::GetFitter();
 Double_t aminazi,aminazi_l,chi2azi,chi2azi_l;
 Double_t amin,edm,errdef;
 Int_t nvpar,nparx;

 //fitterazi->PrintResults(2,0.);
 fitterazi->GetStats(aminazi,edm,errdef,nvpar,nparx);
 chi2azi = aminazi; //<------------ no need of multiplying by 2 as the minimum fcn value is already the pearson chi2 value since fit is by minimizing chi2!! 
 //cout << "chi2azi = "<< chi2azi;  

 hAlt->Fit("pol2","N");
 a0 = pol2->GetParameter(0);
 double y1 = pol2->GetParameter(1)/a0;
 double y2 = pol2->GetParameter(2)/a0;
 double y3 = 0;
 
 double chi2alt = pol2->GetChisquare();
 cout << "--------- 1D FITS --------  "<<endl;
 cout << "--------------------------  "<<endl;
 cout << "chi2azi_ = " << chi2azi<<endl;
 cout << "chi2alt = " << chi2alt<<endl;
 cout << "p(azi) = " << TMath::Prob(chi2azi,Nbin2-(4))<<endl;
 cout << "p(alt) = " << TMath::Prob(chi2alt,Nbin2-(4))<<endl;



 gROOT->SetStyle("Plain");
 gStyle->SetOptFit(1111);
 gStyle->SetOptStat("nemr");

 
 const Int_t npar = 11;
 double amp = 1.;
 double SigFit = sig; 
    

    const Int_t np20 = 9;
    Double_t f20params[np20] = {a0,SigFit,amp,y1,y2,y3,x1,x2,x3};
    TF2 *f20 = new TF2("f20",gaus20,-MaxGrad,MaxGrad,-MaxGrad,MaxGrad, np20);
    f20->SetParameters(f20params);
    //f20->SetParameters(100,0.1,1.,0.001,0.001,0.001,0.001,0.001,0.001);
    f20->SetParNames("a0","sigma","amp","y1","y2","y3","x1","x2","x3");    
    f20->SetParLimits(1,0.1,1.00);
    f20->SetParLimits(2,0.1,2.00);

    const Int_t npbg = 7;
    Double_t b2params[npbg] = {a0,y1,y2,y3,x1,x2,x3};
    TF2 *bg2 = new TF2("bg2",bkg,-MaxGrad,MaxGrad,-MaxGrad,MaxGrad, npbg);
    bg2->SetParameters(b2params);
    //bg2->SetParameters(100,0.01,0.001,0.001,0.001,0.001,0.001);
    bg2->SetParNames("a0","y1","y2","y3","x1","x2","x3");



   // ---------- 2Dim Fit on h2 histogram of event distribution
        
    // fit of MS at nominal position
    h2_shad->Fit("f20","NL");
    double amp0 = f20->GetParameter(2);
    double Sig0 = f20->GetParameter(1);

    TVirtualFitter *fitter20 = TVirtualFitter::GetFitter();
    Double_t amin20;
    double chi2Fit0;
    //fitter20->PrintResults(2,0.);
    fitter20->GetStats(amin20,edm,errdef,nvpar,nparx);
    chi2Fit0 = amin20*2;

    // background only fit
    
    h2_shad->Fit("bg2","QNL");
    double chi2Fitbg;
    TVirtualFitter *fitterbg = TVirtualFitter::GetFitter();
    Double_t aminbg;
    fitterbg->GetStats(aminbg,edm,errdef,nvpar,nparx);
    chi2Fitbg = aminbg*2;

   // delta chi2
   double lamb = chi2Fit0 - chi2Fitbg;
   
   cout << "PARAMETERS for MC : ";
   //cout << "Events selected: "<<nh2<<endl;
   cout << "--------- NOMINAL FIT (xs=0, ys=0) --------  "<<endl;
   cout << "Sigma: "<< Sig0<<"  with error: "<<f20->GetParError(1)<<endl; // comment  if sigma fixed for nominal fit
   cout << "amp(0,0) = " << amp0 <<"with error: "<<f20->GetParError(2)<< endl;
   cout << "lambda = " << lamb <<endl;
   
   ofstream myfile("csv/arca8/arca8rbr_mu_pre+anoise+moon+combined_shadow_likes.txt",ios::app);
   if(myfile){
     myfile << like_cut <<" " << lamb <<" "<< Sig0 <<" "<< amp0 << endl;
   }
   
   like_cut+=like_step;

   h2_shad->Reset();
   h2_noshad->Reset();
   hAzi->Reset();
   hAlt->Reset(); 

   myfile.close();
   }
 cout<<"finished"<<endl;
 return 0;
 }
