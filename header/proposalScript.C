# include <iostream>
# include <TFile.h>
# include <TPad.h>
# include <TLegend.h>
# include <TTree.h>
# include <TCut.h>
# include <TH1F.h>
# include "TF1.h"
# include "TH1D.h"
# include <TH1.h>
# include <TH2.h>

using namespace std;



//===========================================FITTING==================================================
	Double_t background(Double_t *x, Double_t *par) {
	  return (par[0] + par[1]*x[0] + par[2]*x[0]*x[0]);
	}

	Double_t backgroundExp(Double_t *x, Double_t *par) {
	  return (par[0] + TMath::Exp(-(par[1] + par[2]*x[0])));
	}

	Double_t GaussianPeak(Double_t *x, Double_t *par)
	// The signal function: a gaussian
	{
	  
  Double_t arg = 0;
  if (par[2]<0) par[2]=-par[2];  // par[2]: sigma
  if (par[2] != 0) arg = (x[0] - par[1])/par[2];  // par[1]: mean
 
 //return par[0]*BIN_SIZE*TMath::Exp(-0.5*arg*arg)/
  //   (TMath::Sqrt(2*TMath::Pi())*par[2]); 
   return (par[0]*TMath::Exp(-0.5*arg*arg)/
     (TMath::Sqrt(2*TMath::Pi())*par[2])); // par[0] is constant
 

	}


	Double_t Gauss(Double_t *x, Double_t *par)
	// The signal function: a gaussian
	{
	  
  Double_t arg = 0;
  if (par[2]<0) par[2]=-par[2];  // par[2]: sigma
  if (par[2] != 0) arg = (x[0] - par[1])/par[2];  // par[1]: mean
 
 //return par[0]*BIN_SIZE*TMath::Exp(-0.5*arg*arg)/
  //   (TMath::Sqrt(2*TMath::Pi())*par[2]); 
   return (par[0]*TMath::Exp(-0.5*arg*arg)); // par[0] is constant
 

	}

	// Lorenzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) /
    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
   + .25*par[1]*par[1]);

}

	// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  return  lorentzianPeak(x,&par[3]) + background(x,par);
}

// Sum of background and peak function
Double_t fitFunction2(Double_t *x, Double_t *par) {
  return backgroundExp(x,par) + Gauss(x,&par[3]);
}

// combine function with gaus and polynomial
Double_t the_gausppar(double* x, double* par){

     return (par[0]*TMath::Gaus(x[0],par[1],par[2])+
         par[3]+par[4]*x[0]+par[5]*x[0]*x[0]);
     }

void format_line(TAttLine* line,int col,int sty){
	
     line->SetLineWidth(5); line->SetLineColor(col);
     line->SetLineStyle(sty);
 }


void proposalScript()
{
  

    // gStyle->SetOptStat(0);
     gStyle->SetOptFit(1100);
  //	gStyle->SetOptStat("e");

	// gStyle->SetTitleFontSize(0.1);
	// data before missing mass and momentum cut,including timing and vertex cut 
	//TFile *file = TFile::Open("kpppim_kinfit_old.root");
	//TTree *tf = (TTree*)file->Get("output");
    TFile *file = TFile::Open("kpppim_kinfit.root");
	TTree *tf = (TTree*)file->Get("output");
  
	
	TFile *outFile = new TFile("histoKppmpim.root", "RECREATE");
  	TCut mmCut = "abs(sqrt(mm2pkp) - 0.139)<0.08 ";
  	//TCut IM = "abs(sqrt(invm2ppim) - 1.115)<0.04";
  	TCut mmKP = "sqrt(mm2kp) < 1.16";
 	//TCut mm = "sqrt(mm2pkp) < 0.22";  // fitting is implemented here
 	TCut mm = "abs(sqrt(mm2pkp) - 0.1311)< 0.036";
 	TCut nPlus = "trackPositive == 2";
 	TCut nNeu = "trackNeutral == 0";

	/*
	// event selection based on missing mass, missing momentum and time of flight cut only
     
	TCanvas *c = new TCanvas("c" ,"", 20, 20, 1000, 1000);
	c -> Divide(2,2);
  
	TH1F *g = new TH1F("g", "", 100, 1.0, 1.5);
	TH1F *g1 = new TH1F("g1", "", 100, 0.2,0.8);
	TH1F *g2 = new TH1F("g2", "", 100, 1.0, 1.5);
	TH1F *g3 = new TH1F("g3", "",100, 0.0, 0.6);
    
	c -> cd(1);
	tf -> Draw("sqrt(invm2ppim) >>g");
	g -> SetFillColor(3);
	g -> GetXaxis() -> SetTitle("mass of proton-pion");
	c -> cd(2);
	tf -> Draw("sqrt(mm2ppim) >>g1");
	g1 -> SetFillColor(3);
	g1 -> GetXaxis() -> SetTitle("mass of kaon");
	c -> cd(3);
	tf -> Draw("sqrt(mm2kp) >>g2");
	g2 -> SetFillColor(3);
	g2 -> GetXaxis() -> SetTitle("massing mass of kaon");
	c -> cd(4);
	tf ->Draw("sqrt(mm2pkp)>>g3");
	g3 -> SetFillColor(3);
	g3 -> GetXaxis() -> SetTitle("mass of pion");
	
	// event selection based on missing mass, missing momentum, tof, tofKO, fiducial cut, the only difference is mass selection cut of kaon "kcut"
	// bunch of 2d plot
  
	TCanvas *cp = new TCanvas("cp" ,"", 20, 20, 1000, 1000);
	cp -> Divide(2,2);
	
	TH2F *gc = new TH2F("gc", "", 100, 1.0, 1.5,100,1.0,1.5);
	TH2F *g1c = new TH2F("g1c", "", 100, 1.0,1.9,100, 0.0,0.6);
	TH2F *g2c = new TH2F("g2c", "", 100, 0.5, 5.6,100,0.8,1.4);
	TH1F *g3c = new TH1F("g3c", "",100, 0.1, 2.0);
    
	cp -> cd(2);
	tf -> Draw("sqrt(mm2kp):sqrt(invm2ppim) >>gc");
	gc -> GetYaxis() -> SetTitle("missing mass of kaon");
	gc -> GetXaxis() -> SetTitle("invariant mass of lambda");
	gc -> Draw("colz");

	 // measured mass of kaon and proton
	cp -> cd(1);
	TH2F *h1 = new TH2F("h1", " ", 200, 0.35, 0.7, 200, 0.8,1.4);
	h1 -> GetYaxis()-> SetTitle("proton mass");
	h1 -> GetXaxis()-> SetTitle("kaon mass");
	
	tf -> Draw("sqrt(pow(pp,2)*(1 - pow(tofbetap,2)))/(tofbetap ):sqrt(pow(kpp,2)*(1 - pow(tofbetakp,2)))/(tofbetakp ) >> h1"); 
     
	h1 -> Draw("colz");
	
	cp -> cd(3);
	tf -> Draw("mm2ppim:mm2kp >> g1c");
	g1c -> GetYaxis() -> SetTitle("missing mass square of proton-pion");
	g1c -> GetXaxis() -> SetTitle("missing mass square of kaon");
	g1c -> Draw("colz");

	cp -> cd(4);
	tf -> Draw("sqrt(mm2kp):ebeam >> g2c");
	g2c -> GetYaxis() -> SetTitle(" missing mass of kaon");
	g2c -> GetXaxis() -> SetTitle("ebeam");
	g2c -> Draw("colz");
    
	/*
	// script for inset plot of invariant mass of proton and pion
  

	TH1F *h1 = new TH1F("h1", "", 100, 1.09,1.16);
	TH1F *h = new TH1F("h", "", 100, 1.09,1.5);
  

	tT -> Draw("mppim >>h1",tof&&fid&&TOFKO);
	h1 -> GetYaxis()-> SetTitle("Events");
	h1 -> GetXaxis()-> SetTitle("Mass(p#pi^{-}) GeV/c^{2}");
	tf -> Draw("mppim >>h");
	TH1F *h2 = (TH1F*)h -> Clone();
  
	// h2 -> GetXaxis()-> SetLabelOffset(999);
	//  h2 -> GetXaxis() -> SetLabelSize(0);
 
  
	//  h -> GetYaxis()-> SetTitle("No. of events");
	//  h -> GetXaxis()-> SetTitle("Mass(p#pi^{-})");
   
  
	TCanvas *c2p = new TCanvas("c2p","",20,20,1000,1000);
	c2p -> cd();
	h1 -> Draw();
	TPad *npad = new TPad("npad", "", 0.5051, 0.5032, 0.9068,0.8564);
	npad -> SetFillStyle(0);
	npad -> Draw();
	npad -> cd();
	h2 -> Draw("same");
	npad ->Draw();
  
	// timing plot for proton, pion and kaon in the same plot 

	
    
     
      
     
      
      
      TCanvas *c2 = new TCanvas("c2","",20,20,1000,1000);
       c2 -> Divide(2,2);
   
      c2 -> cd(1);
     
      TH2F *h5 = new TH2F("h5", "", 200, 0,5 , 200, -20,20);
      TH2F *h4p = new TH2F("h4p", "", 200, 0,5 , 200, -20,20);
      TH2F *h4pi = new TH2F("h4pi", "", 200, 0,5 , 200, -20,20);
      h5 -> GetXaxis()-> SetTitle("Momentums");
      h5 -> GetYaxis()-> SetTitle("#Deltat");
     
      tf -> Draw("(vtime - (sctkp - sclkp/(betakp*29.99))):kpp>> h5"); 
      h5 -> Draw("colz");
      c2->cd(2);
      tf -> Draw("(vtime - (sctkp - sclkp/(kpp/sqrt(0.938**2+kpp**2)*29.99))):kpp>> h4p"); 
      c2->cd(3);
      tf -> Draw("(vtime - (sctkp - sclkp/(kpp/sqrt(0.139**2+kpp**2)*29.99))):kpp>> h4pi"); 
      c2->cd(4);
      h5->Draw("COLZ");
      h4p->Draw("COLSAME");
      h4pi->Draw("COLSAME");
    
    
      
      /*
      //Add all the histograms in this root macro
      outFile->Append(g);
      
       outFile->Append(g1);
        outFile->Append(g2);
	 outFile->Append(g3);
	  outFile->Append(gc);
       outFile->Append(g1c);
        outFile->Append(g2c);
	 outFile->Append(g3c);
	  outFile->Append(h);
	   outFile->Append(h1);
	    outFile->Append(f1);
	     outFile->Append(f1);
	      outFile->Append(f2);
	       outFile->Append(f3);
	        outFile->Append(f4);
      
      outFile->Write();
      outFile -> Close();
      c -> Close();	
      cp -> Close();
      c2p -> Close();
      c1 -> Close();


	g -> Delete();
	g1 -> Delete();
	g2 -> Delete();
	g3 -> Delete();
	gc -> Delete();
	g1c -> Delete();
	g2c -> Delete();
	g3c -> Delete();
	h -> Delete();
	h1 -> Delete();
	f1 -> Delete();
	f2 -> Delete();
	f3 -> Delete();
	f4 -> Delete();



	

  	TCanvas *c1 = new TCanvas("c1","",20,20,1000,1000);
	c1 -> cd();
     
	         
	TH2F *f1 = new TH2F("f1", "", 200, 0,5, 200, 0,1.2);
	f1 -> GetXaxis()-> SetTitle("");
	f1 -> GetYaxis()-> SetTitle("");
	tf -> Draw("tofbetakp:kpp>> f1"); 
	 
      
	 
	TH2F *f3 = new TH2F("f3", "", 200, 0,5 , 200, 0,1.2);
	f3 -> GetXaxis()-> SetTitle("p [GeV/c]");
	f3 -> GetYaxis()-> SetTitle("#beta");
	tf -> Draw("tofbetap:pp>> f3");
	f3 -> GetYaxis() -> CenterTitle();
	f3 -> GetXaxis() -> CenterTitle();
	f3 -> Draw("col");
	f1 -> Draw("colsame");


	

         //    Top three after now are for delta t vs momentum plot for kp, p, pim seperately
      //                       without any particle selection
      TCanvas *c1 = new TCanvas("c1","",20,20,1000,1000);
       	        
      TH2F *f1 = new TH2F("f1", "", 200, 0,5 , 200, -20,20);
      f1 -> GetXaxis()-> SetTitle("Momentum(k^{+})");
      f1 -> GetYaxis()-> SetTitle("#Deltat");
      tf -> Draw("(sclkp/(tofbetakp*29.99)*(1 - sqrt((pow(0.938,2)+ pow(kpp,2))/(pow(tofmasskp,2) + pow(kpp,2))))):kpp>> f1"); 
      
       
      
      c1 -> cd();   
      TH2F *f2 = new TH2F("f2", "", 200, 0,5 , 200, -20,20);
      f2 -> GetXaxis()-> SetTitle("Momentum(p)");
      f2 -> GetYaxis()-> SetTitle("#Deltat");
      tf -> Draw("(sclkp/(tofbetakp*29.99)*(1 - sqrt((pow(0.4936,2)+ pow(kpp,2))/(pow(tofmasskp,2) + pow(kpp,2))))):kpp>> f2");
      f2 -> Draw("col");
      f1 -> Draw("colsame");

*/
      TCanvas *cPull = new TCanvas("cPull"," ", 20,20,1000,1000);
      cPull -> Divide(2,4);
      TH1F *l1 = new TH1F("l1","pull[0] ",100,-4,4);
      TH1F *l2 = new TH1F("l2","pull[1] ",100,-4,4);
      TH1F *l3 = new TH1F("l3","pull[2] ",100,-4,4);
      TH1F *l4 = new TH1F("l4","pull[3] ",100,-4,4);
      TH1F *l5 = new TH1F("l5","pull[4] ",100,-4,4);
      TH1F *l6 = new TH1F("l6","pull[5] ",100,-4,4);
      TH1F *l7 = new TH1F("l7"," pull[6]",100,-4,4);
      TH1F *l8 = new TH1F("l8","prob",100,0,1);

      cPull -> cd(2);
      tf -> Draw("pull[0]>>l1",mm&&mmKP&&nPlus&&nNeu);
      cPull -> cd(3);
      tf -> Draw("pull[1]>>l2",mm&&mmKP&&nPlus&&nNeu);
      cPull -> cd(4);
      tf -> Draw("pull[2]>>l3",mm&&mmKP&&nPlus&&nNeu);
      cPull -> cd(5);
      tf -> Draw("pull[3]>>l4",mm&&mmKP&&nPlus&&nNeu);
      cPull -> cd(6);
      tf -> Draw("pull[4]>>l5",mm&&mmKP&&nPlus&&nNeu);
      cPull -> cd(7);
      tf -> Draw("pull[5]>>l6",mm&&mmKP&&nPlus&&nNeu);
      cPull -> cd(8);
      tf -> Draw("pull[6]>>l7",mm&&mmKP&&nPlus&&nNeu);

      cPull -> cd(1);
      gPad -> SetLogy();
      tf -> Draw("prob>>l8",mm&&mmKP&&nPlus&&nNeu);
      l8 -> GetXaxis() -> SetTitle("prob");
	 cPull -> Print("pullDistributionAfterMK.pdf");
   

	 TCanvas *cPull = new TCanvas("cPull"," ", 20,20,1000,1000);
      cPull -> Divide(2,4);
      TH1F *l1 = new TH1F("l1","pull[0] ",100,-4,4);
      TH1F *l2 = new TH1F("l2","pull[1] ",100,-4,4);
      TH1F *l3 = new TH1F("l3","pull[2] ",100,-4,4);
      TH1F *l4 = new TH1F("l4","pull[3] ",100,-4,4);
      TH1F *l5 = new TH1F("l5","pull[4] ",100,-4,4);
      TH1F *l6 = new TH1F("l6","pull[5] ",100,-4,4);
      TH1F *l7 = new TH1F("l7"," pull[6]",100,-4,4);
      TH1F *l8 = new TH1F("l8","prob",100,0,1);

      cPull -> cd(2);
      tf -> Draw("pull[0]>>l1");
      cPull -> cd(3);
      tf -> Draw("pull[1]>>l2");
      cPull -> cd(4);
      tf -> Draw("pull[2]>>l3");
      cPull -> cd(5);
      tf -> Draw("pull[3]>>l4");
      cPull -> cd(6);
      tf -> Draw("pull[4]>>l5");
      cPull -> cd(7);
      tf -> Draw("pull[5]>>l6");
      cPull -> cd(8);
      tf -> Draw("pull[6]>>l7");

      cPull -> cd(1);
      gPad -> SetLogy();
      tf -> Draw("prob>>l8");
      l8 -> GetXaxis() -> SetTitle("prob");
	 // cPull -> Print("pullDistributionAfterMK.pdf");



*/
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	 // missing mass vs proton momentum
	 TCanvas *cprobMass = new TCanvas("cprobMass","",20,20,1000,1000);
	 cprobMass -> cd();
	 cprobMass -> Divide(2,3);


	 cprobMass -> cd(1);
	TH2F *mmVsmom = new TH2F("mmVsmom", "", 200, 0,5 , 200, 0.05,0.25);
	mmVsmom -> GetXaxis()-> SetTitle("p [GeV/c]");
	mmVsmom -> GetYaxis()-> SetTitle("mmkpp");
	tf -> Draw("sqrt(mm2pkp):pp>> mmVsmom");
	mmVsmom -> GetYaxis() -> CenterTitle();
	mmVsmom -> GetXaxis() -> CenterTitle();
	mmVsmom -> Draw("col");

	cprobMass -> cd(2);
	//gPad -> SetLogy();
	TH2F *mmVsProb = new TH2F("mmVsProb", "",500, -0.01,1.01, 500,0,0.3);
	mmVsProb -> GetYaxis()-> SetTitle("mmkpp");
	mmVsProb -> GetXaxis()-> SetTitle("prob");
	tf -> Draw("sqrt(mm2pkp):prob>> mmVsProb");
	mmVsProb -> GetYaxis() -> CenterTitle();
	mmVsProb -> GetXaxis() -> CenterTitle();
	//mmVsProb -> Draw("colz");

	cprobMass -> cd(3);
	TH2F *mm_kp_Vs_prob = new TH2F("mm_kp_Vs_prob", "", 200, -0.01,1.01 , 500, 1.08,1.15);
	mm_kp_Vs_prob -> GetXaxis()-> SetTitle("prob");
	mm_kp_Vs_prob -> GetYaxis()-> SetTitle("mmkp");
	tf -> Draw("sqrt(mm2kp):prob>> mm_kp_Vs_prob");
	mm_kp_Vs_prob -> GetYaxis() -> CenterTitle();
	mm_kp_Vs_prob -> GetXaxis() -> CenterTitle();
	mm_kp_Vs_prob -> Draw("col");

	cprobMass -> cd(4);
	//gPad -> SetLogy();
	TH2F *mm_ppim_vs_prob = new TH2F("mm_ppim_vs_prob", "",500, -0.01,1.01, 500,0.0,0.9);
	mm_ppim_vs_prob -> GetYaxis()-> SetTitle("mmppim");
	mm_ppim_vs_prob -> GetXaxis()-> SetTitle("prob");
	tf -> Draw("sqrt(mm2ppim):prob>> mm_ppim_vs_prob");
	mm_ppim_vs_prob -> GetYaxis() -> CenterTitle();
	mm_ppim_vs_prob -> GetXaxis() -> CenterTitle();
	mm_ppim_vs_prob -> Draw("col");

	cprobMass -> cd(5);
	//gPad -> SetLogy();
	TH2F *mm_kppim_vs_prob = new TH2F("mm_kppim_vs_prob", "",500, -0.01,1.01, 500,0.7,1.1);
	mm_kppim_vs_prob -> GetYaxis()-> SetTitle("mmppim");
	mm_kppim_vs_prob -> GetXaxis()-> SetTitle("prob");
	tf -> Draw("sqrt(mm2kppim):prob>> mm_kppim_vs_prob");
	mm_kppim_vs_prob -> GetYaxis() -> CenterTitle();
	mm_kppim_vs_prob -> GetXaxis() -> CenterTitle();
	mm_kppim_vs_prob -> Draw("col");




	// 2 dimensional plots 
/*

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	 // missing mass vs proton momentum
	 TCanvas *cLambda_mass = new TCanvas("cLambda_mass","",20,20,1000,1000);
	 cLambda_mass -> cd();
	 cLambda_mass -> Divide(2,2);


	cLambda_mass -> cd(1);
	TH2F *mass_Lambda_kpp = new TH2F("mass_Lambda_kpp", "", 200, 0,5.5 , 200, 1.0,1.2);
	mass_Lambda_kpp -> GetXaxis()-> SetTitle("kp_momentum [GeV/c]");
	mass_Lambda_kpp -> GetYaxis()-> SetTitle("mmkp");
	tf -> Draw("sqrt(mm2kp):kpp>> mass_Lambda_kpp");
	mass_Lambda_kpp -> GetYaxis() -> CenterTitle();
	mass_Lambda_kpp -> GetXaxis() -> CenterTitle();
	mass_Lambda_kpp -> Draw("col");

	cLambda_mass -> cd(2);
	//gPad -> SetLogy();
	TH2F *inv_mass_lambda = new TH2F("inv_mass_lambda", "",500, 0,5, 500,1.02,1.3);
	inv_mass_lambda -> GetYaxis()-> SetTitle("mmkp");
	inv_mass_lambda -> GetXaxis()-> SetTitle("kp_momentum");
	tf -> Draw("sqrt(invm2ppim):kpp >> inv_mass_lambda");
	inv_mass_lambda -> GetYaxis() -> CenterTitle();
	inv_mass_lambda -> GetXaxis() -> CenterTitle();
	inv_mass_lambda -> Draw("colz");


	cLambda_mass -> cd(3);
	//gPad -> SetLogy();
	TH2F *inv_mm_mLambda = new TH2F("inv_mm_mLambda", "",500, 1.0,1.2, 500,1.02,1.3);
	inv_mm_mLambda -> GetYaxis()-> SetTitle("sqrt(invm2ppim)");
	inv_mm_mLambda -> GetXaxis()-> SetTitle("sqrt(mm2kp)");
	tf -> Draw("sqrt(invm2ppim):sqrt(mm2kp) >> inv_mm_mLambda");
	inv_mm_mLambda -> GetYaxis() -> CenterTitle();
	inv_mm_mLambda -> GetXaxis() -> CenterTitle();
	inv_mm_mLambda -> Draw("colz");


	cLambda_mass -> cd(4);
	//gPad -> SetLogy();
	TH2F *energy_mLambda = new TH2F("energy_mLambda", "",500, 0.0,5.5, 500,1.02,1.2);
	energy_mLambda -> GetYaxis()-> SetTitle("sqrt(mm2kp)");
	energy_mLambda -> GetXaxis()-> SetTitle("beam energy");
	tf -> Draw("sqrt(mm2kp):ebeam >> energy_mLambda");
	energy_mLambda -> GetYaxis() -> CenterTitle();
	energy_mLambda -> GetXaxis() -> CenterTitle();
	energy_mLambda -> Draw("colz");


/*

      TCanvas *cFit = new TCanvas("cFit","", 20,20,1000,1000);
      cFit -> Divide(1,2);
      cFit -> cd(1);

      TH1F *hPiMass = new TH1F("hPiMass","", 1000, 0.0, 0.3);
      tf ->Draw("sqrt(mm2pkp)>>hPiMass",mm);
      
      TF1 *signal = new TF1("signal", , 0.02,0.22,3);
      TF1 *backgroundfit = new TF1("backgroundfit", backgroundExp,0.02,0.22,3);
      
      TF1 *parabola = new TF1("parabola","[0]+[1]*x+[2]*x**2",0.02,0.22);
      //format_line(&parabola,kBlue,2);

      TF1 *gaussian = new TF1("gaussian","[0]*TMath::Gaus(x,[1],[2])",0.02,0.22);
      //format_line(&gaussian,kRed,2);


      TF1 *fitter = new TF1("fitter",the_gausppar, 0.02, 0.22,6);  
	  
	  fitter->SetParName(0, "Constant");
	  fitter->SetParName(1, "Mean");
	  fitter->SetParName(2, "#sigma");
	  fitter->SetParName(3, "const term");
	  fitter->SetParName(4, "linear term");
	  fitter->SetParName(5, "quad term");
	  //fitter->SetParName(6, "cube term");
	  //fitter->SetParameters(1, 1., 1, 1, 1,1);
    
	  fitter->SetParameter(0,5000);
	  fitter->SetParameter(1,0.139);
	  fitter->SetParameter(2,0.000142);
	  fitter->SetParameter(3,-10.0);
	  fitter->SetParameter(4,0.5);
	  fitter->SetParameter(5,100.0);
     // fitter->SetParameter(6,1.0);
/*
      TF1 *fitter1 = new TF1("fitter1",GaussianPeak, 0.11, 0.16,3);  
	  fitter1->SetParName(0, "Norm");
	  fitter1->SetParName(1, "Factor");
	  fitter1->SetParName(2, "#sigma");
	  
	  //fitter->SetParName(6, "cube term");
	  fitter1->SetParameters(1, 1., 1);
    
	  fitter1->SetParameter(0,1000);
	  fitter1->SetParameter(1,0.134);
	  fitter1->SetParameter(2,0.001);
	  

	  hPiMass->Fit("fitter1","QR");
*/
	  

      
      
/*
	  hPiMass->Fit("fitter","R");

	  Double_t param[7];
	  fitter -> GetParameters(param);
	  gaussian -> SetParameters(&param[0]);
	  parabola -> SetParameters(&param[3]);


	  TH1F *hisSignal = new TH1F(*hPiMass);
	  hisSignal -> Sumw2();
	  hisSignal -> Add(parabola,-1);
	  TH1F *PureSignal = (TH1F*) hisSignal -> Clone();
	 // hisSignal -> Draw();
	//  hPiMass -> Draw("e");
	  hisSignal -> Draw("SAME");
	//  gaussian -> Draw("SAME");
	  parabola -> Draw("SAME");

	 
	  cFit -> cd(2);
	  PureSignal -> Draw();
*/



}







