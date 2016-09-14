#define KinFitTest_cxx
#include "KinFitTest.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

#include "../KinFit/KinLine.C"
#include "../KinFit/Kstream.C"

using namespace std;

void KinFitTest::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L KinFitTest.C
//      Root > KinFitTest t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   
   TFile* outputFile = new TFile("kpppim_kinfit.root","recreate");
   TTree* nTree = fChain->CloneTree(0);

   Long64_t nentries = fChain->GetEntriesFast();

   const Float_t PROTON_MASS = .938272;
   const Float_t PION_MASS = .139570;
   const Float_t KAON_MASS = .493667;
   string _exp = "g12";

   //Booleans for kinFit
   bool multi = true; 
   bool mc = false;

   //Number of particles for kinFit
   const int num_parts = 2;
   const int dim = 3*num_parts + 1;

   //Covariance matrix to be inserted to KinFit
   TMatrixD covTrack(dim,dim);
   TMatrixD covMatrix(dim,dim);
   vector<bool> set(num_parts,false); // tracks to be constrained are true.
   vector<string> n_particle(num_parts); // name of the particles
   vector<TVector3> vertex(num_parts); // vertices of the particles
   vector<TLorentzVector> outP(num_parts); // output of kinFitter
   vector<TLorentzVector> inP(num_parts); // input 4 vectors


   //Lab frame vectors
   TLorentzVector Beam, Beam_fit;
   TLorentzVector target, target_fit;
   TLorentzVector kp, kp_fit;
   TLorentzVector pim, pim_fit, miss_fit;
   TLorentzVector p, p_fit;

   //Post-Kinematic Fit
   Float_t s_fit;
   Float_t tkp;
   Float_t ebeam_fit;
   Float_t mpkp_fit;
   Float_t mppim_fit;
   Float_t mmkp_fit;
   Float_t mkppim_fit;
   Float_t mm2kp_fit;
   Float_t mm2kpp_fit;
   Float_t mm2p_fit;
   Float_t mm2ppim_fit;
   Float_t mm2kppim_fit;

   Float_t E_pfit;
   Float_t p_pfit;
   Float_t theta_pfit;
   Float_t phi_pfit;

   Float_t E_kpfit;
   Float_t p_kpfit;
   Float_t theta_kpfit;
   Float_t phi_kpfit;

   Float_t E_miss_fit;
   Float_t p_miss_fit;
   Float_t theta_miss_fit;
   Float_t phi_miss_fit;

   Float_t E_pimfit;
   Float_t p_pimfit;
   Float_t theta_pimfit;
   Float_t phi_pimfit;


   Float_t pull[dim];
   Float_t mpull[9];

    //Variables to add
   Float_t chi2;
   Float_t prob;

   //Branches
   TBranch* b_chi2 = nTree->Branch("chi2",&chi2,"chi2/F");
   TBranch* b_prob = nTree->Branch("prob",&prob,"prob/F");

   TBranch* b_s_fit = nTree->Branch("s",&s_fit,"s/F");
   TBranch* b_tkp = nTree->Branch("tkp",&tkp,"tkp/F");
   
   TBranch* b_ebeam_fit = nTree->Branch("ebeam_fit",&ebeam_fit,"ebeam_fit/F");
   TBranch* b_mpkp_fit = nTree->Branch("mpkp_fit",&mpkp_fit,"mpkp_fit/F");
   TBranch* b_mppim_fit = nTree->Branch("mppim_fit",&mppim_fit,"mppim_fit/F");
   TBranch* b_mkppim_fit = nTree->Branch("mkppim_fit",&mkppim_fit,"mkppim_fit/F");
   TBranch* b_mm2kp_fit = nTree->Branch("mm2kp_fit",&mm2kp_fit,"mm2kp_fit/F");
   TBranch* b_mm2kpp_fit = nTree->Branch("mm2kpp_fit",&mm2kpp_fit,"mm2kpp_fit/F");
   TBranch* b_mm2p_fit = nTree->Branch("mm2p_fit",&mm2p_fit,"mm2p_fit/F");
   TBranch* b_mm2ppim_fit = nTree->Branch("mm2ppim_fit",&mm2ppim_fit,"mm2ppim_fit/F");
   TBranch* b_mm2kppim_fit = nTree->Branch("mm2kppim_fit",&mm2kppim_fit,"mm2kppim_fit/F");

   TBranch* b_E_pfit = nTree->Branch("E_pfit",&E_pfit,"E_pfit/F");
   TBranch* b_p_pfit = nTree->Branch("p_pfit",&p_pfit,"p_pfit/F");
   TBranch* b_theta_pfit = nTree->Branch("theta_pfit",&theta_pfit,"theta_pfit/F");
   TBranch* b_phi_pfit = nTree->Branch("phi_pfit",&phi_pfit,"phi_pfit/F");

   TBranch* b_E_kpfit = nTree->Branch("E_kpfit",&E_kpfit,"E_kpfit/F");
   TBranch* b_p_kpfit = nTree->Branch("p_kpfit",&p_kpfit,"p_kpfit/F");
   TBranch* b_theta_kpfit = nTree->Branch("theta_kpfit",&theta_kpfit,"theta_kpfit/F");
   TBranch* b_phi_kpfit = nTree->Branch("phi_kpfit",&phi_kpfit,"phi_kpfit/F");

   TBranch* b_E_pimfit = nTree->Branch("E_pimfit",&E_pimfit,"E_pimfit/F");
   TBranch* b_p_pimfit = nTree->Branch("p_pimfit",&p_pimfit,"p_pimfit/F");
   TBranch* b_theta_pimfit = nTree->Branch("theta_pimfit",&theta_pimfit,"theta_pimfit/F");
   TBranch* b_phi_pimfit = nTree->Branch("phi_pimfit",&phi_pimfit,"phi_pimfit/F");


   TBranch* b_pull = nTree->Branch("pull",&pull,"pull[7]/F");
   TBranch* b_mpull = nTree->Branch("mpull",&mpull,"mpull[9]/F");



   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      Beam.SetXYZM(0.0,0.0,ebeam,0.0);
      target.SetXYZM(0.0,0.0,0.0,PROTON_MASS);
      kp.SetXYZM(kppx,kppy,kppz,KAON_MASS);
      p.SetXYZM(ppx,ppy,ppz,PROTON_MASS);

      target_fit = target;

      inP[0] = p;
      inP[1] = kp;

      //Setting Covariance matrix
      covTrack(0,0) = pow((0.001*5.715),2.0)/3.0;
      covTrack(1,1) = pow(p.P(),4.0) * c11_p; 
      covTrack(1,2) = covTrack(2,1) = -pow(p.P(),2.0) * c12_p;
      covTrack(1,3) = covTrack(3,1) = -pow(p.P(),2.0) * c13_p;
      covTrack(2,2) = c22_p;
      covTrack(2,3) = covTrack(3,2) = c23_p;
      covTrack(3,3) = c33_p;

      covTrack(4,4) = pow(kp.P(),4.0) * c11_kp;
      covTrack(4,5) = covTrack(5,4) = -pow(kp.P(),2.0) * c12_kp;
      covTrack(4,6) = covTrack(6,4) = -pow(kp.P(),2.0) * c13_kp;
      covTrack(5,5) = c22_kp;
      covTrack(5,6) = covTrack(6,5) = c23_kp;
      covTrack(6,6) = c33_kp;

      //Setting Vertex
      vertex[0].SetX(vx);
      vertex[0].SetY(vy);
      vertex[0].SetZ(vz);
      vertex[1].SetX(vx);
      vertex[1].SetY(vy);
      vertex[1].SetZ(vz);
      

	  //Setting Name of Particles
      n_particle[0] = "p";
      n_particle[1] = "k+";      


      //Setting Constraint bools
      set[0] = false;
      set[1] = false;

      bool include = false; // include missing particle
      covMatrix = CorrectCLAS_V(covTrack,n_particle,inP,vertex,multi,mc,_exp);
      Kstream testFit;
      testFit.StringNames(n_particle);
      testFit.FitInput(ebeam,inP,covMatrix,PROTON_MASS);
      testFit.Fit(0.13957);

      // testFit.Fit("pi-",set,include); // testFit.Fit("pi-");

      chi2 = testFit.Chi2();
      prob = testFit.Prob();

      ebeam_fit = testFit.FitPhotonEnergy();

      Beam_fit.SetXYZT(0.0,0.0,ebeam_fit,ebeam_fit);
      p_fit = testFit.FitP4(0);
      kp_fit = testFit.FitP4(1);

      miss_fit = Beam_fit + target_fit - p_fit - kp_fit;

      E_pfit = p_fit.T();
      p_pfit = p_fit.Vect().Mag();
      theta_pfit = p_fit.Vect().Theta();
      phi_pfit = p_fit.Vect().Phi();

      E_kpfit = kp_fit.T();
      p_kpfit = kp_fit.Vect().Mag();
      theta_kpfit = kp_fit.Vect().Theta();
      phi_kpfit = kp_fit.Vect().Phi();

      E_miss_fit = miss_fit.T();
      p_miss_fit = miss_fit.Vect().Mag();
      theta_miss_fit = miss_fit.Vect().Theta();
      phi_miss_fit = miss_fit.Vect().Phi();

      pim_fit.SetXYZM(miss_fit.X(),miss_fit.Y(),miss_fit.Z(),PION_MASS);
      E_pimfit = pim_fit.T();
      p_pimfit = pim_fit.Vect().Mag();
      theta_pimfit = pim_fit.Vect().Theta();
      phi_pimfit = pim_fit.Vect().Phi();



      s_fit = (Beam_fit + target).M2();
      tkp = (Beam_fit - kp_fit).M2();

      mpkp_fit = (p_fit + kp_fit).M();
      mppim_fit = (p_fit + pim_fit).M();
      mkppim_fit = (kp_fit + pim_fit).M();

      // missing masses
      mm2kp_fit = (Beam_fit + target - kp_fit).M2();
      mm2kpp_fit = (Beam_fit + target - p_fit - kp_fit).M2();
      mm2p_fit = (Beam_fit + target - p_fit).M2();
      mm2ppim_fit = (Beam_fit + target - p_fit - pim_fit).M2();
      mm2kppim_fit = (Beam_fit + target - pim_fit - kp_fit).M2();

      if (testFit.Prob() > 0.05){
      for(int i = 0; i < dim; i++)
         pull[i] = testFit.GetPull(i);
      
      mpull[0] = ebeam_fit - ebeam;
      mpull[1] = p_fit.T() - p.T();
      mpull[2] = p_fit.X() - p.X();
      mpull[3] = p_fit.Y() - p.Y();
      mpull[4] = p_fit.Z() - p.Z();
      mpull[5] = kp_fit.T() - kp.T();
      mpull[6] = kp_fit.X() - kp.X();
      mpull[7] = kp_fit.Y() - kp.Y();
      mpull[8] = kp_fit.Z() - kp.Z();
      
      }
		nTree->Fill();
      
      	if(jentry%10000==0)
         cout << (float)jentry / (float) nentries * 100.0 << "%" << "     \r" << flush;
      
   }
   cout << 100.0 << "%\n     \r" << flush;
   //nTree->Print();
   nTree->Write();
   outputFile->Close();
}
