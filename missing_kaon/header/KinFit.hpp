#define KinFit_cxx
#include "KinFit.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "../KinFit/KinLine.C"
#include "../KinFit/Kstream.C"

using namespace std;

void KinFit::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L KinFit.C
//      Root > KinFit t
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

   TFile *output = new TFile("kpppim_kinfit.root", "recreate");
   TTree *nTree = fChain -> CloneTree(0);
   
   Long64_t nentries = fChain->GetEntriesFast();

   const Float_t PROTON_MASS = .938272;
   const Float_t PION_MASS = .13957;
   const Float_t KAON_MASS = .493667;
   string _exp = "g12";
   
   //Booleans for kinFit
   bool multi = false;// true; 
   bool mc = false;

   //number of particle and dimension 
   const int num_parts = 2;
   const int dim= 3*num_parts +1;
   //Covariance matrix to be inserted to KinFit
   TMatrixD covTrack(dim,dim);  // matrix will be use
   TMatrixD covMatrix(dim,dim);  // matrix from CLAS TBER
   vector<bool> set(num_parts, false);   // track to be constrained are true.
   vector<string> n_particle(num_parts);   // name of the particles
   vector<TVector3> vertex(num_parts); // vertices of the particles
   vector<TLorentzVector> outP(num_parts); // output of kinFitter
   vector<TLorentzVector> inP(num_parts); // input 4 vectors
   
   // Lab frame vectors 

   TLorentzVector Beam, Beam_fit;
   TLorentzVector target, target_fit;
  
   TLorentzVector p, p_fit;
   TLorentzVector pim, pim_fit;



   //Variables to add
   Float_t chi2;
   Float_t prob;
  

   
   Float_t ebeam_fit;
   Float_t E_pfit;
   Float_t p_pfit;
   Float_t theta_pfit;
   Float_t phi_pfit;
   Float_t E_pimfit;
   Float_t p_pimfit;
   Float_t theta_pimfit;
   Float_t phi_pimfit;


   Float_t pull[dim];
   Float_t mpull[7];

   // defind energy for three particles
   Float_t ebeam_p;
   Float_t ebeam_pim;



   TBranch* b_chi2 = nTree->Branch("chi2",&chi2,"chi2/F");
   TBranch* b_prob = nTree->Branch("prob",&prob,"prob/F");
 

  
  
   TBranch* b_E_pfit = nTree->Branch("E_pfit",&E_pfit,"E_pfit/F");
   TBranch* b_p_pfit = nTree->Branch("p_pfit",&p_pfit,"p_pfit/F");
   TBranch* b_theta_pfit = nTree->Branch("theta_pfit",&theta_pfit,"theta_pfit/F");
   TBranch* b_phi_pfit = nTree->Branch("phi_pfit",&phi_pfit,"phi_pfit/F");
   TBranch* b_E_pimfit = nTree->Branch("E_pimfit",&E_pimfit,"E_pimfit/F");
   TBranch* b_p_pimfit = nTree->Branch("p_pimfit",&p_pimfit,"p_pimfit/F");
   TBranch* b_theta_pimfit = nTree->Branch("theta_pimfit",&theta_pimfit,"theta_pimfit/F");
   TBranch* b_phi_pimfit = nTree->Branch("phi_pimfit",&phi_pimfit,"phi_pimfit/F");


   TBranch* b_pull = nTree->Branch("pull",&pull,"pull[7]/F");
   TBranch* b_mpull = nTree->Branch("mpull",&mpull,"mpull[7]/F");
 
   // loop start here

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;
   
     // calculating the energy 
      
     ebeam_p = sqrt(pow(pp,2) + pow(PROTON_MASS,2));
     ebeam_pim = sqrt(pow(ppim,2) + pow(PION_MASS,2));


     //Setting 4vectors

     Beam.SetXYZT(0.0,0.0,beam,beam);
     target.SetXYZM(0.0,0.0,0.0,PROTON_MASS);
     target_fit = target;
      
     
      
     
     p.SetX(pxp);
     p.SetY(pyp);
     p.SetZ(pzp);
     p.SetT(ebeam_p);

     pim.SetX(pxpim);
     pim.SetY(pypim);
     pim.SetZ(pzpim);
     pim.SetT(ebeam_pim);

     
     inP[0] = p;
     inP[1] = pim;



     //Setting Covariance matrix
     



     covTrack(0,0) = pow((0.001*5.715),2.0)/3.0;
     covTrack(1,1) = pow(pp,4.0) * c11_p; 
     covTrack(1,2) = covTrack(2,1) = -pow(pp,2.0) * c12_p;
     covTrack(1,3) = covTrack(3,1) = -pow(pp,2.0) * c13_p;
     covTrack(2,2) = c22_p;
     covTrack(2,3) = covTrack(3,2) = c23_p;
     covTrack(3,3) = c33_p;


     covTrack(4,4) = pow(ppim,4.0) * c11_pim;
     covTrack(4,5) = covTrack(5,4) = -pow(ppim,2.0) * c12_pim;
     covTrack(4,6) = covTrack(6,4) = -pow(ppim,2.0) * c13_pim;
     covTrack(5,5) = c22_pim;
     covTrack(5,6) = covTrack(6,5) = c23_pim;
     covTrack(6,6) = c33_pim;



     //Setting Vertex
     vertex[0].SetX(vx);
     vertex[0].SetY(vy);
     vertex[0].SetZ(vz);
     vertex[1].SetX(vx);
     vertex[1].SetY(vy);    
     vertex[1].SetZ(vz);

     //Setting Name of particles
     n_particle[0] = "p";
     n_particle[1] = "pi-";
    

     //Setting Constraint bools
     set[0] = false;
     set[1] = false;
    
     Kstream testFit; 
     covMatrix = CorrectCLAS_V(covTrack,n_particle,inP,vertex,multi,mc,_exp);
     testFit.StringNames(n_particle);
     testFit.FitInput(beam,inP,covMatrix,PROTON_MASS);
     bool include = true;
     testFit.Fit(0.493667);

     chi2 = testFit.Chi2();
     prob = testFit.Prob();

     ebeam_fit = testFit.FitPhotonEnergy();
     Beam_fit.SetXYZT(0.0,0.0,ebeam_fit,ebeam_fit);
     p_fit = testFit.FitP4(0);
     pim_fit = testFit.FitP4(1);
    


     E_pfit = p_fit.T();
     p_pfit = p_fit.Vect().Mag();
     theta_pfit = p_fit.Vect().Theta();
     phi_pfit = p_fit.Vect().Phi();


     E_pimfit = pim_fit.T();
     p_pimfit = pim_fit.Vect().Mag();
     theta_pimfit = pim_fit.Vect().Theta();
     phi_pimfit = pim_fit.Vect().Phi();

      

     for(int i = 0; i < dim; i++)
       {pull[i] = testFit.GetPull(i);}

     mpull[0] = (ebeam_fit- beam)/(float)beam;
     mpull[1] = (p_pfit*sin(theta_pfit)*cos(phi_pfit) - pxp)/(float)pxp;
     mpull[2] = (p_pfit*sin(theta_pfit)*sin(phi_pfit) - pyp)/(float)pyp;
     mpull[3] = (p_pfit*cos(theta_pfit) - pzp)/(float)pzp;
     mpull[4] = (p_pimfit*sin(theta_pimfit)*cos(phi_pimfit) - pxpim)/(float)pxpim;
     mpull[5] = (p_pimfit*sin(theta_pimfit)*sin(phi_pimfit) - pypim)/(float)pypim;
     mpull[6] = (p_pimfit*cos(theta_pimfit) - pzpim)/(float)pzpim;



     nTree->Fill();
     if(jentry%10000==0)
       cout << (float)jentry / (float) nentries * 100.0 << "%" << "\r" << flush;
      

   }
   cout << 100.0 << "%\n\r" << flush;
   //nTree->Print();
   nTree->Write();
   output ->Close();

}

     
