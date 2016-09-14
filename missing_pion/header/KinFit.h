//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep  1 12:56:45 2016 by ROOT version 5.34/05
// from TTree output/data
// found on file: kpppim.root
//////////////////////////////////////////////////////////

#ifndef KinFit_h
#define KinFit_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class KinFit {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         beam_no_ecor;
   Float_t         mmppim_no_ecor;
   Float_t         mmkpppim_no_ecor;
   Float_t         mmppim_no_pcor;
   Float_t         mmkpppim_no_pcor;
   Float_t         mmppim_no_eloss;
   Float_t         mmkpppim_no_eloss;
   Float_t         mmppim_no_corrections;
   Float_t         mmkpppim_no_corrections;
   Float_t         mmppim;
   Float_t         mmkpppim;
   Float_t         run;
   Float_t         evt;
   Float_t         scv;
   Float_t         stv;
   Float_t         vtime;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         nkp;
   Float_t         npim;
   Float_t         np;
   Float_t         tofbetakp;
   Float_t         tofbetap;
   Float_t         tofbetapim;
   Float_t         sctkp;
   Float_t         sctp;
   Float_t         sctpim;
   Float_t         sclkp;
   Float_t         sclp;
   Float_t         sclpim;
   Float_t         tofmasskp;
   Float_t         tofmassp;
   Float_t         tofmasspim;
   Float_t         stvkp;
   Float_t         stvp;
   Float_t         stvpim;
   Float_t         goodSTkp;
   Float_t         goodSTp;
   Float_t         goodSTpim;
   Float_t         stid_kp;
   Float_t         stid_p;
   Float_t         stid_pim;
   Float_t         scidkp;
   Float_t         sectkp;
   Float_t         tofkp;
   Float_t         tofexpkp;
   Float_t         qkp;
   Float_t         scidp;
   Float_t         sectp;
   Float_t         tofp;
   Float_t         tofexpp;
   Float_t         qp;
   Float_t         scidpim;
   Float_t         sectpim;
   Float_t         tofpim;
   Float_t         tofexppim;
   Float_t         qpim;
   Float_t         dkppim;
   Float_t         dppim;
   Float_t         dkpp;
   Float_t         beam;
   Float_t         w;
   Float_t         tlam;
   Float_t         uppim;
   Float_t         tkp;
   Float_t         mm2kpp;
   Float_t         mm2kp;
   Float_t         mm2p;
   Float_t         mm2kppim;
   Float_t         mm2ppim;
   Float_t         mm2kpppim;
   Float_t         mkpgamma;
   Float_t         mppimgamma;
   Float_t         mlamgamma;
   Float_t         mkpp;
   Float_t         mkppim;
   Float_t         mppim;
   Float_t         mmkp;
   Float_t         mmp;
   Float_t         pkp;
   Float_t         pxkp;
   Float_t         pykp;
   Float_t         pzkp;
   Float_t         pp;
   Float_t         pxp;
   Float_t         pyp;
   Float_t         pzp;
   Float_t         ppim;
   Float_t         pxpim;
   Float_t         pypim;
   Float_t         pzpim;
   Float_t         missx;
   Float_t         missy;
   Float_t         missz;
   Float_t         missp;
   Float_t         misse;
   Float_t         nphoton_nb0;
   Float_t         nphoton_nb1;
   Float_t         nphoton_nb2;
   Float_t         betakp;
   Float_t         betap;
   Float_t         betapim;
   Float_t         tkpp;
   Float_t         c11_kp;
   Float_t         c12_kp;
   Float_t         c13_kp;
   Float_t         c14_kp;
   Float_t         c15_kp;
   Float_t         c22_kp;
   Float_t         c23_kp;
   Float_t         c24_kp;
   Float_t         c25_kp;
   Float_t         c33_kp;
   Float_t         c34_kp;
   Float_t         c35_kp;
   Float_t         c44_kp;
   Float_t         c45_kp;
   Float_t         c55_kp;
   Float_t         c11_p;
   Float_t         c12_p;
   Float_t         c13_p;
   Float_t         c14_p;
   Float_t         c15_p;
   Float_t         c22_p;
   Float_t         c23_p;
   Float_t         c24_p;
   Float_t         c25_p;
   Float_t         c33_p;
   Float_t         c34_p;
   Float_t         c35_p;
   Float_t         c44_p;
   Float_t         c45_p;
   Float_t         c55_p;
   Float_t         c11_pim;
   Float_t         c12_pim;
   Float_t         c13_pim;
   Float_t         c14_pim;
   Float_t         c15_pim;
   Float_t         c22_pim;
   Float_t         c23_pim;
   Float_t         c24_pim;
   Float_t         c25_pim;
   Float_t         c33_pim;
   Float_t         c34_pim;
   Float_t         c35_pim;
   Float_t         c44_pim;
   Float_t         c45_pim;
   Float_t         c55_pim;
   Float_t         COS_pi_x;
   Float_t         COS_pi_y;
   Float_t         COS_pi_z;
   Float_t         COS_p_x;
   Float_t         COS_p_y;
   Float_t         COS_p_z;
   Float_t         kplabtheta;
   Float_t         kplabphi;
   Float_t         plabtheta;
   Float_t         plabphi;
   Float_t         pilabtheta;
   Float_t         pilabphi;
   Float_t         Xlabtheta;
   Float_t         Xlabphi;
   Float_t         hel;
   Float_t         COSth_pi_x;
   Float_t         COSth_pi_y;
   Float_t         COSth_pi_z;
   Float_t         COSth_p_x;
   Float_t         COSth_p_y;
   Float_t         COSth_p_z;
   Float_t         costheta;
   Float_t         COSth_k_z_cm;
   Float_t         beamPol;

   // List of branches
   TBranch        *b_beam_no_ecor;   //!
   TBranch        *b_mmppim_no_ecor;   //!
   TBranch        *b_mmkpppim_no_ecor;   //!
   TBranch        *b_mmppim_no_pcor;   //!
   TBranch        *b_mmkpppim_no_pcor;   //!
   TBranch        *b_mmppim_no_eloss;   //!
   TBranch        *b_mmkpppim_no_eloss;   //!
   TBranch        *b_mmppim_no_corrections;   //!
   TBranch        *b_mmkpppim_no_corrections;   //!
   TBranch        *b_mmppim;   //!
   TBranch        *b_mmkpppim;   //!
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_scv;   //!
   TBranch        *b_stv;   //!
   TBranch        *b_vtime;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_nkp;   //!
   TBranch        *b_npim;   //!
   TBranch        *b_np;   //!
   TBranch        *b_tofbetakp;   //!
   TBranch        *b_tofbetap;   //!
   TBranch        *b_tofbetapim;   //!
   TBranch        *b_sctkp;   //!
   TBranch        *b_sctp;   //!
   TBranch        *b_sctpim;   //!
   TBranch        *b_sclkp;   //!
   TBranch        *b_sclp;   //!
   TBranch        *b_sclpim;   //!
   TBranch        *b_tofmasskp;   //!
   TBranch        *b_tofmassp;   //!
   TBranch        *b_tofmasspim;   //!
   TBranch        *b_stvkp;   //!
   TBranch        *b_stvp;   //!
   TBranch        *b_stvpim;   //!
   TBranch        *b_goodSTkp;   //!
   TBranch        *b_goodSTp;   //!
   TBranch        *b_goodSTpim;   //!
   TBranch        *b_stid_kp;   //!
   TBranch        *b_stid_p;   //!
   TBranch        *b_stid_pim;   //!
   TBranch        *b_scidkp;   //!
   TBranch        *b_sectkp;   //!
   TBranch        *b_tofkp;   //!
   TBranch        *b_tofexpkp;   //!
   TBranch        *b_qkp;   //!
   TBranch        *b_scidp;   //!
   TBranch        *b_sectp;   //!
   TBranch        *b_tofp;   //!
   TBranch        *b_tofexpp;   //!
   TBranch        *b_qp;   //!
   TBranch        *b_scidpim;   //!
   TBranch        *b_sectpim;   //!
   TBranch        *b_tofpim;   //!
   TBranch        *b_tofexppim;   //!
   TBranch        *b_qpim;   //!
   TBranch        *b_dkppim;   //!
   TBranch        *b_dppim;   //!
   TBranch        *b_dkpp;   //!
   TBranch        *b_beam;   //!
   TBranch        *b_w;   //!
   TBranch        *b_tlam;   //!
   TBranch        *b_uppim;   //!
   TBranch        *b_tkp;   //!
   TBranch        *b_mm2kpp;   //!
   TBranch        *b_mm2kp;   //!
   TBranch        *b_mm2p;   //!
   TBranch        *b_mm2kppim;   //!
   TBranch        *b_mm2ppim;   //!
   TBranch        *b_mm2kpppim;   //!
   TBranch        *b_mkpgamma;   //!
   TBranch        *b_mppimgamma;   //!
   TBranch        *b_mlamgamma;   //!
   TBranch        *b_mkpp;   //!
   TBranch        *b_mkppim;   //!
   TBranch        *b_mppim;   //!
   TBranch        *b_mmkp;   //!
   TBranch        *b_mmp;   //!
   TBranch        *b_pkp;   //!
   TBranch        *b_pxkp;   //!
   TBranch        *b_pykp;   //!
   TBranch        *b_pzkp;   //!
   TBranch        *b_pp;   //!
   TBranch        *b_pxp;   //!
   TBranch        *b_pyp;   //!
   TBranch        *b_pzp;   //!
   TBranch        *b_ppim;   //!
   TBranch        *b_pxpim;   //!
   TBranch        *b_pypim;   //!
   TBranch        *b_pzpim;   //!
   TBranch        *b_missx;   //!
   TBranch        *b_missy;   //!
   TBranch        *b_missz;   //!
   TBranch        *b_missp;   //!
   TBranch        *b_misse;   //!
   TBranch        *b_nphoton_nb0;   //!
   TBranch        *b_nphoton_nb1;   //!
   TBranch        *b_nphoton_nb2;   //!
   TBranch        *b_betakp;   //!
   TBranch        *b_betap;   //!
   TBranch        *b_betapim;   //!
   TBranch        *b_tkpp;   //!
   TBranch        *b_c11_kp;   //!
   TBranch        *b_c12_kp;   //!
   TBranch        *b_c13_kp;   //!
   TBranch        *b_c14_kp;   //!
   TBranch        *b_c15_kp;   //!
   TBranch        *b_c22_kp;   //!
   TBranch        *b_c23_kp;   //!
   TBranch        *b_c24_kp;   //!
   TBranch        *b_c25_kp;   //!
   TBranch        *b_c33_kp;   //!
   TBranch        *b_c34_kp;   //!
   TBranch        *b_c35_kp;   //!
   TBranch        *b_c44_kp;   //!
   TBranch        *b_c45_kp;   //!
   TBranch        *b_c55_kp;   //!
   TBranch        *b_c11_p;   //!
   TBranch        *b_c12_p;   //!
   TBranch        *b_c13_p;   //!
   TBranch        *b_c14_p;   //!
   TBranch        *b_c15_p;   //!
   TBranch        *b_c22_p;   //!
   TBranch        *b_c23_p;   //!
   TBranch        *b_c24_p;   //!
   TBranch        *b_c25_p;   //!
   TBranch        *b_c33_p;   //!
   TBranch        *b_c34_p;   //!
   TBranch        *b_c35_p;   //!
   TBranch        *b_c44_p;   //!
   TBranch        *b_c45_p;   //!
   TBranch        *b_c55_p;   //!
   TBranch        *b_c11_pim;   //!
   TBranch        *b_c12_pim;   //!
   TBranch        *b_c13_pim;   //!
   TBranch        *b_c14_pim;   //!
   TBranch        *b_c15_pim;   //!
   TBranch        *b_c22_pim;   //!
   TBranch        *b_c23_pim;   //!
   TBranch        *b_c24_pim;   //!
   TBranch        *b_c25_pim;   //!
   TBranch        *b_c33_pim;   //!
   TBranch        *b_c34_pim;   //!
   TBranch        *b_c35_pim;   //!
   TBranch        *b_c44_pim;   //!
   TBranch        *b_c45_pim;   //!
   TBranch        *b_c55_pim;   //!
   TBranch        *b_COS_pi_x;   //!
   TBranch        *b_COS_pi_y;   //!
   TBranch        *b_COS_pi_z;   //!
   TBranch        *b_COS_p_x;   //!
   TBranch        *b_COS_p_y;   //!
   TBranch        *b_COS_p_z;   //!
   TBranch        *b_kplabtheta;   //!
   TBranch        *b_kplabphi;   //!
   TBranch        *b_plabtheta;   //!
   TBranch        *b_plabphi;   //!
   TBranch        *b_pilabtheta;   //!
   TBranch        *b_pilabphi;   //!
   TBranch        *b_Xlabtheta;   //!
   TBranch        *b_Xlabphi;   //!
   TBranch        *b_hel;   //!
   TBranch        *b_COSth_pi_x;   //!
   TBranch        *b_COSth_pi_y;   //!
   TBranch        *b_COSth_pi_z;   //!
   TBranch        *b_COSth_p_x;   //!
   TBranch        *b_COSth_p_y;   //!
   TBranch        *b_COSth_p_z;   //!
   TBranch        *b_costheta;   //!
   TBranch        *b_COSth_k_z_cm;   //!
   TBranch        *b_beamPol;   //!

   KinFit(TTree *tree=0);
   virtual ~KinFit();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef KinFit_cxx
KinFit::KinFit(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/kpppim.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../data/kpppim.root");
      }
      f->GetObject("output",tree);

   }
   Init(tree);
}

KinFit::~KinFit()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t KinFit::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t KinFit::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void KinFit::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("beam_no_ecor", &beam_no_ecor, &b_beam_no_ecor);
   fChain->SetBranchAddress("mmppim_no_ecor", &mmppim_no_ecor, &b_mmppim_no_ecor);
   fChain->SetBranchAddress("mmkpppim_no_ecor", &mmkpppim_no_ecor, &b_mmkpppim_no_ecor);
   fChain->SetBranchAddress("mmppim_no_pcor", &mmppim_no_pcor, &b_mmppim_no_pcor);
   fChain->SetBranchAddress("mmkpppim_no_pcor", &mmkpppim_no_pcor, &b_mmkpppim_no_pcor);
   fChain->SetBranchAddress("mmppim_no_eloss", &mmppim_no_eloss, &b_mmppim_no_eloss);
   fChain->SetBranchAddress("mmkpppim_no_eloss", &mmkpppim_no_eloss, &b_mmkpppim_no_eloss);
   fChain->SetBranchAddress("mmppim_no_corrections", &mmppim_no_corrections, &b_mmppim_no_corrections);
   fChain->SetBranchAddress("mmkpppim_no_corrections", &mmkpppim_no_corrections, &b_mmkpppim_no_corrections);
   fChain->SetBranchAddress("mmppim", &mmppim, &b_mmppim);
   fChain->SetBranchAddress("mmkpppim", &mmkpppim, &b_mmkpppim);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("scv", &scv, &b_scv);
   fChain->SetBranchAddress("stv", &stv, &b_stv);
   fChain->SetBranchAddress("vtime", &vtime, &b_vtime);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("nkp", &nkp, &b_nkp);
   fChain->SetBranchAddress("npim", &npim, &b_npim);
   fChain->SetBranchAddress("np", &np, &b_np);
   fChain->SetBranchAddress("tofbetakp", &tofbetakp, &b_tofbetakp);
   fChain->SetBranchAddress("tofbetap", &tofbetap, &b_tofbetap);
   fChain->SetBranchAddress("tofbetapim", &tofbetapim, &b_tofbetapim);
   fChain->SetBranchAddress("sctkp", &sctkp, &b_sctkp);
   fChain->SetBranchAddress("sctp", &sctp, &b_sctp);
   fChain->SetBranchAddress("sctpim", &sctpim, &b_sctpim);
   fChain->SetBranchAddress("sclkp", &sclkp, &b_sclkp);
   fChain->SetBranchAddress("sclp", &sclp, &b_sclp);
   fChain->SetBranchAddress("sclpim", &sclpim, &b_sclpim);
   fChain->SetBranchAddress("tofmasskp", &tofmasskp, &b_tofmasskp);
   fChain->SetBranchAddress("tofmassp", &tofmassp, &b_tofmassp);
   fChain->SetBranchAddress("tofmasspim", &tofmasspim, &b_tofmasspim);
   fChain->SetBranchAddress("stvkp", &stvkp, &b_stvkp);
   fChain->SetBranchAddress("stvp", &stvp, &b_stvp);
   fChain->SetBranchAddress("stvpim", &stvpim, &b_stvpim);
   fChain->SetBranchAddress("goodSTkp", &goodSTkp, &b_goodSTkp);
   fChain->SetBranchAddress("goodSTp", &goodSTp, &b_goodSTp);
   fChain->SetBranchAddress("goodSTpim", &goodSTpim, &b_goodSTpim);
   fChain->SetBranchAddress("stid_kp", &stid_kp, &b_stid_kp);
   fChain->SetBranchAddress("stid_p", &stid_p, &b_stid_p);
   fChain->SetBranchAddress("stid_pim", &stid_pim, &b_stid_pim);
   fChain->SetBranchAddress("scidkp", &scidkp, &b_scidkp);
   fChain->SetBranchAddress("sectkp", &sectkp, &b_sectkp);
   fChain->SetBranchAddress("tofkp", &tofkp, &b_tofkp);
   fChain->SetBranchAddress("tofexpkp", &tofexpkp, &b_tofexpkp);
   fChain->SetBranchAddress("qkp", &qkp, &b_qkp);
   fChain->SetBranchAddress("scidp", &scidp, &b_scidp);
   fChain->SetBranchAddress("sectp", &sectp, &b_sectp);
   fChain->SetBranchAddress("tofp", &tofp, &b_tofp);
   fChain->SetBranchAddress("tofexpp", &tofexpp, &b_tofexpp);
   fChain->SetBranchAddress("qp", &qp, &b_qp);
   fChain->SetBranchAddress("scidpim", &scidpim, &b_scidpim);
   fChain->SetBranchAddress("sectpim", &sectpim, &b_sectpim);
   fChain->SetBranchAddress("tofpim", &tofpim, &b_tofpim);
   fChain->SetBranchAddress("tofexppim", &tofexppim, &b_tofexppim);
   fChain->SetBranchAddress("qpim", &qpim, &b_qpim);
   fChain->SetBranchAddress("dkppim", &dkppim, &b_dkppim);
   fChain->SetBranchAddress("dppim", &dppim, &b_dppim);
   fChain->SetBranchAddress("dkpp", &dkpp, &b_dkpp);
   fChain->SetBranchAddress("beam", &beam, &b_beam);
   fChain->SetBranchAddress("w", &w, &b_w);
   fChain->SetBranchAddress("tlam", &tlam, &b_tlam);
   fChain->SetBranchAddress("uppim", &uppim, &b_uppim);
   fChain->SetBranchAddress("tkp", &tkp, &b_tkp);
   fChain->SetBranchAddress("mm2kpp", &mm2kpp, &b_mm2kpp);
   fChain->SetBranchAddress("mm2kp", &mm2kp, &b_mm2kp);
   fChain->SetBranchAddress("mm2p", &mm2p, &b_mm2p);
   fChain->SetBranchAddress("mm2kppim", &mm2kppim, &b_mm2kppim);
   fChain->SetBranchAddress("mm2ppim", &mm2ppim, &b_mm2ppim);
   fChain->SetBranchAddress("mm2kpppim", &mm2kpppim, &b_mm2kpppim);
   fChain->SetBranchAddress("mkpgamma", &mkpgamma, &b_mkpgamma);
   fChain->SetBranchAddress("mppimgamma", &mppimgamma, &b_mppimgamma);
   fChain->SetBranchAddress("mlamgamma", &mlamgamma, &b_mlamgamma);
   fChain->SetBranchAddress("mkpp", &mkpp, &b_mkpp);
   fChain->SetBranchAddress("mkppim", &mkppim, &b_mkppim);
   fChain->SetBranchAddress("mppim", &mppim, &b_mppim);
   fChain->SetBranchAddress("mmkp", &mmkp, &b_mmkp);
   fChain->SetBranchAddress("mmp", &mmp, &b_mmp);
   fChain->SetBranchAddress("pkp", &pkp, &b_pkp);
   fChain->SetBranchAddress("pxkp", &pxkp, &b_pxkp);
   fChain->SetBranchAddress("pykp", &pykp, &b_pykp);
   fChain->SetBranchAddress("pzkp", &pzkp, &b_pzkp);
   fChain->SetBranchAddress("pp", &pp, &b_pp);
   fChain->SetBranchAddress("pxp", &pxp, &b_pxp);
   fChain->SetBranchAddress("pyp", &pyp, &b_pyp);
   fChain->SetBranchAddress("pzp", &pzp, &b_pzp);
   fChain->SetBranchAddress("ppim", &ppim, &b_ppim);
   fChain->SetBranchAddress("pxpim", &pxpim, &b_pxpim);
   fChain->SetBranchAddress("pypim", &pypim, &b_pypim);
   fChain->SetBranchAddress("pzpim", &pzpim, &b_pzpim);
   fChain->SetBranchAddress("missx", &missx, &b_missx);
   fChain->SetBranchAddress("missy", &missy, &b_missy);
   fChain->SetBranchAddress("missz", &missz, &b_missz);
   fChain->SetBranchAddress("missp", &missp, &b_missp);
   fChain->SetBranchAddress("misse", &misse, &b_misse);
   fChain->SetBranchAddress("nphoton_nb0", &nphoton_nb0, &b_nphoton_nb0);
   fChain->SetBranchAddress("nphoton_nb1", &nphoton_nb1, &b_nphoton_nb1);
   fChain->SetBranchAddress("nphoton_nb2", &nphoton_nb2, &b_nphoton_nb2);
   fChain->SetBranchAddress("betakp", &betakp, &b_betakp);
   fChain->SetBranchAddress("betap", &betap, &b_betap);
   fChain->SetBranchAddress("betapim", &betapim, &b_betapim);
   fChain->SetBranchAddress("tkpp", &tkpp, &b_tkpp);
   fChain->SetBranchAddress("c11_kp", &c11_kp, &b_c11_kp);
   fChain->SetBranchAddress("c12_kp", &c12_kp, &b_c12_kp);
   fChain->SetBranchAddress("c13_kp", &c13_kp, &b_c13_kp);
   fChain->SetBranchAddress("c14_kp", &c14_kp, &b_c14_kp);
   fChain->SetBranchAddress("c15_kp", &c15_kp, &b_c15_kp);
   fChain->SetBranchAddress("c22_kp", &c22_kp, &b_c22_kp);
   fChain->SetBranchAddress("c23_kp", &c23_kp, &b_c23_kp);
   fChain->SetBranchAddress("c24_kp", &c24_kp, &b_c24_kp);
   fChain->SetBranchAddress("c25_kp", &c25_kp, &b_c25_kp);
   fChain->SetBranchAddress("c33_kp", &c33_kp, &b_c33_kp);
   fChain->SetBranchAddress("c34_kp", &c34_kp, &b_c34_kp);
   fChain->SetBranchAddress("c35_kp", &c35_kp, &b_c35_kp);
   fChain->SetBranchAddress("c44_kp", &c44_kp, &b_c44_kp);
   fChain->SetBranchAddress("c45_kp", &c45_kp, &b_c45_kp);
   fChain->SetBranchAddress("c55_kp", &c55_kp, &b_c55_kp);
   fChain->SetBranchAddress("c11_p", &c11_p, &b_c11_p);
   fChain->SetBranchAddress("c12_p", &c12_p, &b_c12_p);
   fChain->SetBranchAddress("c13_p", &c13_p, &b_c13_p);
   fChain->SetBranchAddress("c14_p", &c14_p, &b_c14_p);
   fChain->SetBranchAddress("c15_p", &c15_p, &b_c15_p);
   fChain->SetBranchAddress("c22_p", &c22_p, &b_c22_p);
   fChain->SetBranchAddress("c23_p", &c23_p, &b_c23_p);
   fChain->SetBranchAddress("c24_p", &c24_p, &b_c24_p);
   fChain->SetBranchAddress("c25_p", &c25_p, &b_c25_p);
   fChain->SetBranchAddress("c33_p", &c33_p, &b_c33_p);
   fChain->SetBranchAddress("c34_p", &c34_p, &b_c34_p);
   fChain->SetBranchAddress("c35_p", &c35_p, &b_c35_p);
   fChain->SetBranchAddress("c44_p", &c44_p, &b_c44_p);
   fChain->SetBranchAddress("c45_p", &c45_p, &b_c45_p);
   fChain->SetBranchAddress("c55_p", &c55_p, &b_c55_p);
   fChain->SetBranchAddress("c11_pim", &c11_pim, &b_c11_pim);
   fChain->SetBranchAddress("c12_pim", &c12_pim, &b_c12_pim);
   fChain->SetBranchAddress("c13_pim", &c13_pim, &b_c13_pim);
   fChain->SetBranchAddress("c14_pim", &c14_pim, &b_c14_pim);
   fChain->SetBranchAddress("c15_pim", &c15_pim, &b_c15_pim);
   fChain->SetBranchAddress("c22_pim", &c22_pim, &b_c22_pim);
   fChain->SetBranchAddress("c23_pim", &c23_pim, &b_c23_pim);
   fChain->SetBranchAddress("c24_pim", &c24_pim, &b_c24_pim);
   fChain->SetBranchAddress("c25_pim", &c25_pim, &b_c25_pim);
   fChain->SetBranchAddress("c33_pim", &c33_pim, &b_c33_pim);
   fChain->SetBranchAddress("c34_pim", &c34_pim, &b_c34_pim);
   fChain->SetBranchAddress("c35_pim", &c35_pim, &b_c35_pim);
   fChain->SetBranchAddress("c44_pim", &c44_pim, &b_c44_pim);
   fChain->SetBranchAddress("c45_pim", &c45_pim, &b_c45_pim);
   fChain->SetBranchAddress("c55_pim", &c55_pim, &b_c55_pim);
   fChain->SetBranchAddress("COS_pi_x", &COS_pi_x, &b_COS_pi_x);
   fChain->SetBranchAddress("COS_pi_y", &COS_pi_y, &b_COS_pi_y);
   fChain->SetBranchAddress("COS_pi_z", &COS_pi_z, &b_COS_pi_z);
   fChain->SetBranchAddress("COS_p_x", &COS_p_x, &b_COS_p_x);
   fChain->SetBranchAddress("COS_p_y", &COS_p_y, &b_COS_p_y);
   fChain->SetBranchAddress("COS_p_z", &COS_p_z, &b_COS_p_z);
   fChain->SetBranchAddress("kplabtheta", &kplabtheta, &b_kplabtheta);
   fChain->SetBranchAddress("kplabphi", &kplabphi, &b_kplabphi);
   fChain->SetBranchAddress("plabtheta", &plabtheta, &b_plabtheta);
   fChain->SetBranchAddress("plabphi", &plabphi, &b_plabphi);
   fChain->SetBranchAddress("pilabtheta", &pilabtheta, &b_pilabtheta);
   fChain->SetBranchAddress("pilabphi", &pilabphi, &b_pilabphi);
   fChain->SetBranchAddress("Xlabtheta", &Xlabtheta, &b_Xlabtheta);
   fChain->SetBranchAddress("Xlabphi", &Xlabphi, &b_Xlabphi);
   fChain->SetBranchAddress("hel", &hel, &b_hel);
   fChain->SetBranchAddress("COSth_pi_x", &COSth_pi_x, &b_COSth_pi_x);
   fChain->SetBranchAddress("COSth_pi_y", &COSth_pi_y, &b_COSth_pi_y);
   fChain->SetBranchAddress("COSth_pi_z", &COSth_pi_z, &b_COSth_pi_z);
   fChain->SetBranchAddress("COSth_p_x", &COSth_p_x, &b_COSth_p_x);
   fChain->SetBranchAddress("COSth_p_y", &COSth_p_y, &b_COSth_p_y);
   fChain->SetBranchAddress("COSth_p_z", &COSth_p_z, &b_COSth_p_z);
   fChain->SetBranchAddress("costheta", &costheta, &b_costheta);
   fChain->SetBranchAddress("COSth_k_z_cm", &COSth_k_z_cm, &b_COSth_k_z_cm);
   fChain->SetBranchAddress("beamPol", &beamPol, &b_beamPol);
   Notify();
}

Bool_t KinFit::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void KinFit::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t KinFit::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef KinFit_cxx
