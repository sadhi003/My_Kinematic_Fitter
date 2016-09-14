//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 31 19:59:57 2016 by ROOT version 5.34/05
// from TTree output/data
// found on file: kppmpimCut_v4.root
//////////////////////////////////////////////////////////

#ifndef KinFitTest_h
#define KinFitTest_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class KinFitTest {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         runno;
   Float_t         eventno;
   Float_t         ebeam_no_ecor;
   Float_t         kp_no_p;
   Float_t         kp_no_p_e;
   Float_t         kp_no_px;
   Float_t         kp_no_py;
   Float_t         kp_no_pz;
   Float_t         p_no_p;
   Float_t         p_no_p_e;
   Float_t         p_no_px;
   Float_t         p_no_py;
   Float_t         p_no_pz;
   Float_t         miss_no_p;
   Float_t         miss_no_p_e;
   Float_t         miss_no_px;
   Float_t         miss_no_py;
   Float_t         miss_no_pz;
   Float_t         ebeam;
   Float_t         w;
   Float_t         kpp;
   Float_t         kpE;
   Float_t         kppx;
   Float_t         kppy;
   Float_t         kppz;
   Float_t         pp;
   Float_t         pE;
   Float_t         ppx;
   Float_t         ppy;
   Float_t         ppz;
   Float_t         missp;
   Float_t         missE;
   Float_t         misspx;
   Float_t         misspy;
   Float_t         misspz;
   Float_t         misspion_p;
   Float_t         misspion_e;
   Float_t         misspion_px;
   Float_t         misspion_py;
   Float_t         misspion_pz;
   Float_t         mm2pkp;
   Float_t         mm2p;
   Float_t         mm2kp;
   Float_t         mm2ppim;
   Float_t         mm2kppim;
   Float_t         invm2kppim;
   Float_t         invm2ppim;
   Float_t         invm2kpp;
   Float_t         s;
   Float_t         t;
   Float_t         u;
   Float_t         missp2;
   Float_t         misspZ2;
   Float_t         misspT2;
   Float_t         scv;
   Float_t         stv;
   Float_t         vtime;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         sctkp;
   Float_t         sctp;
   Float_t         sclkp;
   Float_t         sclp;
   Float_t         tofbetakp;
   Float_t         tofbetap;
   Float_t         betakp;
   Float_t         betap;
   Float_t         tofmasskp;
   Float_t         tofmassp;
   Float_t         goodSTkp;
   Float_t         goodSTp;
   Float_t         stid_kp;
   Float_t         stid_p;
   Float_t         nphoton_nb0;
   Float_t         nphoton_nb1;
   Float_t         nphoton_nb2;
   Float_t         trackPositive;
   Float_t         trackNegative;
   Float_t         trackNeutral;
   Float_t         NKp;
   Float_t         NKm;
   Float_t         NPip;
   Float_t         NPim;
   Float_t         NGamma;
   Float_t         tkpp;
   Float_t         COSth_kp_lab;
   Float_t         Phi_kp_lab;
   Float_t         COSth_p_lab;
   Float_t         Phi_p_lab;
   Float_t         COSth_pim_lab;
   Float_t         Phi_pim_lab;
   Float_t         hel;
   Float_t         trigbit;
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
   Float_t         kpfidL;
   Float_t         kpfidN;
   Float_t         kpfidT;
   Float_t         pfidL;
   Float_t         pfidN;
   Float_t         pfidT;
   Float_t         kptofKO;
   Float_t         ptofKO;
   Float_t         theta_kp;
   Float_t         theta_p;
   Float_t         phi_kp;
   Float_t         phi_p;

   // List of branches
   TBranch        *b_runno;   //!
   TBranch        *b_eventno;   //!
   TBranch        *b_ebeam_no_ecor;   //!
   TBranch        *b_kp_no_p;   //!
   TBranch        *b_kp_no_p_e;   //!
   TBranch        *b_kp_no_px;   //!
   TBranch        *b_kp_no_py;   //!
   TBranch        *b_kp_no_pz;   //!
   TBranch        *b_p_no_p;   //!
   TBranch        *b_p_no_p_e;   //!
   TBranch        *b_p_no_px;   //!
   TBranch        *b_p_no_py;   //!
   TBranch        *b_p_no_pz;   //!
   TBranch        *b_miss_no_p;   //!
   TBranch        *b_miss_no_p_e;   //!
   TBranch        *b_miss_no_px;   //!
   TBranch        *b_miss_no_py;   //!
   TBranch        *b_miss_no_pz;   //!
   TBranch        *b_ebeam;   //!
   TBranch        *b_w;   //!
   TBranch        *b_kpp;   //!
   TBranch        *b_kpE;   //!
   TBranch        *b_kppx;   //!
   TBranch        *b_kppy;   //!
   TBranch        *b_kppz;   //!
   TBranch        *b_pp;   //!
   TBranch        *b_pE;   //!
   TBranch        *b_ppx;   //!
   TBranch        *b_ppy;   //!
   TBranch        *b_ppz;   //!
   TBranch        *b_missp;   //!
   TBranch        *b_missE;   //!
   TBranch        *b_misspx;   //!
   TBranch        *b_misspy;   //!
   TBranch        *b_misspz;   //!
   TBranch        *b_misspion_p;   //!
   TBranch        *b_misspion_e;   //!
   TBranch        *b_misspion_px;   //!
   TBranch        *b_misspion_py;   //!
   TBranch        *b_misspion_pz;   //!
   TBranch        *b_mm2pkp;   //!
   TBranch        *b_mm2p;   //!
   TBranch        *b_mm2kp;   //!
   TBranch        *b_mm2ppim;   //!
   TBranch        *b_mm2kppim;   //!
   TBranch        *b_invm2kppim;   //!
   TBranch        *b_invm2ppim;   //!
   TBranch        *b_invm2kpp;   //!
   TBranch        *b_s;   //!
   TBranch        *b_t;   //!
   TBranch        *b_u;   //!
   TBranch        *b_missp2;   //!
   TBranch        *b_misspZ2;   //!
   TBranch        *b_misspT2;   //!
   TBranch        *b_scv;   //!
   TBranch        *b_stv;   //!
   TBranch        *b_vtime;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_sctkp;   //!
   TBranch        *b_sctp;   //!
   TBranch        *b_sclkp;   //!
   TBranch        *b_sclp;   //!
   TBranch        *b_tofbetakp;   //!
   TBranch        *b_tofbetap;   //!
   TBranch        *b_betakp;   //!
   TBranch        *b_betap;   //!
   TBranch        *b_tofmasskp;   //!
   TBranch        *b_tofmassp;   //!
   TBranch        *b_goodSTkp;   //!
   TBranch        *b_goodSTp;   //!
   TBranch        *b_stid_kp;   //!
   TBranch        *b_stid_p;   //!
   TBranch        *b_nphoton_nb0;   //!
   TBranch        *b_nphoton_nb1;   //!
   TBranch        *b_nphoton_nb2;   //!
   TBranch        *b_trackPositive;   //!
   TBranch        *b_trackNegative;   //!
   TBranch        *b_trackNeutral;   //!
   TBranch        *b_NKp;   //!
   TBranch        *b_NKm;   //!
   TBranch        *b_NPip;   //!
   TBranch        *b_NPim;   //!
   TBranch        *b_NGamma;   //!
   TBranch        *b_tkpp;   //!
   TBranch        *b_COSth_kp_lab;   //!
   TBranch        *b_Phi_kp_lab;   //!
   TBranch        *b_COSth_p_lab;   //!
   TBranch        *b_Phi_p_lab;   //!
   TBranch        *b_COSth_pim_lab;   //!
   TBranch        *b_Phi_pim_lab;   //!
   TBranch        *b_hel;   //!
   TBranch        *b_trigbit;   //!
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
   TBranch        *b_kpfidL;   //!
   TBranch        *b_kpfidN;   //!
   TBranch        *b_kpfidT;   //!
   TBranch        *b_pfidL;   //!
   TBranch        *b_pfidN;   //!
   TBranch        *b_pfidT;   //!
   TBranch        *b_kptofKO;   //!
   TBranch        *b_ptofKO;   //!
   TBranch        *b_theta_kp;   //!
   TBranch        *b_theta_p;   //!
   TBranch        *b_phi_kp;   //!
   TBranch        *b_phi_p;   //!

   KinFitTest(TTree *tree=0);
   virtual ~KinFitTest();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef KinFitTest_cxx
KinFitTest::KinFitTest(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/kppmpimAfterCut.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../data/kppmpimAfterCut.root");
      }
      f->GetObject("output",tree);

   }
   Init(tree);
}

KinFitTest::~KinFitTest()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t KinFitTest::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t KinFitTest::LoadTree(Long64_t entry)
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

void KinFitTest::Init(TTree *tree)
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

   fChain->SetBranchAddress("runno", &runno, &b_runno);
   fChain->SetBranchAddress("eventno", &eventno, &b_eventno);
   fChain->SetBranchAddress("ebeam_no_ecor", &ebeam_no_ecor, &b_ebeam_no_ecor);
   fChain->SetBranchAddress("kp_no_p", &kp_no_p, &b_kp_no_p);
   fChain->SetBranchAddress("kp_no_p_e", &kp_no_p_e, &b_kp_no_p_e);
   fChain->SetBranchAddress("kp_no_px", &kp_no_px, &b_kp_no_px);
   fChain->SetBranchAddress("kp_no_py", &kp_no_py, &b_kp_no_py);
   fChain->SetBranchAddress("kp_no_pz", &kp_no_pz, &b_kp_no_pz);
   fChain->SetBranchAddress("p_no_p", &p_no_p, &b_p_no_p);
   fChain->SetBranchAddress("p_no_p_e", &p_no_p_e, &b_p_no_p_e);
   fChain->SetBranchAddress("p_no_px", &p_no_px, &b_p_no_px);
   fChain->SetBranchAddress("p_no_py", &p_no_py, &b_p_no_py);
   fChain->SetBranchAddress("p_no_pz", &p_no_pz, &b_p_no_pz);
   fChain->SetBranchAddress("miss_no_p", &miss_no_p, &b_miss_no_p);
   fChain->SetBranchAddress("miss_no_p_e", &miss_no_p_e, &b_miss_no_p_e);
   fChain->SetBranchAddress("miss_no_px", &miss_no_px, &b_miss_no_px);
   fChain->SetBranchAddress("miss_no_py", &miss_no_py, &b_miss_no_py);
   fChain->SetBranchAddress("miss_no_pz", &miss_no_pz, &b_miss_no_pz);
   fChain->SetBranchAddress("ebeam", &ebeam, &b_ebeam);
   fChain->SetBranchAddress("w", &w, &b_w);
   fChain->SetBranchAddress("kpp", &kpp, &b_kpp);
   fChain->SetBranchAddress("kpE", &kpE, &b_kpE);
   fChain->SetBranchAddress("kppx", &kppx, &b_kppx);
   fChain->SetBranchAddress("kppy", &kppy, &b_kppy);
   fChain->SetBranchAddress("kppz", &kppz, &b_kppz);
   fChain->SetBranchAddress("pp", &pp, &b_pp);
   fChain->SetBranchAddress("pE", &pE, &b_pE);
   fChain->SetBranchAddress("ppx", &ppx, &b_ppx);
   fChain->SetBranchAddress("ppy", &ppy, &b_ppy);
   fChain->SetBranchAddress("ppz", &ppz, &b_ppz);
   fChain->SetBranchAddress("missp", &missp, &b_missp);
   fChain->SetBranchAddress("missE", &missE, &b_missE);
   fChain->SetBranchAddress("misspx", &misspx, &b_misspx);
   fChain->SetBranchAddress("misspy", &misspy, &b_misspy);
   fChain->SetBranchAddress("misspz", &misspz, &b_misspz);
   fChain->SetBranchAddress("misspion_p", &misspion_p, &b_misspion_p);
   fChain->SetBranchAddress("misspion_e", &misspion_e, &b_misspion_e);
   fChain->SetBranchAddress("misspion_px", &misspion_px, &b_misspion_px);
   fChain->SetBranchAddress("misspion_py", &misspion_py, &b_misspion_py);
   fChain->SetBranchAddress("misspion_pz", &misspion_pz, &b_misspion_pz);
   fChain->SetBranchAddress("mm2pkp", &mm2pkp, &b_mm2pkp);
   fChain->SetBranchAddress("mm2p", &mm2p, &b_mm2p);
   fChain->SetBranchAddress("mm2kp", &mm2kp, &b_mm2kp);
   fChain->SetBranchAddress("mm2ppim", &mm2ppim, &b_mm2ppim);
   fChain->SetBranchAddress("mm2kppim", &mm2kppim, &b_mm2kppim);
   fChain->SetBranchAddress("invm2kppim", &invm2kppim, &b_invm2kppim);
   fChain->SetBranchAddress("invm2ppim", &invm2ppim, &b_invm2ppim);
   fChain->SetBranchAddress("invm2kpp", &invm2kpp, &b_invm2kpp);
   fChain->SetBranchAddress("s", &s, &b_s);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("u", &u, &b_u);
   fChain->SetBranchAddress("missp2", &missp2, &b_missp2);
   fChain->SetBranchAddress("misspZ2", &misspZ2, &b_misspZ2);
   fChain->SetBranchAddress("misspT2", &misspT2, &b_misspT2);
   fChain->SetBranchAddress("scv", &scv, &b_scv);
   fChain->SetBranchAddress("stv", &stv, &b_stv);
   fChain->SetBranchAddress("vtime", &vtime, &b_vtime);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("sctkp", &sctkp, &b_sctkp);
   fChain->SetBranchAddress("sctp", &sctp, &b_sctp);
   fChain->SetBranchAddress("sclkp", &sclkp, &b_sclkp);
   fChain->SetBranchAddress("sclp", &sclp, &b_sclp);
   fChain->SetBranchAddress("tofbetakp", &tofbetakp, &b_tofbetakp);
   fChain->SetBranchAddress("tofbetap", &tofbetap, &b_tofbetap);
   fChain->SetBranchAddress("betakp", &betakp, &b_betakp);
   fChain->SetBranchAddress("betap", &betap, &b_betap);
   fChain->SetBranchAddress("tofmasskp", &tofmasskp, &b_tofmasskp);
   fChain->SetBranchAddress("tofmassp", &tofmassp, &b_tofmassp);
   fChain->SetBranchAddress("goodSTkp", &goodSTkp, &b_goodSTkp);
   fChain->SetBranchAddress("goodSTp", &goodSTp, &b_goodSTp);
   fChain->SetBranchAddress("stid_kp", &stid_kp, &b_stid_kp);
   fChain->SetBranchAddress("stid_p", &stid_p, &b_stid_p);
   fChain->SetBranchAddress("nphoton_nb0", &nphoton_nb0, &b_nphoton_nb0);
   fChain->SetBranchAddress("nphoton_nb1", &nphoton_nb1, &b_nphoton_nb1);
   fChain->SetBranchAddress("nphoton_nb2", &nphoton_nb2, &b_nphoton_nb2);
   fChain->SetBranchAddress("trackPositive", &trackPositive, &b_trackPositive);
   fChain->SetBranchAddress("trackNegative", &trackNegative, &b_trackNegative);
   fChain->SetBranchAddress("trackNeutral", &trackNeutral, &b_trackNeutral);
   fChain->SetBranchAddress("NKp", &NKp, &b_NKp);
   fChain->SetBranchAddress("NKm", &NKm, &b_NKm);
   fChain->SetBranchAddress("NPip", &NPip, &b_NPip);
   fChain->SetBranchAddress("NPim", &NPim, &b_NPim);
   fChain->SetBranchAddress("NGamma", &NGamma, &b_NGamma);
   fChain->SetBranchAddress("tkpp", &tkpp, &b_tkpp);
   fChain->SetBranchAddress("COSth_kp_lab", &COSth_kp_lab, &b_COSth_kp_lab);
   fChain->SetBranchAddress("Phi_kp_lab", &Phi_kp_lab, &b_Phi_kp_lab);
   fChain->SetBranchAddress("COSth_p_lab", &COSth_p_lab, &b_COSth_p_lab);
   fChain->SetBranchAddress("Phi_p_lab", &Phi_p_lab, &b_Phi_p_lab);
   fChain->SetBranchAddress("COSth_pim_lab", &COSth_pim_lab, &b_COSth_pim_lab);
   fChain->SetBranchAddress("Phi_pim_lab", &Phi_pim_lab, &b_Phi_pim_lab);
   fChain->SetBranchAddress("hel", &hel, &b_hel);
   fChain->SetBranchAddress("trigbit", &trigbit, &b_trigbit);
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
   fChain->SetBranchAddress("kpfidL", &kpfidL, &b_kpfidL);
   fChain->SetBranchAddress("kpfidN", &kpfidN, &b_kpfidN);
   fChain->SetBranchAddress("kpfidT", &kpfidT, &b_kpfidT);
   fChain->SetBranchAddress("pfidL", &pfidL, &b_pfidL);
   fChain->SetBranchAddress("pfidN", &pfidN, &b_pfidN);
   fChain->SetBranchAddress("pfidT", &pfidT, &b_pfidT);
   fChain->SetBranchAddress("kptofKO", &kptofKO, &b_kptofKO);
   fChain->SetBranchAddress("ptofKO", &ptofKO, &b_ptofKO);
   fChain->SetBranchAddress("theta_kp", &theta_kp, &b_theta_kp);
   fChain->SetBranchAddress("theta_p", &theta_p, &b_theta_p);
   fChain->SetBranchAddress("phi_kp", &phi_kp, &b_phi_kp);
   fChain->SetBranchAddress("phi_p", &phi_p, &b_phi_p);
   Notify();
}

Bool_t KinFitTest::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void KinFitTest::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t KinFitTest::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef KinFitTest_cxx
