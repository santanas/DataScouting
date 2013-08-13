//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun 26 12:38:59 2013 by ROOT version 5.32/00
// from TTree DSComp/
// found on file: ntuple_HLT_vs_RECO__JetHT_Run2012B-v1__370kevents.root
//////////////////////////////////////////////////////////

#ifndef AnalysisClass_h
#define AnalysisClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class AnalysisClass {
    public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    Int_t           nDSJets;
    Float_t         dsJetPt[128];   //[nDSJets]
    Float_t         dsJetEta[128];   //[nDSJets]
    Float_t         dsJetPhi[128];   //[nDSJets]
    Float_t         dsJetE[128];   //[nDSJets]
    Float_t         dsJetFracHad[128];   //[nDSJets]
    Int_t           dsJetMatchIndex[128];   //[nDSJets]
    Float_t         dsRho;
    Float_t         dsMetPt;
    Float_t         dsMetPhi;
    Float_t         dsMetCleanPt;
    Float_t         dsMetCleanPhi;
    Int_t           nRECOJets;
    Float_t         recoJetPt[256];   //[nRECOJets]
    Float_t         recoJetEta[256];   //[nRECOJets]
    Float_t         recoJetPhi[256];   //[nRECOJets]
    Float_t         recoJetE[256];   //[nRECOJets]
    Float_t         recoRho;
    Float_t         recoMetPt;
    Float_t         recoMetPhi;
    Float_t         recoMetCleanPt;
    Float_t         recoMetCleanPhi;
    Char_t          ECALTPFilterFlag;
    Char_t          HBHENoiseFilterResultFlag;
    Char_t          hcalLaserEventFilterFlag;
    Char_t          eeBadScFilterFlag;
    Int_t           ECALDeadDRFilterFlag;
    Int_t           ECALBoundaryDRFilterFlag;

    // List of branches
    TBranch        *b_nDSJets;   //!
    TBranch        *b_dsJetPt;   //!
    TBranch        *b_dsJetEta;   //!
    TBranch        *b_dsJetPhi;   //!
    TBranch        *b_dsJetE;   //!
    TBranch        *b_dsJetFracHad;   //!
    TBranch        *b_dsJetMatchIndex;   //!
    TBranch        *b_dsRho;   //!
    TBranch        *b_dsMetPt;   //!
    TBranch        *b_dsMetPhi;   //!
    TBranch        *b_dsMetCleanPt;   //!
    TBranch        *b_dsMetCleanPhi;   //!
    TBranch        *b_nRECOJets;   //!
    TBranch        *b_recoJetPt;   //!
    TBranch        *b_recoJetEta;   //!
    TBranch        *b_recoJetPhi;   //!
    TBranch        *b_recoJetE;   //!
    TBranch        *b_recoRho;   //!
    TBranch        *b_recoMetPt;   //!
    TBranch        *b_recoMetPhi;   //!
    TBranch        *b_recoMetCleanPt;   //!
    TBranch        *b_recoMetCleanPhi;   //!
    TBranch        *b_ECALTPFilterFlag;   //!
    TBranch        *b_HBHENoiseFilterResultFlag;   //!
    TBranch        *b_hcalLaserEventFilterFlag;   //!
    TBranch        *b_eeBadScFilterFlag;   //!
    TBranch        *b_ECALDeadDRFilterFlag;   //!
    TBranch        *b_ECALBoundaryDRFilterFlag;   //!

    AnalysisClass(TTree *tree=0);
    virtual ~AnalysisClass();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual int      Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalysisClass_cxx
AnalysisClass::AnalysisClass(TTree *tree) : fChain(0) 
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ntuple_HLT_vs_RECO__JetHT_Run2012B-v1__370kevents.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("ntuple_HLT_vs_RECO__JetHT_Run2012B-v1__370kevents.root");
        }
        f->GetObject("DSComp",tree);

    }
    Init(tree);
}

AnalysisClass::~AnalysisClass()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t AnalysisClass::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
    }
    Long64_t AnalysisClass::LoadTree(Long64_t entry)
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

void AnalysisClass::Init(TTree *tree)
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

    fChain->SetBranchAddress("nDSJets", &nDSJets, &b_nDSJets);
    fChain->SetBranchAddress("dsJetPt", dsJetPt, &b_dsJetPt);
    fChain->SetBranchAddress("dsJetEta", dsJetEta, &b_dsJetEta);
    fChain->SetBranchAddress("dsJetPhi", dsJetPhi, &b_dsJetPhi);
    fChain->SetBranchAddress("dsJetE", dsJetE, &b_dsJetE);
    fChain->SetBranchAddress("dsJetFracHad", dsJetFracHad, &b_dsJetFracHad);
    fChain->SetBranchAddress("dsJetMatchIndex", dsJetMatchIndex, &b_dsJetMatchIndex);
    fChain->SetBranchAddress("dsRho", &dsRho, &b_dsRho);
    fChain->SetBranchAddress("dsMetPt", &dsMetPt, &b_dsMetPt);
    fChain->SetBranchAddress("dsMetPhi", &dsMetPhi, &b_dsMetPhi);
    fChain->SetBranchAddress("dsMetCleanPt", &dsMetCleanPt, &b_dsMetCleanPt);
    fChain->SetBranchAddress("dsMetCleanPhi", &dsMetCleanPhi, &b_dsMetCleanPhi);
    fChain->SetBranchAddress("nRECOJets", &nRECOJets, &b_nRECOJets);
    fChain->SetBranchAddress("recoJetPt", recoJetPt, &b_recoJetPt);
    fChain->SetBranchAddress("recoJetEta", recoJetEta, &b_recoJetEta);
    fChain->SetBranchAddress("recoJetPhi", recoJetPhi, &b_recoJetPhi);
    fChain->SetBranchAddress("recoJetE", recoJetE, &b_recoJetE);
    //    fChain->SetBranchAddress("recoJetE", recoJetE, &b_recoJetE);
    fChain->SetBranchAddress("recoRho", &recoRho, &b_recoRho);
    fChain->SetBranchAddress("recoMetPt", &recoMetPt, &b_recoMetPt);
    fChain->SetBranchAddress("recoMetPhi", &recoMetPhi, &b_recoMetPhi);
    fChain->SetBranchAddress("recoMetCleanPt", &recoMetCleanPt, &b_recoMetCleanPt);
    fChain->SetBranchAddress("recoMetCleanPhi", &recoMetCleanPhi, &b_recoMetCleanPhi);
    fChain->SetBranchAddress("ECALTPFilterFlag", &ECALTPFilterFlag, &b_ECALTPFilterFlag);
    fChain->SetBranchAddress("HBHENoiseFilterResultFlag", &HBHENoiseFilterResultFlag, &b_HBHENoiseFilterResultFlag);
    fChain->SetBranchAddress("hcalLaserEventFilterFlag", &hcalLaserEventFilterFlag, &b_hcalLaserEventFilterFlag);
    fChain->SetBranchAddress("eeBadScFilterFlag", &eeBadScFilterFlag, &b_eeBadScFilterFlag);
    fChain->SetBranchAddress("ECALDeadDRFilterFlag", &ECALDeadDRFilterFlag, &b_ECALDeadDRFilterFlag);
    fChain->SetBranchAddress("ECALBoundaryDRFilterFlag", &ECALBoundaryDRFilterFlag, &b_ECALBoundaryDRFilterFlag);
    Notify();
}

Bool_t AnalysisClass::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void AnalysisClass::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
    }
    Int_t AnalysisClass::Cut(Long64_t entry)
    {
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef AnalysisClass_cxx
