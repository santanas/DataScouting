#define AnalysisClass_cxx
#include "AnalysisClass.h"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
#include <CBShape.C>
using std::cout;
using std::endl;
using std::pair;
using std::vector; 

//--------------------------------------------------------------------------//
//-------------------------GENERAL SETUP------------------------------------//                                                                        
//--------------------------------------------------------------------------//

// value of pi
const double PI = 2*acos(0.0);

// because I'm a newbie programmer, I do debugging by using print statements
void DEBUG_CHECK() { cout << "*** DEBUG CHECK ***\n"; }

// enable verbose = 1 to print info on events as they are processed
int verbose = 0;

// enable previewMode = 1 to run through only the first 100 events
int previewMode = 0;

// Set up some vectors to hold information gained from looping through the jets.
vector<float> dsMass, recoMass;
vector<TLorentzVector> dsJetList, recoJetList;
vector<pair<TLorentzVector, TLorentzVector> > matchedJets;
vector<float> matchedJetDeltaR;

// Set up general options for plotting and saving histograms
string overviewHistogramFile = "OverviewHists.root";
string binHistogramFile = "BinHists.root";
string fitHistogramFile = "FitHists.root";
string summaryHistogramFile = "meanSigmaHists.root";
int canvasWidth = 800;
int canvasHeight = 600;

//---------------------------------------------------------------------------//
//-----------------USEFUL FUNCTION DEFINITIONS-------------------------------//
//---------------------------------------------------------------------------//

float DeltaR(TLorentzVector *pVec1, TLorentzVector *pVec2) {
    // Returns the value of R (from eta and phi) between two Lorentz vectors
    return sqrt(pow(pVec1->Eta()-pVec2->Eta(),2) + pow(pVec1->Phi()-pVec2->Phi(),2));
}

float DeltaEta(TLorentzVector *pVec1, TLorentzVector *pVec2) {
    // Returns the absolute value of DeltaEta
    return fabs(pVec1->Eta()-pVec2->Eta());
}

float DeltaPhi(TLorentzVector *pVec1, TLorentzVector *pVec2) {
    // Returns the absolute value of DeltaPhi
    return fabs(pVec1->Phi()-pVec2->Phi());
}

float DijetMass(TLorentzVector *pVec1, TLorentzVector *pVec2) {
    // Returns the dijet mass, given Lorentz vectors of two jets
    // calculate (E1+E2)**2 as normal for two scalars
    float energy1 = pVec1->E(), energy2 = pVec2->E();
    float energySum = energy1 + energy2;
    
    // calculate (p1+p2)**2 as dot product with itself
    TVector3 vectorSum = pVec1->Vect() + pVec2->Vect();
    return sqrt(energySum*energySum - vectorSum*vectorSum);
}

float ColSum (TH2F * hist, int xbin) {
    // get the sum of squares of the histogram hist along xbin'th column
    float sum = 0;
    int Nbins = hist->GetYaxis()->GetNbins();
    for (int ybin = 0; ybin<Nbins; ++ybin) {
        sum += hist->GetBinContent(xbin,ybin);
    }
    return sum;
}

float RowSum (TH2F * hist, int ybin) {
    // get the sum of squares of the histogram hist along ybin'th row
    float sum = 0;
    int Nbins = hist->GetXaxis()->GetNbins();
    for (int xbin = 0; xbin<Nbins; ++xbin) {
        sum += hist->GetBinContent(xbin,ybin);
    }
    return sum;
}

int ColNumNonzero (TH2F * hist, int xbin) {
    // get the number of nonzero bins from hist along xbin'th column
    int nonzero = 0;
    int Nbins = hist->GetYaxis()->GetNbins();
    for (int ybin = 0; ybin<Nbins; ++ybin) {
        if (hist->GetBinContent(xbin,ybin) != 0) {
            nonzero++;
        }
    }
    return nonzero;
}

int RowNumNonzero (TH2F * hist, int ybin) {
    // get the number of nonzero bins from hist along ybin'th row
    int nonzero = 0;
    int Nbins = hist->GetXaxis()->GetNbins();
    for (int xbin = 0; xbin<Nbins; ++xbin) {
        if (hist->GetBinContent(xbin,ybin) != 0) {
            nonzero++;
        }
    }
    return nonzero;
}

vector<TLorentzVector> GetFirstJetVectorFromPair (
    vector<pair<TLorentzVector, TLorentzVector> > * pMatchedJets) {
    // Given a vector of pairs, returns just the items representing the first jet vector
    vector<TLorentzVector> outvector;
    for (unsigned int i=0; i<pMatchedJets->size(); ++i) {
        outvector.push_back((pMatchedJets->at(i)).first);
    }
    return outvector; 
}

vector<TLorentzVector> GetSecondJetVectorFromPair (
    vector<pair<TLorentzVector, TLorentzVector> > * pMatchedJets) {
    // Given a vector of pairs, returns just the items representing the second jet vector
    vector<TLorentzVector> outvector;
    for (unsigned int i=0; i<pMatchedJets->size(); ++i) {
        outvector.push_back((pMatchedJets->at(i)).second);
    }
    return outvector; 
}

//----------------------------------------------------------------------------//
//---------------------------CENTRAL LOOP CODE--------------------------------//
//----------------------------------------------------------------------------//

int AnalysisClass::Loop() {
    // Loops over all events and extracts key information contained within.
    // In particular, does the following:
    //   extract lists of DS and RECO jets
    //   applies threshold cuts as defined below
    //   combines into wide jets
    //   matches up each DS jet with a corresponding RECO jet
    //   calculates the dijet mass from the wide jets
    
    if (fChain == 0) { return 0; }

    // Get number of entries
    Long64_t nentries = fChain->GetEntriesFast();
    if (previewMode) { nentries = 100; }
    if (verbose) { cout << "There are " << nentries << " entries.\n"; }

    // Set threshold values at which events are excluded
    float minPtThreshold = 30.; // GeV
    float maxAbsEtaThreshold = 2.5;
    float maxEtaSepThreshold = 2.0;
    float maxRThreshold = 1.1;
    float minPhiSepThreshold = PI/3.;
    float minDeltaRThreshold = 0.1;

    // Enable only the relevant branches to speed up looping (has a minor effect)
    fChain->SetBranchStatus("*",0);
    fChain->SetBranchStatus("n*",1);
    fChain->SetBranchStatus("n*",1);
    fChain->SetBranchStatus("*Pt",1);
    fChain->SetBranchStatus("*Eta",1);
    fChain->SetBranchStatus("*Phi",1);
    fChain->SetBranchStatus("*E",1);
    fChain->SetBranchStatus("*Flag",1);
    
    // Loop through entries
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        // obtain event info and print event number       
        if (verbose) { cout << "\nEvent #" << jentry << endl; }
        Long64_t ientry = LoadTree(jentry); 
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry); 
        nbytes += nb;
        
        // for each event we'll define a vector to hold the information from the DS and RECO jets for comparison   vector<float> recoMass;  
        vector<TLorentzVector> dsJetList, recoJetList;
        
        // Discards (continues on) events that fail the Pt or Eta thresholds for at least one of the two leading jets. Because jets are sorted according to Pt, we just use indices [0] and [1] to access them.
        int atLeastTwoJetsCheck = (nDSJets > 1) && (nRECOJets > 1);
        if (!atLeastTwoJetsCheck) {
            if (verbose) { cout << "Not enough jets\n"; }
            continue; 
        }          
        int meetsMinPtThreshold = ((dsJetPt[0] > minPtThreshold) 
                && (dsJetPt[1] > minPtThreshold) 
                && (recoJetPt[0] > minPtThreshold) 
                && (recoJetPt[1] > minPtThreshold));
        int meetsMaxAbsEtaThreshold = ((fabs(dsJetEta[0]) < maxAbsEtaThreshold) 
                && (fabs(dsJetEta[1]) < maxAbsEtaThreshold)
                && (fabs(recoJetEta[0]) < maxAbsEtaThreshold)
                && (fabs(recoJetEta[1]) < maxAbsEtaThreshold));
        int allFlagsGood = (ECALTPFilterFlag && HBHENoiseFilterResultFlag 
                && hcalLaserEventFilterFlag && eeBadScFilterFlag
                && ECALDeadDRFilterFlag && ECALBoundaryDRFilterFlag);
        int keepEvent = (meetsMinPtThreshold && meetsMaxAbsEtaThreshold && allFlagsGood);
        if (!keepEvent) {
            if (verbose) { cout << "Event discarded\n"; }
            continue; 
        }

        // Get the Lorentz vectors of the two leading jets for each type
        TLorentzVector leadDsJet1, leadDsJet2, leadRecoJet1, leadRecoJet2;
        leadDsJet1.SetPtEtaPhiE(dsJetPt[0],dsJetEta[0],dsJetPhi[0],dsJetE[0]);
        leadDsJet2.SetPtEtaPhiE(dsJetPt[1],dsJetEta[1],dsJetPhi[1],dsJetE[1]);
        leadRecoJet1.SetPtEtaPhiE(recoJetPt[0],recoJetEta[0],recoJetPhi[0],recoJetE[0]);
        leadRecoJet2.SetPtEtaPhiE(recoJetPt[1],recoJetEta[1],recoJetPhi[1],recoJetE[1]);

        // Set up the some wide jets as zero Lorentz vectors 
        // (the TLorentzVector class initialises with zero values)
        TLorentzVector wideDsJet1, wideDsJet2, wideRecoJet1, wideRecoJet2;
        
        if (verbose) { cout << "DS Jets:\n"; }

        // Step through the DS jets and combines into two wide jets
        for (int ijet=0; ijet<nDSJets; ++ijet) {
            // Excludes jets that do not satisfy the Pt or Eta thresholds.
            int jetMinPtThresholdCheck = (dsJetPt[ijet] > minPtThreshold);
            int jetMaxAbsEtaThresholdCheck = (fabs(dsJetEta[ijet]) < maxAbsEtaThreshold);
            int keepJetCheck = (jetMinPtThresholdCheck && jetMaxAbsEtaThresholdCheck);
            if (!keepJetCheck) {
                if (verbose) { cout << "  Jet failed check\n"; }
                continue;
            }

            // For each jet in the event, get its Lorentz vector and 
            // compare to the leading jets to determine which one to add to
            TLorentzVector thisJet;
            thisJet.SetPtEtaPhiE(dsJetPt[ijet], dsJetEta[ijet], dsJetPhi[ijet], dsJetE[ijet]);
            float distToJet1 = DeltaR(&leadDsJet1, &thisJet);
            float distToJet2 = DeltaR(&leadDsJet2, &thisJet); 

            // Only select jets with R < 1.1, as given in the analysis paper
            int maxRThresholdCheck = ((distToJet1 < maxRThreshold) || (distToJet2 < maxRThreshold));
            if (!maxRThresholdCheck) {
                if (verbose) { cout << "  Jet too far away\n"; }
                continue; 
            }

            // Add together the Lorentz vectors of jets close to leading jet
            if (distToJet1 <= distToJet2) { 
                wideDsJet1 += thisJet; 
                if (verbose) { cout<<"  J0<-"<<ijet<<endl;} 
            } else { 
                wideDsJet2 += thisJet; 
                if (verbose) { cout<<"  J1<-"<<ijet<<endl;} 
            }
            // I realise this has an infintesimal bias towards preferring the
            // first jet but realistically the deltaR isn't going to be equal!
            
            // Finally, add to a list of DS jets for matching later.
            dsJetList.push_back(thisJet);            
        }
        
        if (verbose) { cout << "RECO Jets:\n"; }
        
        // Repeat with the RECO jets (code copied directly - TODO functionalise it)
        for (int ijet=0; ijet<nRECOJets; ++ijet) {
            // Excludes jets that do not satisfy the Pt or Eta thresholds.
            int jetMinPtThresholdCheck = (recoJetPt[ijet] > minPtThreshold);
            int jetMaxAbsEtaThresholdCheck = (fabs(recoJetEta[ijet]) < maxAbsEtaThreshold);
            int keepJetCheck = (jetMinPtThresholdCheck && jetMaxAbsEtaThresholdCheck);
            if (!keepJetCheck) {
                if (verbose) { cout << "  Jet failed check\n"; }
                continue;
            }

            // For each jet in the event, get its Lorentz vector and 
            // compare to the leading jets to determine which one to add to
            TLorentzVector thisJet;
            thisJet.SetPtEtaPhiE(recoJetPt[ijet], recoJetEta[ijet], recoJetPhi[ijet], recoJetE[ijet]);
            float distToJet1 = DeltaR(&leadRecoJet1, &thisJet);
            float distToJet2 = DeltaR(&leadRecoJet2, &thisJet); 

            // Only select jets with R < 1.1, as given in the analysis paper
            int maxRThresholdCheck = ((distToJet1 < maxRThreshold) || (distToJet2 < maxRThreshold));
            if (!maxRThresholdCheck) {
                if (verbose) { cout << "  Jet too far away\n"; }
                continue; 
            }

            // Add together the Lorentz vectors of jets close to leading jet
            if (distToJet1 <= distToJet2) { 
                wideRecoJet1 += thisJet; 
                if (verbose) { cout<<"  J0<-"<<ijet<<endl;} 
            } else { 
                wideRecoJet2 += thisJet; 
                if (verbose) { cout<<"  J1<-"<<ijet<<endl;} 
            }
            // I realise this has an infintesimal bias towards preferring the
            // first jet but realistically the deltaR isn't going to be equal!
            
            // Finally, add to a list of RECO jets for matching later.
            recoJetList.push_back(thisJet); 
        }
        
        // JET COMPARISONS
        // at this point let's do a couple of checks to match up the jets
        if (verbose) { cout << "Trying to match " << dsJetList.size() << " DS jets\n"; }
        for (unsigned int iDs=0; iDs<dsJetList.size(); ++iDs) {
            // for each accepted DS jet, pick the closest accepted RECO jet 
            int minDeltaRIndex=0;
            float minDeltaR=-1;
            for (unsigned int iReco=0; iReco<recoJetList.size(); ++iReco) {
                float dR = DeltaR(&dsJetList[iDs], &recoJetList[iReco]);
                if (dR < minDeltaR || minDeltaR < 0) {
                    // store the new minimum R and the index of the matching jet
                    minDeltaR = dR;
                    minDeltaRIndex = iReco; 
                }
            }
            
            // ignore matches with a separation in R that's too large
            if (minDeltaR > minDeltaRThreshold) {
                if (verbose) { cout << "  R separation too large!\n"; }
                continue;
            }
            
            // the closest RECO jet for each DS jet is matched and stored
            pair<TLorentzVector,TLorentzVector> jetPair; // use make_pair instead?
            jetPair.first = dsJetList[iDs];
            jetPair.second = recoJetList[minDeltaRIndex];
            matchedJets.push_back(jetPair);
            if (verbose) { cout << "  Jets matched!\n"; }
        }
        
        // DIJET MASS CALCULATION

        // Calculate pseudorepidity separation and angular phi separation as exclusion criteria
        float etaSep = DeltaEta(&wideDsJet1, &wideDsJet2);
        float phiSep = DeltaPhi(&wideDsJet1, &wideDsJet2);
        int phiSepCheck = (phiSep > minPhiSepThreshold);
        int etaSepCheck = (etaSep < maxEtaSepThreshold);
        int keepDijetCheck = (phiSepCheck && etaSepCheck);
        if (!keepDijetCheck) {
            if (verbose) { cout << "Dijet failed check\n"; }
            continue; 
        }

        // Finally, if all has gone well, calculate the dijet mass!
        float dsDijetMass = DijetMass(&wideDsJet1, &wideDsJet2); 
        float recoDijetMass = DijetMass(&wideRecoJet1, &wideRecoJet2);
        if (verbose) { 
            cout << "mjj_ds = " << dsDijetMass << endl; 
            cout << "mjj_reco = " << dsDijetMass << endl; 
        }

        // Dijet masses less than 400 GeV are discarded
        if (dsDijetMass < 400) {
            if (verbose) { cout << "DS Dijet mass too low\n"; }
            continue; 
        } else {
            // Write the dijet mass to a storage vector
            dsMass.push_back(dsDijetMass);
            // Finish up.
            if (verbose) { cout << "***DS EVENT KEPT***\n"; } 
        }

        // Dijet masses less than 400 GeV are discarded
        if (recoDijetMass < 400) {
            if (verbose) { cout << "RECO Dijet mass too low\n"; }
            continue; 
        } else {
            // Write the dijet mass to a storage vector
            recoMass.push_back(recoDijetMass);
            // Finish up.
            if (verbose) { cout << "***RECO EVENT KEPT***\n"; } 
        }
    }

    return 1;
}

//----------------------------------------------------------------------------//
//----------------------HISTOGRAMS AND SUCH-----------------------------------//
//----------------------------------------------------------------------------//

int MakeOverviewHistograms(vector<pair<TLorentzVector, TLorentzVector> > *pMatchedJets, 
                    vector<float> *pDsDijetMass, vector<float> *pRecoDijetMass) {
    // Makes a variety of histograms which cover the distributions of some jet statistics
    // e.g. distance between matched jets; distributions of Pt and Eta; energy differences
    
    // Open histogram file
    TFile * file = new TFile(overviewHistogramFile.c_str(),"recreate");
    
    // Distribution of jet Pts
    TCanvas * jetPtsCanvas = new TCanvas("jetPtsCanvas", "Jet P_{t} distribution", canvasWidth, canvasHeight);
    TH1F * recoPtsHist = new TH1F("recoPtsHist", "Jet P_{t} distribution", 100, 30, 2500);
    TH1F * dsPtsHist = new TH1F("dsPtsHist", "Jet P_{t} distribution", 100, 30, 2500);
    
    // Dijet masses for each event (unmatched)
    TCanvas * dijetMassCanvas = new TCanvas("dijetMassCanvas", "DS/RECO Dijet Mass", canvasWidth, canvasHeight);
    TH1F * dsMassHist = new TH1F("dsMassHist", "DS/RECO Dijet Mass", 100, 400, 8000);
    TH1F * recoMassHist = new TH1F("recoMassHist", "DS/RECO Dijet Mass", 100, 400, 8000);
    
    // Differences in R for the matched jets
    TH1F * deltaRHist = new TH1F("deltaRHist", "Matched DS/RECO jet #DeltaR/R", 100,0,0.12);
    
    // Differences in Pt
    TH1F * deltaPtHist = new TH1F("deltaPtHist", "Matched DS/RECO jet #DeltaP_{t}/P_{t}", 100,-0.1, 0.1);
    
    // Differences in eta
    TH1F * deltaEtaHist = new TH1F("deltaEtaHist", "Matched DS/RECO jet #Delta#eta/#eta", 100,-0.1,0.1);
    
    // Differences in E
    TH1F * deltaEHist = new TH1F("deltaEHist", "Matched DS/RECO jet #DeltaE/E", 100,-0.1,0.1);
    
    // Differences in rho
    TH1F * deltaRhoHist = new TH1F("deltaRhoHist", "Matched DS/RECO jet #Delta#rho/#rho", 100,-0.1,0.1);
    
    // Differences in phi 
    TH1F * deltaPhiHist = new TH1F("deltaPhiHist", "Matched DS/RECO jet #Delta#phi/#phi", 100,-0.1,0.1);
    
    // 2D histogram for event numbers in pt/|eta| parameter space (DS jets only)
    TH2F * ptEtaHist = new TH2F("ptEtaHist", "p_{t} and |#eta| jet distributions for DS jets", 20, 30, 1200, 20, 0, 3);
    
    // Fill histograms generated from varibles saved once per jet match
    for (unsigned int i=0;i<(pMatchedJets->size());++i) {        
        // copy a few variables to simplify the names
        // I wish I could figure out how to do this using pointers, but I don't think you can take the memory address of the Pt(), Eta() etc methods of the TLorentzVector
        TLorentzVector * pDs = &(pMatchedJets->at(i).first);
        TLorentzVector * pReco = &(pMatchedJets->at(i).second);  
        
        float dsPt = pDs->Pt();
        float recoPt = pReco->Pt();
        float dsEta = pDs->Eta();
        float absDsEta = fabs(dsEta);
        float recoEta = pReco->Eta();
        float dsE = pDs->E();
        float recoE = pReco->E(); 
        float dsRho = pDs->Rho();
        float recoRho = pReco->Rho(); 
        float dsPhi = pDs->Phi();
        float recoPhi = pReco->Phi();
        
        // fill ALL the things
        dsPtsHist->Fill(dsPt);
        recoPtsHist->Fill(recoPt);
        ptEtaHist->Fill(dsPt,absDsEta);
        deltaRHist->Fill(DeltaR(pDs, pReco));
        deltaPtHist->Fill((dsPt - recoPt) / recoPt);
        deltaEtaHist->Fill((dsEta - recoEta) / recoEta);
        deltaRhoHist->Fill((dsRho - recoRho) / recoRho);
        deltaPhiHist->Fill((dsPhi - recoPhi) / recoPhi);
        deltaEHist->Fill((dsE - recoE) / recoE);
    }
    
    // Fill histograms generated from DS dijets
    for (unsigned int i=0;i<(pDsDijetMass->size());++i) {
        dsMassHist->Fill(pDsDijetMass->at(i)); 
    }
    // Fill histograms generated from RECO dijets
    for (unsigned int i=0;i<(pRecoDijetMass->size());++i) {
        recoMassHist->Fill(pRecoDijetMass->at(i)); 
    }
    
    // Set options for histograms that need to be drawn differently
    ptEtaHist->SetOption("colztext");
    
    // Draw canvas with Pt distribution across matched jets
    jetPtsCanvas->cd();    
    recoPtsHist->Draw();
    dsPtsHist->SetLineColor(kRed);
    dsPtsHist->Draw("same");
    TLegend * jetPtsLegend = new TLegend(0.6,0.7,0.89,0.89);
    jetPtsLegend->AddEntry(dsPtsHist,"DS","l");
    jetPtsLegend->AddEntry(recoPtsHist,"RECO","l");
    jetPtsLegend->Draw();
    jetPtsCanvas->Write();
    
    // Draw canvas with dijet mass distribution across all events
    dijetMassCanvas->cd();
    recoMassHist->Draw();
    dsMassHist->SetLineColor(kRed);
    dsMassHist->Draw("same");
    TLegend * dijetMassLegend = new TLegend(0.6,0.7,0.89,0.89);
    dijetMassLegend->AddEntry(dsMassHist,"DS","l");
    dijetMassLegend->AddEntry(recoPtsHist,"RECO","l");
    dijetMassLegend->Draw();
    dijetMassCanvas->Write();
    
    // finish up
    file->Write();    
    file->Close(); 
    return 1;
}

pair<vector<float>,vector<float> > GenerateBinEdges 
    (float ptMin, float ptMax, float ptSize, 
     float absEtaMin, float absEtaMax, float absEtaSize) {
    /* Given limits on the extent of pt and |eta| to search, will generate intial 
     * bins for |eta| and pt.
     * Returns two vectors of bin edges for the pt and eta directions.
     *
     * - ptMin, ptMax, absEtaMin, absEtaMax: define the dimensions of the parameter space
     * - ptSize, absEtaSize: define the starting size 
     * - xlist, ylist: vectors holding bin edges
     * 
     * Returns a pair of vectors holding bin edges */

    vector<float> xlist, ylist;
    float x=ptMin, y=absEtaMin;
    while (x < (ptMax+ptSize)) {
        xlist.push_back(x);
        x += ptSize;
    }
    while (y < (absEtaMax+absEtaSize)) {
        ylist.push_back(y);
        y += absEtaSize;
    }
    
    pair <vector<float>,vector<float> > outbins = make_pair (xlist,ylist);
    return outbins;
}

pair<vector<float>,vector<float> > GetBinEdges (TH2F * myHist) {
    // Return an (x,y) pair containing vector<float>s with the edges of a given TH2F
    TAxis * xaxis = myHist->GetXaxis();
    TAxis * yaxis = myHist->GetYaxis();
    vector<float> xlist, ylist;
    for (int i=1; i < xaxis->GetNbins()+2; ++i) {
        xlist.push_back(xaxis->GetBinLowEdge(i));
    }
    for (int j=1; j < yaxis->GetNbins()+2; ++j) {
        ylist.push_back(yaxis->GetBinLowEdge(j));
    }
    return make_pair(xlist,ylist);
} 

TH2F * MakeBinHist (pair<vector<float>,vector<float> > binedges,
                    vector<pair<TLorentzVector, TLorentzVector> > * pMatchedJets) {
    // Make a histogram of the pt/eta distribution of matched DS jets.
    // binedges: (x,y) pair of vectors containing the edges of the desired bins
    // pMatchedJets: pointer to the matched jet pair list    
    
    // Get bin edges
    vector<float> xBins = binedges.first;
    vector<float> yBins = binedges.second;
    unsigned int nXBins = xBins.size() - 1;
    unsigned int nYBins = yBins.size() - 1;
    float xBinsArray [nXBins];
    float yBinsArray [nYBins];
    for (unsigned int i=0; i<xBins.size(); ++i) {
        xBinsArray[i] = xBins[i];
    }
    for (unsigned int i=0; i<yBins.size(); ++i) {
        yBinsArray[i] = yBins[i];
    }
    
    // Get the matched DS jets only
    vector<TLorentzVector> ds = GetFirstJetVectorFromPair(pMatchedJets);
    
    // Histogram setup and fill (overwrites previous histograms)
    TFile * file = new TFile(binHistogramFile.c_str(),"recreate");
    TH2F * binHist = new TH2F("binHist", "Binning in the pt/|#eta| space", nXBins, xBinsArray, nYBins, yBinsArray);
    binHist->SetOption("colztext");
    for (unsigned int i=0; i<ds.size(); ++i) {
        binHist->Fill(ds.at(i).Pt(), fabs(ds.at(i).Eta()));
    }
    
    // finish up
    file->Write();   
    binHist->SetDirectory(0);
    file->Close();    
    return binHist;
}

int RebinHistogram (TH2F * hist, pair<vector<float>,vector<float> > binedges) {
    // Rebins a 2D histogram according to the bin edges given in binPair
    // the edges have to be the same, bar deletions in the rebinned case
    // you CAN rebin {0,5,10,15,20}x{0,2,4,6,8,10} to {0,10,20}x{0,4,6,10}...
    // but you CAN'T rebin {0,5,10,15,20}x{0,2,4,6,8,10} to {5,15,25}x{0,5,9}, for example.

    // Get bin edges (coped from above - TODO functionalise this?)
    vector<float> xBins = binedges.first;
    vector<float> yBins = binedges.second;
    unsigned int nXBins = xBins.size() - 1;
    unsigned int nYBins = yBins.size() - 1;
    float xBinsArray [nXBins];
    float yBinsArray [nYBins];
    for (unsigned int i=0; i<xBins.size(); ++i) {
        xBinsArray[i] = xBins[i];
    }
    for (unsigned int i=0; i<yBins.size(); ++i) {
        yBinsArray[i] = yBins[i];
    }
    
    // Histogram setup (creates a new file, for now)    
    TFile * file = new TFile(binHistogramFile.c_str(),"recreate");
    TH2F * newHist = new TH2F("newHist", hist->GetTitle(), nXBins, xBinsArray, nYBins, yBinsArray);
    TAxis *xaxis = hist->GetXaxis();
    TAxis *yaxis = hist->GetYaxis();
    
    for (int j=1;j<=yaxis->GetNbins();j++) {
        for (int i=1;i<=xaxis->GetNbins();i++) {
            newHist->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),hist->GetBinContent(i,j));
        }
    }
    
    newHist->SetOption("colztext");
    
    // finish up
    newHist->Write();
    *hist = *newHist;
    hist->SetDirectory(0);
    newHist->SetDirectory(0);
    file->Close();
    return 1;
} 

int AutoRebin (TH2F * hist) {
    // Automatically rebins a TH2F by grouping together columns and rows, alternating between the two. 
    // The decision as to whether to rebin is made by seeing whether the average number of counts in nonzero cells within a given row or column is higher than some threshold value
    
    // threshold values
    float rowSumThreshold = 200;
    float columnSumThreshold = 200;
    
    // initially we mark a rebin as required
    int rebinPt = 1;
    int rebinEta = 1;   
    int rebin = 1;
    
    // if a rebin is required...
    while (rebin == 1) {
        rebin = 0;
        
        // rebin in the Pt (x) direction
        if (rebinPt) {
            // setup variables
            rebinPt = 0;
            vector<float> xEdges, yEdges;
            TAxis * xaxis = hist->GetXaxis();
            TAxis * yaxis = hist->GetYaxis();
            int xNbins = xaxis->GetNbins();
            int yNbins = yaxis->GetNbins();
            if (verbose) { cout << "There are " << xNbins << " columns.\n"; }
            
            // put on the first bin edge
            xEdges.push_back(xaxis->GetBinLowEdge(1));

            for (int x = 1; x < (xNbins + 1); ++x) {
                if (rebinPt) {
                    // if we just rebinned, skip this column completely unless it's the last
                    if (rebinPt == 2) { 
                        rebinPt = 1;
                        if (x == xNbins) {
                            // last column - so make sure we get the upper bin edge
                            xEdges.push_back(xaxis->GetBinUpEdge(x));
                        }
                        if (verbose) { cout << "(Column " << x << " is gone)" << endl; }
                        continue;
                    }
                    // if we have already rebinned on this pass, just copy the bin edge
                    xEdges.push_back(xaxis->GetBinUpEdge(x));
                    if (verbose) { cout << "Skipping column " << x << endl; }
                    rebinPt = 1;
                    continue;
                }

                // print some info
                if (verbose) { 
                    cout << "Column " << x << endl; 
                    cout << "Sum: " << ColSum(hist, x) << endl;
                    cout << "Range: " << xaxis->GetBinLowEdge(x) << "-" << xaxis->GetBinUpEdge(x) << endl;
                }

                // If the average bin number is less than our threshold, needs rebinning
                if (ColSum(hist, x) >= (columnSumThreshold * ColNumNonzero(hist, x))) {
                    // if the column is good, copy the bin edge
                    xEdges.push_back(xaxis->GetBinUpEdge(x)); 
                    if (verbose) { cout << "Keeping column " << x << endl; }
                } else {
                    if (x == xNbins) {
                        // if it's the last bin, merge with the bin before it
                        xEdges.pop_back();
                        xEdges.push_back(xaxis->GetBinUpEdge(x));
                        if (verbose) { cout << "Rebinning last columns " << x-1 << "+" << x << endl; }
                    } else {
                        // otherwise just don't copy across the bin edge
                        if (verbose) { cout << "Rebinning columns " << x << "+" << x+1 << endl; }
                    }
                    rebinPt = 2; // mark for a repeat pass
                }
            }
            // don't change the y axis at this point
            for (int y = 1; y < (yNbins + 2); ++y) {
                yEdges.push_back(yaxis->GetBinLowEdge(y));
            }
            
            // remake the histogram
            pair<vector<float>,vector<float> > newBinEdges = make_pair (xEdges,yEdges);
            RebinHistogram(hist, newBinEdges);
        }

        // Do it again... in the |eta| (y) direction! 
        if (rebinEta) {
            // setup variables
            rebinEta = 0;
            vector<float> xEdges, yEdges;
            TAxis * xaxis = hist->GetXaxis();
            TAxis * yaxis = hist->GetYaxis();
            int xNbins = xaxis->GetNbins();
            int yNbins = yaxis->GetNbins();
            if (verbose) { cout << "There are " << xNbins << " rows.\n"; }
            
            // put on the first bin edge
            yEdges.push_back(yaxis->GetBinLowEdge(1));

            for (int y = 1; y < (yNbins + 1); ++y) {
                if (rebinEta) {
                    // if we just rebinned, skip this row completely unless it's the last
                    if (rebinEta == 2) { 
                        rebinEta = 1;
                        if (y == yNbins) {
                            // last row - so make sure we get the upper bin edge
                            yEdges.push_back(yaxis->GetBinUpEdge(y));
                        } 
                        if (verbose) { cout << "(Row " << y << " is gone)" << endl; }
                        continue;
                    }
                    // if we have already rebinned on this pass, just copy the bin edge
                    yEdges.push_back(yaxis->GetBinUpEdge(y));
                    if (verbose) { cout << "Skipping row " << y << endl; }
                    rebinEta = 1;
                    continue;
                }

                // print some info
                if (verbose) { 
                    cout << "Row " << y << endl; 
                    cout << "Sum: " << RowSum(hist, y) << endl;
                    cout << "Range: " << yaxis->GetBinLowEdge(y) << "-" << yaxis->GetBinUpEdge(y) << endl;
                }

                // If the average bin number is less than our threshold, needs rebinning
                if (RowSum(hist, y) >= (rowSumThreshold * RowNumNonzero(hist, y))) {
                    // if the row is good, copy the bin edge
                    yEdges.push_back(yaxis->GetBinUpEdge(y)); 
                    if (verbose) { cout << "Keeping row " << y << endl; }
                } else {
                    if (y == yNbins) {
                        // if it's the last bin, merge with the bin before it
                        yEdges.pop_back();
                        yEdges.push_back(yaxis->GetBinUpEdge(y));
                        if (verbose) { cout << "Rebinning last rows " << y-1 << "+" << y << endl; }
                    } else {
                        // otherwise just don't copy across the bin edge
                        if (verbose) { cout << "Rebinning rows " << y << "+" << y+1 << endl; }
                    }
                    rebinEta = 2; // mark for a repeat pass
                }
            }
            // don't change the x axis at this point
            for (int x = 1; x < (xNbins + 2); ++x) {
                xEdges.push_back(xaxis->GetBinLowEdge(x));
            }
            
            // remake the histogram
            pair<vector<float>,vector<float> > newBinEdges = make_pair (xEdges,yEdges);
            RebinHistogram(hist, newBinEdges);
        }
        
        // flag for a repeat if either was rebinned last time
        if (rebinPt || rebinEta) { rebin = 1; }
    }
    
    return 1;
}

TF1 * FitShapeToBin(TH2F * hist, vector<pair<TLorentzVector, TLorentzVector> > *pMatchedJets, int xbin, int ybin, bool useCBShape=true) {
    /* For a given bin within the histogram hist, plot the delta pt / pt distribution within that bin.
     * Then fit a Gaussian/Crystal Ball to it and extract its mean and standard deviation.
     * - hist: the histogram with the bin sizes to be calculated.
     * - pMatchedJets: pointer to the jets extracted from AnalysisClass::Loop.
     * - xbin, ybin: x and y numbers of the bin, range [1, nBins] 
     * - useCBShape: if 1, fits using a crystal ball function instead of a gaussian */
    
    // plot options
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(11111111);
    
    // Get the bin limits
    TAxis * xaxis = hist->GetXaxis();
    TAxis * yaxis = hist->GetYaxis();
    float xmid = xaxis->GetBinCenter(xbin);
    float xmax = xaxis->GetBinUpEdge(xbin), xmin = xaxis->GetBinLowEdge(xbin);
    float ymid = yaxis->GetBinCenter(ybin);
    float ymax = yaxis->GetBinUpEdge(ybin), ymin = yaxis->GetBinLowEdge(ybin);
    
    // Set up the histogram (adds to existing file, for now)
    TFile * file = new TFile(fitHistogramFile.c_str(),"update");
    std::ostringstream s;
    s << "Fit(" << xmid << "," << ymid << ")";
    string histNameStr = s.str(); // need to assign .str() method to a temporary variable
    const char * histName = histNameStr.c_str(); // and THEN do the char string conversion
    TCanvas * fitCanvas = new TCanvas(histName, "#DeltaP_{t}/P_{t}", canvasWidth, canvasHeight);  
    fitCanvas->cd();
    TH1F * deltaPtFitHist = new TH1F(histName, "#DeltaP_{t}/P_{t}", 50,-0.3,0.1);
    
    // Step through the jet pairs and fill in the histogram
    for (unsigned int i=0;i<(pMatchedJets->size());++i) {
        // copy a few variables to simplify the names
        TLorentzVector * pDs = &(pMatchedJets->at(i).first);
        TLorentzVector * pReco = &(pMatchedJets->at(i).second);
        
        float dsPt = pDs->Pt();
        float recoPt = pReco->Pt();
        float dsEta = pDs->Eta();
        float absDsEta = fabs(dsEta);
        
        // fill delta Pt / Pt histogram        
        if ((dsPt > xmin) && (dsPt < xmax) && (absDsEta > ymin) && (absDsEta < ymax)) {
            deltaPtFitHist->Fill((dsPt - recoPt) / recoPt);
        }
    }
    
    // get the maximum value for use in the parameters
    float maxVal = deltaPtFitHist->GetBinContent(deltaPtFitHist->GetMaximumBin());
    if (useCBShape) {
        // load the crystal ball function
        TF1 *cbfit = new TF1(histName,CrystalBallFunction,-0.3,0.1,5);
        cbfit->SetParameters(maxVal,1,3,-0.01,0.04);
        cbfit->SetParNames("A","alpha","n","x0","sigma");        
        //cbfit->SetParLimits(0,maxVal*0.9,maxVal*1.1);
        cbfit->SetParLimits(1,0,5);
        cbfit->SetParLimits(3,-0.1,0);
        cbfit->SetParLimits(4,0,0.1);  

        // Fit and get parameters
        deltaPtFitHist->Fit(cbfit,"MRQ"); // MRQ (M=improve fit, R=use fn range, Q=quiet)
        
        // Draw the histogram and fit
        deltaPtFitHist->Draw();
        fitCanvas->Write();
        cbfit->Write();
        
        // finish up
        file->Write();
        file->Close();
        return cbfit;
    } else {
        // Define a Gaussian with some reasonable initial parameters
        TF1 * gausfit = new TF1(histName, "gaus", 0, 0.1);
        gausfit->SetParNames("alpha","x0","sigma");
        gausfit->SetParameters(maxVal, 0, 0.1);
       
        // Fit and get parameters
        deltaPtFitHist->Fit(gausfit);
        
        // Draw the histogram and fit
        deltaPtFitHist->Draw();
        fitCanvas->Write();
        gausfit->Write();
        
        // finish up
        file->Write();
        //file->Close();
        return gausfit;
    }        
}

int GetBinFitStats (TH2F * myHist,
                    vector<float> * mean, vector<float> * sigma, 
                    vector<float> * pt, vector<float> *absEta,
                    pair<vector<float>,vector<float> > binEdges,
                    bool useCBShape=true) {
    // For each of the bins within myHist, fit a defined shape to the pt histogram
    // Store each fit in a file and return some stats about each fit.
    
    // get axes    
    TAxis * xaxis = myHist->GetXaxis();
    TAxis * yaxis = myHist->GetYaxis();
    
    // perform the fit and save to variables
    for (int x=1; x<=xaxis->GetNbins(); ++x){
        for (int y=1; y<=yaxis->GetNbins(); ++y){
            if (myHist->GetBinContent(x,y) > 100) {
                TF1 * newFit = FitShapeToBin(myHist,&matchedJets,x,y,useCBShape);
                mean->push_back(newFit->GetParameter(3));
                sigma->push_back(newFit->GetParameter(4));
                pt->push_back(xaxis->GetBinCenter(x));
                absEta->push_back(yaxis->GetBinCenter(y));
            }
        }
    }

    // Get bin edges - copied, TODO functionalise
    vector<float> xBins = binEdges.first;
    vector<float> yBins = binEdges.second;
    unsigned int nXBins = xBins.size() - 1;
    unsigned int nYBins = yBins.size() - 1;
    float xBinsArray [nXBins];
    float yBinsArray [nYBins];
    for (unsigned int i=0; i<xBins.size(); ++i) {
        xBinsArray[i] = xBins[i];
    }
    for (unsigned int i=0; i<yBins.size(); ++i) {
        yBinsArray[i] = yBins[i];
    }
    
    // open file for writing and make histograms 
    TFile * file = new TFile(summaryHistogramFile.c_str(),"recreate");
    TH2F * meanDiff = new TH2F("meanDiff","Mean",nXBins, xBinsArray, nYBins, yBinsArray);
    TH2F * sigmaError = new TH2F("sigmaError","Sigma",nXBins, xBinsArray, nYBins, yBinsArray);
    meanDiff->SetOption("colztext");
    sigmaError->SetOption("colztext");
    
    // axis labels
    meanDiff->GetXaxis()->SetTitle("p_{T} / GeV");
    meanDiff->GetYaxis()->SetTitle("|#eta|");
    sigmaError->GetXaxis()->SetTitle("p_{T} / GeV");
    sigmaError->GetYaxis()->SetTitle("|#eta|");
    
    // fill histograms
    for (unsigned int i=0; i<pt->size(); ++i) {
        cout << pt->at(i) << " " << absEta->at(i) << "   " << mean->at(i) << " " << sigma->at(i) << " " << endl; 
        meanDiff->Fill(pt->at(i),absEta->at(i),mean->at(i));
        sigmaError->Fill(pt->at(i),absEta->at(i),sigma->at(i));
    }
    
    // finish up
    file->Write();
    file->Close();
    return 1;
}