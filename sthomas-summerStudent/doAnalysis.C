#include <AnalysisClass.C>
#include <TFile.h>
#include <TBrowser.h>

TH2F * myHist;

int doAnalysis () {
    // remove this line to draw graphs to the screen, hopefully
    gROOT->SetBatch();
    
    AnalysisClass t;

    // Extract and match 4-vectors from the data file
    t.Loop();
    MakeOverviewHistograms(&matchedJets, &dsMass, &recoMass);
    
    pair<vector<float>,vector<float> > binEdges;
    binEdges = GenerateBinEdges(30,1200,100,0,2.5,0.5);
    myHist = MakeBinHist(binEdges, &matchedJets);
    AutoRebin(myHist);
    binEdges = GetBinEdges(myHist);
    
    vector<float> mean;
    vector<float> sigma;
    vector<float> pt;
    vector<float> absEta;
    
    GetBinFitStats(myHist, &mean, &sigma, &pt, &absEta, binEdges, true);
    
    return 1;
}