#include <TFile.h>
#include <string.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <sstream>
#include <TStyle.h>

Double_t BackgroundShape (Double_t *xx, Double_t *par) {
    // variable x
    Double_t x = xx[0]/8000.;
    
    // parameters nPar, P0, P1, P2, P3, ...
    int nPar = int(par[0]+0.5); // round to int
    Double_t P [nPar];
    for (int i=0; i<=nPar; ++i) {
        P[i] = par[i+1];
    }
    
    Double_t numerator = P[0] * pow(1-x,P[1]);
    Double_t denominator = pow(x,P[2]+P[3]*log(x));
    Double_t basefunc = numerator / denominator;
    
    Double_t polynomial = 1.0;
    
    if (nPar == 4) {
        return basefunc;
    } else {
        for (int pNum = 4; pNum<nPar; ++pNum) {
            polynomial = polynomial + P[pNum]*pow(x,pNum-3);
        }
    return basefunc * polynomial;
    }
}

Double_t GetResidualSumSquares (TH1D * hist, TF1 * fitfunc, Double_t * pars) {
    int nBins = hist->GetXaxis()->GetNbins();
    Double_t RSS = 0; 
    for (int i=1; i<=nBins; ++i) {
        Double_t binLow = hist->GetBinLowEdge();
        Double_t binHigh = hist->GetBinUpEdge();
        Double_t binVal = hist->GetBinContents(i);
        funcInt = fitfunc->Integral(binLow, binHigh, pars);
        RSS += pow(binVal - funcInt, 2);
    }
    return RSS;
}                
                
                
int NormaliseByBinWidth (TH1D * hist) {
    // divide each bin by its width (assume zero error on bin widths)
    int nbins = hist->GetNbinsX();
    for (int i=1; i<=nbins; ++i) {
        double newbin = hist->GetBinContent(i) / hist->GetBinWidth(i);
        double newerror = hist->GetBinError(i) / hist->GetBinWidth(i);
        hist->SetBinContent(i, newbin);
        hist->SetBinError(i, newerror);
    }
    return 1;
}

int NormaliseByLuminosity (TH1D * hist, double lum) {
    // divide each bin by the luminosity (assume zero luminosity error)
    hist->Scale(1/lum);
    return 1;
}

int doHistAnalysis () {
    
    // plot options
    gStyle->SetOptFit(111111111);
    gStyle->SetOptStat(11111111);
    
    // paths to files
    string Bfile_path = "rootfiles/DataScouting_V00-01-06_Run2012B_runrange_193752-197044_dijet_alfaT_razor_dijetpairs_trijetpairs.root";
    string Cfile_path = "rootfiles/DataScouting_V00-01-06_Run2012C_runrange_197885-203755_dijet_alfaT_razor_dijetpairs_trijetpairs.root";
    string Dfile_path = "rootfiles/DataScouting_V00-01-06_Run2012D_runrange_203773-208686_dijet_alfaT_razor_dijetpairs_trijetpairs.root";
    string BCDfile_path = "rootfiles/DataScouting_V00-01-06_Run2012B_Run2012C_Run2012D_runrange_193752-208686_dijet_alfaT_razor_dijetpairs_trijetpairs.root";
    string histPath = "DQMData/Run 999999/DataScouting/Run summary/DiJet";
    //string mjj_histName = "h1_MjjWide_finalSel_varbin;1";
    string mjj_histName = "h1_MjjWide_finalSel;1";
    
    double BLum_fb = 4.445, CLum_fb = 6.806, DLum_fb = 7.384; // in inverse fb
    double BCDLum_fb = BLum_fb + CLum_fb + DLum_fb;
    double BLum_pb = BLum_fb * 1e3;
    double CLum_pb = CLum_fb * 1e3;
    double DLum_pb = DLum_fb * 1e3;
    double BCDLum_pb = BCDLum_fb * 1e3; // converts to inverse pb
    
    // read in the files and get histograms
    TFile * Bfile = new TFile (Bfile_path.c_str());
    TDirectory * BHistDir = Bfile->GetDirectory(histPath.c_str());
    TH1D * BmjjHist = (TH1D*)BHistDir->Get(mjj_histName.c_str());
    BmjjHist->SetName("RunB dijet mass");
    BmjjHist->SetTitle("RunB dijet mass");
    
    TFile * Cfile = new TFile (Cfile_path.c_str());
    TDirectory * CHistDir = Cfile->GetDirectory(histPath.c_str());
    TH1D * CmjjHist = (TH1D*)CHistDir->Get(mjj_histName.c_str());
    CmjjHist->SetName("RunC dijet mass");
    CmjjHist->SetTitle("RunC dijet mass");
    
    TFile * Dfile = new TFile (Dfile_path.c_str());
    TDirectory * DHistDir = Dfile->GetDirectory(histPath.c_str());
    TH1D * DmjjHist = (TH1D*)DHistDir->Get(mjj_histName.c_str());
    DmjjHist->SetName("RunD dijet mass");
    DmjjHist->SetTitle("RunD dijet mass");
    
    TFile * BCDfile = new TFile (BCDfile_path.c_str());
    TDirectory * BCDHistDir = BCDfile->GetDirectory(histPath.c_str());
    TH1D * BCDmjjHist = (TH1D*)BCDHistDir->Get(mjj_histName.c_str());
    BCDmjjHist->SetName("RunsBCD dijet mass");
    BCDmjjHist->SetTitle("RunsBCD dijet mass");
    
    // read out into memory and close the files
    TFile * writefile = new TFile("outputHists.root","recreate");
    BmjjHist->SetDirectory(writefile);
    CmjjHist->SetDirectory(writefile);
    DmjjHist->SetDirectory(writefile);
    BCDmjjHist->SetDirectory(writefile);
    Bfile->Close();
    Cfile->Close();
    Dfile->Close();
    BCDfile->Close();
    
    // normalise histograms
    NormaliseByLuminosity(BmjjHist, BLum_pb);
    NormaliseByBinWidth(BmjjHist);
    NormaliseByLuminosity(CmjjHist, CLum_pb);
    NormaliseByBinWidth(CmjjHist);
    NormaliseByLuminosity(DmjjHist, DLum_pb);
    NormaliseByBinWidth(DmjjHist);
    NormaliseByLuminosity(BCDmjjHist, BCDLum_pb);
    NormaliseByBinWidth(BCDmjjHist);
    
    // make ratio histograms
    TH1D * BCratio = (TH1D*)BmjjHist->Clone();
    BCratio->Divide(CmjjHist);
    BCratio->SetName("B/C ratio");
    BCratio->SetTitle("B/C ratio");
    
    TH1D * BDratio = (TH1D*)BmjjHist->Clone();
    BDratio->Divide(DmjjHist);
    BDratio->SetName("B/D ratio");
    BDratio->SetTitle("B/D ratio");
    
    TH1D * CDratio = (TH1D*)CmjjHist->Clone();
    CDratio->Divide(DmjjHist);
    CDratio->SetName("C/D ratio");
    CDratio->SetTitle("C/D ratio");
    
    // set plot options
    BmjjHist->SetAxisRange(500, 3000, "X");
    CmjjHist->SetAxisRange(500, 3000, "X");
    DmjjHist->SetAxisRange(500, 3000, "X");
    BCDmjjHist->SetAxisRange(270, 3000, "X");
    
    // add extra parameters to the fit
    // an extra parameter is added to hold the number of parameters. tricky
    int nBasePars = 4;
    int nExtraPars = 3;
    int nTotalPars = 1 + nBasePars + nExtraPars;
    
    std::ostringstream s;
    s << "FTest" << nExtraPars;
    string fitNameStr = s.str();

    // set up the background fit
    TF1 *bkfit = new TF1(fitNameStr.c_str(),BackgroundShape, 270,3000, nTotalPars);
    bkfit->SetParName(0,"nPar");
    for (int i=1;i<nTotalPars;++i) {
        std::ostringstream p;
        p << "P" << i-1;
        string parNameStr = p.str();
        bkfit->SetParName(i,parNameStr.c_str());
    }
    
    // set the initial parameter values
    bkfit->FixParameter(0,nTotalPars*1.); // number of parameters does not change!
    bkfit->SetParameter(0+1,5.6e-06); // P0
    bkfit->SetParameter(1+1,7.9); // P1
    bkfit->SetParameter(2+1,6.3); // P2
    bkfit->SetParameter(3+1,0.3); // P3
    for (int i=nBasePars+1;i<nTotalPars;++i) {
        // all others to zero
        bkfit->SetParameter(i+1,0.0);
    }
    
    // make a canvas for the fit
    TCanvas * fitCanvas = new TCanvas(fitNameStr.c_str(), fitNameStr.c_str(),800, 600);  
    
    // Fit and get parameters
    BCDmjjHist->Fit(bkfit,"I"); // MRQ (M=improve fit, R=use fn range, Q=quiet)
    Double_t * params = bkfit->GetParameters();
    fitCanvas->Write();
    BCDmjjHist->Draw();
    
    GetResidual
    
    // write those histograms to file
    writefile->Write();
    BCDmjjHist->SetDirectory(0);
    writefile->Close();
    
    return 1;
}