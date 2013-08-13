#include <TFile.h>
#include <string.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <sstream>
#include <TStyle.h>
#include <TMath.h>
#include <TVector.h>
#include <TGraph.h>

// standard variables
int canvas_width = 800, canvas_height = 600;

//----------------------------------------------------------------------------//
//----------------------------SOME USEFUL FUNCTIONS---------------------------//
//----------------------------------------------------------------------------//

Double_t BackgroundShape (Double_t *xx, Double_t *par) {
    // variable x
    Double_t x = xx[0]/8000.;
    
    // parameters  are nPar, P0, P1, P2, P3, ...
    int nPar = int(par[0]+0.5); // round to int
    Double_t P [nPar];
    for (int i=0; i<=nPar; ++i) {
        P[i] = par[i+1];
    }
    
    // the primary part of the function
    Double_t numerator = P[0] * pow(1-x,P[1]);
    Double_t denominator = pow(x,P[2]+P[3]*log(x));
    Double_t basefunc = numerator / denominator;
    
    // the polynomial (P4 + P5 x + P6 x^2 + ... )
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

Double_t GetResidualSumSquares (TH1D * hist, TF1 * fitfunc, Double_t * pars, 
                                Double_t xmin, Double_t xmax) {
    int nBins = hist->GetXaxis()->GetNbins();
    
    Double_t RSS = 0; 
    for (int i=1; i<=nBins; ++i) {
        Double_t binLow = hist->GetBinLowEdge(i);
        Double_t binHigh = hist->GetBinLowEdge(i+1);
        Double_t binVal = hist->GetBinContent(i);
        Double_t funcInt = fitfunc->Integral(binLow, binHigh, pars);
        Double_t funcAvg = funcInt / (binHigh - binLow);
        
        if ((binLow > xmin) && (binHigh < xmax)) {
            cout << "Bin " << binLow << "-" << binHigh << ": " << binVal << "/" << funcAvg << endl;
            RSS += pow(binVal - funcAvg, 2);
        } else {
            cout << "Bin " << binLow << "-" << binHigh << ": skipped" << endl;
        }
    }
    return RSS;
}

int GetNBinsInRange (TH1D * hist, Double_t xmin, Double_t xmax) {
    int nBins = hist->GetXaxis()->GetNbins();
    int nBinsInRange = 0;    
    
    for (int i=1; i<=nBins; ++i) {
        Double_t binLow = hist->GetBinLowEdge(i);
        Double_t binHigh = hist->GetBinLowEdge(i+1);
        
        if ((binLow > xmin) && (binHigh < xmax)) {
            nBinsInRange += 1;
        }
    }
    return nBinsInRange;
}
                
int NormaliseByBinWidth (TH1D * hist) {
    // divide each bin by its width (assume zero error on bin widths)
    int nbins = hist->GetNbinsX();
    for (int i=1; i<=nbins; ++i) {
        double newbin = hist->GetBinContent(i) / hist->GetBinWidth(i);
        // scale the error too so that relative uncertainty remains untouched
        double newerror = hist->GetBinError(i) / hist->GetBinWidth(i);
        hist->SetBinContent(i, newbin);
        hist->SetBinError(i, newerror);
    }
    return 1;
}

int NormaliseByLuminosity (TH1D * hist, double lum) {
    // divide each bin by the luminosity (assume zero luminosity error)
    hist->Scale(1/lum); // scale keeps relative uncertainty constant
    return 1;
}

//----------------------------------------------------------------------------//
//-----------------------------MAIN ANALYSIS----------------------------------//
//----------------------------------------------------------------------------//

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
    string mjj_histName = "h1_MjjWide_finalSel_varbin;1"; // use VARIABLE BINS
    //string mjj_histName = "h1_MjjWide_finalSel;1"; // use 1 GEV BINS
    
    // luminosities, as given on the twiki
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
    TFile * writefile = new TFile("fTestHists.root","recreate");
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
    
    // choose the range to do the fit in
    Float_t fitRangeLow = 270.;
    Float_t fitRangeHigh = 3000.;
    
    // set plot ranges
    BmjjHist->SetAxisRange(500, 3000, "X");
    CmjjHist->SetAxisRange(500, 3000, "X");
    DmjjHist->SetAxisRange(500, 3000, "X");
    BCDmjjHist->SetAxisRange(fitRangeLow, fitRangeHigh, "X");
    
    // write these histograms to file
    writefile->Write();
    
    // F-test loop
    // repeat until the F-test is satisfied
    int nBasePars = 4;
    int nExtraPars = 0;
    float RSSold = 0;
    float RSS = 0;
    int stopTest = 0;
    while (stopTest == 0) {
        // get number of fit parameters
        // (one of the parameters holds the number of parameters... tricky!)
        int nTotalPars = 1 + nBasePars + nExtraPars;

        // set the name of the fit/histogram
        std::ostringstream s;
        s << "FTest" << nExtraPars;
        string fitNameStr = s.str();
        
        // make a canvas to put fit and residual on
        TCanvas* canvas = new TCanvas((fitNameStr + "-canv").c_str(), (fitNameStr + "-canv").c_str(), 800, 600);
        canvas->cd();
        
        //create 3 pads in the canvas
        TPad* fPads1 = NULL;
        TPad* fPads2 = NULL;

        fPads1 = new TPad("pad1", "", 0.00, 0.20, 0.99, 0.99);
        fPads2 = new TPad("pad3", "", 0.00, 0.00, 0.99, 0.20);
        fPads1->SetFillColor(0);
        fPads1->SetLineColor(0);
        fPads2->SetFillColor(0);
        fPads2->SetLineColor(0);
        fPads1->Draw();
        fPads2->Draw();
        
        // make a new histogram to do the fit on
        fPads1->cd();
        TH1D * fitHist = new TH1D;
        fitHist = BCDmjjHist;
        fitHist->SetName(fitNameStr.c_str());

        // set up the background fit
        TF1 *bkfit = new TF1(fitNameStr.c_str(),BackgroundShape, fitRangeLow, fitRangeHigh, nTotalPars);
        bkfit->SetParName(0,"nPar");
        for (int i=1;i<nTotalPars;++i) {
            std::ostringstream p;
            p << "P" << i-1;
            string parNameStr = p.str();
            bkfit->SetParName(i,parNameStr.c_str());
        }
        
        // set the initial parameter values
        bkfit->FixParameter(0,nTotalPars*1.); // number of parameters does not change!
        bkfit->SetParameter(0+1,2.6e-06); // P0
        bkfit->SetParameter(1+1,5.8); // P1
        bkfit->SetParameter(2+1,6.7); // P2
        bkfit->SetParameter(3+1,0.3); // P3
        for (int i=nBasePars+1;i<nTotalPars;++i) {
            // all others to zero
            bkfit->SetParameter(i+1,0.0);
        }
        
        // Fit and get parameters
        fitHist->Fit(bkfit,"IM"); // MRQ (M=improve, R=use fn range, Q=quiet, I=integral, L=likelihood)
        Double_t * params = bkfit->GetParameters();
        fitHist->Write();
        fitHist->SetDirectory(0);
        
        // Make new TGraph for the residual plot
        fPads2->cd();
        TVectorD nsigma_x(fitHist->GetNbinsX()); 
        TVectorD nsigma_y(fitHist->GetNbinsX()); 
      
        for(int i = 1; i <= fitHist->GetNbinsX(); ++i) {
            double binLow = fitHist->GetXaxis()->GetBinLowEdge(i);
            double binHigh = fitHist->GetXaxis()->GetBinUpEdge(i);
            double exp = fitHist->GetBinContent(i);
            double obs = bkfit->Integral(binLow, binHigh, params) / (binHigh - binLow);
            double error = fitHist->GetBinError(i);
            double x = fitHist->GetBinCenter(i);
            
            double diff = obs - exp;
            double sigma = sqrt(error*error);
        
            if (sigma != 0.0 && obs != 0.0 ) {
                nsigma_x[i] = x;
                nsigma_y[i] = diff / sigma;
            } else {
                nsigma_x[i] = +999999;
                nsigma_y[i] = 0;
            }       
        }
        
        // plot the graph of residuals 
        if (nsigma_x.GetNoElements() != 0 ) {
            TGraph *nsigmaGraph = new TGraph(nsigma_x,nsigma_y);
            nsigmaGraph->SetName((fitNameStr+"-resid").c_str());
            nsigmaGraph->SetTitle("");
            nsigmaGraph->GetYaxis()->SetRangeUser(-5,5);
            nsigmaGraph->GetYaxis()->SetTitle("(Data-Fit)/#sigma");
            nsigmaGraph->GetYaxis()->SetTitleOffset(0.1);
            nsigmaGraph->GetYaxis()->SetTitleSize(0.15);
            nsigmaGraph->GetYaxis()->SetLabelSize(0.07);
            nsigmaGraph->GetXaxis()->SetTitle("");
            nsigmaGraph->GetXaxis()->SetLimits(fitRangeLow, fitRangeHigh);
            nsigmaGraph->GetXaxis()->SetRangeUser(fitRangeLow, fitRangeHigh);
            nsigmaGraph->GetXaxis()->SetTitleOffset(0.01);
            nsigmaGraph->GetXaxis()->SetLabelSize(0.09);
            nsigmaGraph->SetMarkerStyle(8);
            nsigmaGraph->SetMarkerSize(0.8);
            nsigmaGraph->Draw("AP");
            nsigmaGraph->Write();
        }
        
        canvas->Write();
        
        // Get degrees of freedom
        int dof = GetNBinsInRange(fitHist, fitRangeLow, fitRangeHigh) - nTotalPars + 1;
        cout << dof << endl;
        
        // Get the RSS value between fit and histogram
        RSSold = RSS;
        RSS = GetResidualSumSquares(BCDmjjHist, bkfit, params, fitRangeLow, fitRangeHigh);
        cout << nExtraPars << endl;
        cout << RSSold << endl;
        cout << RSS << endl;
        
        if (RSSold==0) {
            // we must do it again with an extra parameter before we can do the f-test
            nExtraPars += 1;
            continue;
        }
        
        // otherwise, do a F-test
        const double alphaFTest=0.05;
        cout << "Starting F-Test evaluation with alpha = " << alphaFTest << endl;
        double fTestVal = (RSSold-RSS) * dof / RSS;
        if (fTestVal < 0) {
            // the fit got WORSE - try again!
            
            // TODO: probably look at this, as at present it will add an extra parameter
            // any time the fit worsens from the last time, even slightly.
            // Really we should only do this if the fit is significantly worse. 
            
            cout << "Fit worsened; automatically trying extra parameter\n";
            nExtraPars += 1;
            continue;
        }
        cout  << "F-test value = " << fTestVal << endl;
        double goodCL =  1.-TMath::FDistI(fTestVal,1,dof);
        cout  << "F-test CL = " << goodCL << endl;
        
        if (goodCL > alphaFTest) {
            // the fit got better but not by much - stop here
            cout << "Requires " << nExtraPars-1 << " extra parameters.\n";
            stopTest = 1;
        } else if (nExtraPars == 10) {
            // 10 extra parameters is enough for anyone
            cout << "Reached 10 extra parameters without success - aborting\n";
            stopTest = 1;
        } else {
            // it got significantly better - try still more parameters
            cout << "F-test failed; trying more parameters\n";
            nExtraPars += 1;
        }
    }
    
    // finish up
    writefile->Close();
    
    return 1;
}