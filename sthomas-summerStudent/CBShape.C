/*_____________________________________________________________________________
* CBShape.C
* Contains the definition of the Crystal Ball function for fitting Gaussians 
* with low-energy tails.
* 
* Originally by Wouter Verkerke & David Kirkby from RooFit: 
* root.cern.ch/viewcvs/trunk/roofit/roofit/src/RooCBShape.cxx?revision=44507 
_____________________________________________________________________________*/

Double_t CrystalBallFunction(Double_t *xx, Double_t *par) {
    
    // variable x
    Double_t x = xx[0];
    
    // parameters alpha, n, x0, sigma
    Double_t A = par[0];
    Double_t alpha = par[1];
    Double_t n = par[2];
    Double_t x0 = par[3];
    Double_t sigma = par[4];
    
    Double_t t = (x-x0)/sigma;
    if (alpha < 0) t = -t;

    Double_t absAlpha = fabs((Double_t)alpha);

    if (t >= -absAlpha) {
        return A*exp(-0.5*t*t);
    }
    else {
        Double_t a = A*TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
        Double_t b = n/absAlpha - absAlpha; 

        return a/TMath::Power(b - t, n);
    }
}