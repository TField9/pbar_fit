#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TFile.h"
#include <cstdlib>
#include <vector>
#include <string>
#include "TFile.h"
#include "TKey.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TH3F.h"
#include "TMath.h"
#include "RooCrystalBall.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "RooAddPdf.h"
#include <iostream>
#include <map>
#include <stdio.h>
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooMinimizer.h"

//ExpGausExp pdf in RooFit
class expgausexp : public RooAbsPdf {
    public:
        expgausexp() {}
        expgausexp(const char *name, const char *title,
                   RooAbsReal& _x,
                   RooAbsReal& _x0,
                   RooAbsReal& _sigmaL,
                   RooAbsReal& _alphaL,
                   RooAbsReal& _sigmaR,
                   RooAbsReal& _alphaR) :
            RooAbsPdf(name, title),
            x("x", "x", this, _x),
            x0("x0", "x0", this, _x0),
            sigmaL("sigmaL", "sigmaL", this, _sigmaL),
            alphaL("alphaL", "alphaL", this, _alphaL),
            sigmaR("sigmaR", "sigmaR", this, _sigmaR),
            alphaR("alphaR", "alphaR", this, _alphaR) {}

        expgausexp(const expgausexp& other, const char* name = 0) :
            RooAbsPdf(other, name),
            x("x", this, other.x),
            x0("x0", this, other.x0),
            sigmaL("sigmaL", this, other.sigmaL),
            alphaL("alphaL", this, other.alphaL),
            sigmaR("sigmaR", this, other.sigmaR),
            alphaR("alphaR", this, other.alphaR) {}

        virtual TObject* clone(const char* newname) const { return new expgausexp(*this, newname); }
        inline virtual ~expgausexp() {}

        Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const {
            if (matchArgs(allVars, analVars, x)) return 1;
            return 0;
        }

        Double_t analyticalIntegral(Int_t code, const char* rangeName) const {
            assert(code == 1);
            constexpr double sqrtPiOver2 = 1.2533141373;
            constexpr double sqrt2 = 1.4142135624;
            double xmin = x.min(rangeName);
            double xmax = x.max(rangeName);
            double tmin = (xmin - x0) / (xmin < x0 ? sigmaL : sigmaR);
            double tmax = (xmax - x0) / (xmax < x0 ? sigmaL : sigmaR);
            double sum = 0;
            if (tmin < -alphaL) {
                double a = 0.5 * alphaL * alphaL;
                double lv = tmin;
                double uv = std::min(tmax, -alphaL->getVal());
                sum += (sigmaL / alphaL) * (std::exp(a + alphaL * uv) - std::exp(a + alphaL * lv));
            }
            if (tmax > alphaR) {
                double a = 0.5 * alphaR * alphaR;
                double lv = std::max(tmin, alphaR->getVal());
                double uv = tmax;
                sum += (sigmaR / alphaR) * (std::exp(a - alphaR * lv) - std::exp(a - alphaR * uv));
            }
            if (tmin < alphaR && tmax > -alphaL) {
                double sigmaMin = (tmin < double(0)) ? sigmaL : sigmaR;
                double sigmaMax = (tmax < double(0)) ? sigmaL : sigmaR;
                sum += sqrtPiOver2 * (sigmaMax * std::erf(std::min(tmax, alphaR->getVal()) / sqrt2) - sigmaMin * std::erf(std::max(tmin, -alphaL) / sqrt2));
            }
            return sum;
        }

    protected:
        RooRealProxy x;
        RooRealProxy x0;
        RooRealProxy sigmaL;
        RooRealProxy alphaL;
        RooRealProxy sigmaR;
        RooRealProxy alphaR;

        Double_t evaluate() const {
            constexpr double sqrtPiOver2 = 1.2533141373;
            constexpr double sqrt2 = 1.4142135624;

            double t = (x - x0) / (x < x0 ? sigmaL : sigmaR);
            double v = 0;
            if (t < -alphaL) {
                double a = 0.5 * alphaL * alphaL;
                double b = alphaL * t;
                v = std::exp(a + b);
            } else if (t <= alphaR) {
                v = std::exp(-0.5 * t * t);
            } else {
                double a = 0.5 * alphaR * alphaR;
                double b = alphaR * (-t);
                v = std::exp(a + b);
            }

            return v;
        }
    };


//minimizing function, which minimizes the combined NLL
std::vector<double> RooCombFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooRealVar x("x", "x", -1.0, 3.0);
    int nbins = hsig->GetNbinsX();
    x.setBins(nbins);
    RooRealVar sig_mean("sig_mean", "sig_mean", 0.8, 0.5, 1.1);
    RooRealVar sig_sigmaL("sig_sigmaL", "sig_sigmaL", 1, 0.05, 1.1);
    RooRealVar sig_alphaL("sig_alphaL", "sig_alphaL", 1, 0.01, 5.0);
    RooRealVar sig_sigmaR("sig_sigmaR", "sig_sigmaR", 1, 0.05, 1.1);
    RooRealVar sig_alphaR("sig_alphaR", "sig_alphaR", 1, 0.01, 6.0);

    RooRealVar bkg_mean("bkg_mean", "bkg_mean", 0.1, 0.001, 0.35);
    RooRealVar bkg_sigmaL("bkg_sigmaL", "bkg_sigmaL", 1, 0.05, 1.1);
    RooRealVar bkg_alphaL("bkg_alphaL", "bkg_alphaL", 1, 0.01, 5.0);
    RooRealVar bkg_sigmaR("bkg_sigmaR", "bkg_sigmaR", 1, 0.05, 1.1);
    RooRealVar bkg_alphaR("bkg_alphaR", "bkg_alphaR", 1, 0.01, 5.0);

    RooRealVar *scal_sig = new RooRealVar("scal_sig", "scal_sig", 0.1, 0., 1);

    RooDataHist sigDataHist("sigDataHist", "Signal Data", RooArgSet(x), hsig);
    RooDataHist bkgDataHist("bkgDataHist", "Background Data", RooArgSet(x), hbkg);
    RooDataHist negDataHist("negDataHist", "Negative Data", RooArgSet(x), hneg);

    expgausexp sigPdf("sigPdf", "Signal PDF", x, sig_mean, sig_sigmaL, sig_alphaL, sig_sigmaR, sig_alphaR);
    expgausexp bkgPdf("bkgPdf", "Background PDF", x, bkg_mean, bkg_sigmaL, bkg_alphaL, bkg_sigmaR, bkg_alphaR);

    RooArgList pdfList(sigPdf, bkgPdf);
    RooArgList coefList(*scal_sig);

    RooAddPdf combinedPdf("combinedPdf", "Combined PDF", pdfList, coefList);
    RooAbsReal* nllSig = sigPdf.createNLL(sigDataHist);  // NLL for signal data
    RooAbsReal* nllBkg = bkgPdf.createNLL(bkgDataHist);  // NLL for background data
    RooAbsReal* nllNeg = combinedPdf.createNLL(negDataHist);  // NLL for negative data

    // combine NLL
    RooAddition totalNLL("totalNLL", "Total NLL", RooArgList(*nllSig, *nllBkg, *nllNeg));
 // Use RooMinimizer to minimize the NLL
    RooMinimizer minimizer(totalNLL);
    minimizer.setStrategy(1);  
    minimizer.setPrintLevel(-1);  

    // Perform the minimization using MIGRAD
    minimizer.migrad();  
    minimizer.hesse();   


    RooFitResult* fitResult = minimizer.save();  
    
    //Return sig&bkg parameters
    std::vector<double> params;
    params.push_back(sig_mean.getVal());
    params.push_back(sig_sigmaL.getVal());
    params.push_back(sig_alphaL.getVal());
    params.push_back(sig_sigmaR.getVal());
    params.push_back(sig_alphaR.getVal());
    params.push_back(bkg_mean.getVal());
    params.push_back(bkg_sigmaL.getVal());
    params.push_back(bkg_alphaL.getVal());
    params.push_back(bkg_sigmaR.getVal());
    params.push_back(bkg_alphaR.getVal());
    params.push_back(scal_sig->getVal());
    params.push_back(sig_mean.getError());
    params.push_back(sig_sigmaL.getError());
    params.push_back(sig_alphaL.getError());
    params.push_back(sig_sigmaR.getError());
    params.push_back(sig_alphaR.getError());
    params.push_back(bkg_mean.getError());
    params.push_back(bkg_sigmaL.getError());
    params.push_back(bkg_alphaL.getError());
    params.push_back(bkg_sigmaR.getError());
    params.push_back(bkg_alphaR.getError());
    params.push_back(scal_sig->getError());
    return params;
}

void fit_start(int seq, int npi, int low, int high){
    TFile *f = new TFile(Form("/data03/tianye/pbar/root/template_bin25/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    TFile *fOutput = new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc100/binned_roocomb1_seq%04d_npi%04d_%04d.root",seq,npi,low),"recreate");
    TTree *fitParamsTree = new TTree("fitParamsTree", "Fit Parameters Tree");
 
    std::vector<double> fitparams_bin25(22), fitparams_bin100(22);
    fitParamsTree->Branch("fitparams_bin25", &fitparams_bin25);
    fitParamsTree->Branch("fitparams_bin100", &fitparams_bin100);


    for(int j=low;j<high;j++){
    TH1D *hsig_25=(TH1D *)f->Get(Form("data_pr25_%04d",j));
    TH1D *hbkg_25=(TH1D *)f->Get(Form("data_pi25_%04d",j));
    TH1D *hneg_25=(TH1D *)f->Get(Form("data_neg25_%04d",j));
    TH1D *hsig_100=(TH1D *)f->Get(Form("data_pr100_%04d",j));
    TH1D *hbkg_100=(TH1D *)f->Get(Form("data_pi100_%04d",j));
    TH1D *hneg_100=(TH1D *)f->Get(Form("data_neg100_%04d",j));

    fitparams_bin25 = RooCombFit(hsig_25, hbkg_25, hneg_25);
    fitparams_bin100 = RooCombFit(hsig_100, hbkg_100, hneg_100);
    fitParamsTree->Fill();

    std::cout <<"j="<< j << std::endl;
    std::cout <<"sig_yields"<< fitparams_bin25[10] << std::endl;
    std::cout <<"sig_error"<< fitparams_bin25[21] << std::endl;
    }
    fOutput->cd();
    fitParamsTree->Write();
    fOutput->Close();
    f->Close();
}
    
int main(int argc, char *argv[]) {
    fit_start(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
    std::cout << "Finished!" << std::endl;
    return 0;

}
  