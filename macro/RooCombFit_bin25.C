#include <iostream>
#include <vector>
#include <cmath>
#include "TH1D.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooPlot.h"
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
#include "RooFitResult.h"
#include "RooExponential.h"
#include <fstream>
#include <chrono>
#include <thread>
#include "RooFormulaVar.h"


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
void RooCombFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double &nsig, double &nsig_error,  bool create_plots) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooRealVar x("x", "x", -1.0, 3.0);
    int nbins = hsig->GetNbinsX();
    int Ntot = hneg->GetEntries();
    x.setBins(nbins);

    double initial_scal = 1.0;  // 初始为1
    double prev_scal = initial_scal;

    RooRealVar sig_mean("sig_mean", "sig_mean", 0.8, 0.5, 1.1);
    RooRealVar sig_sigmaL("sig_sigmaL", "sig_sigmaL", 0.3, 0.05, 1.1);
    RooRealVar sig_alphaL("sig_alphaL", "sig_alphaL", 1.2, 0.01, 5.0);
    RooRealVar sig_sigmaR("sig_sigmaR", "sig_sigmaR", 0.4, 0.05, 1.1);
    RooRealVar sig_alphaR("sig_alphaR", "sig_alphaR", 0.9, 0.01, 6.0);

    RooRealVar bkg_mean("bkg_mean", "bkg_mean", 0.1, 0.0001, 0.35);
    RooRealVar bkg_sigmaL("bkg_sigmaL", "bkg_sigmaL", 0.2, 0.05, 1.1);
    RooRealVar bkg_alphaL("bkg_alphaL", "bkg_alphaL", 1.5, 0.01, 5.0);
    RooRealVar bkg_sigmaR("bkg_sigmaR", "bkg_sigmaR", 0.3, 0.05, 1.1);
    RooRealVar bkg_alphaR("bkg_alphaR", "bkg_alphaR", 1.2, 0.01, 5.0);
    RooRealVar scal_sig("scal_sig", "scal_sig", initial_scal, 0.1, 1.0);

    RooDataHist sigDataHist("sigDataHist", "Signal Data", RooArgSet(x), hsig);
    RooDataHist bkgDataHist("bkgDataHist", "Background Data", RooArgSet(x), hbkg);
    RooDataHist negDataHist("negDataHist", "Negative Data", RooArgSet(x), hneg);

    expgausexp sigPdf("sigPdf", "Signal PDF", x, sig_mean, sig_sigmaL, sig_alphaL, sig_sigmaR, sig_alphaR);
    expgausexp bkgPdf("bkgPdf", "Background PDF", x, bkg_mean, bkg_sigmaL, bkg_alphaL, bkg_sigmaR, bkg_alphaR);

    RooArgList pdfList(sigPdf, bkgPdf);
    RooArgList coefList(scal_sig);

    RooAddPdf combinedPdf("combinedPdf", "Combined PDF", pdfList, coefList);
    RooAbsReal* nllSig = sigPdf.createNLL(sigDataHist); 
    RooAbsReal* nllBkg = bkgPdf.createNLL(bkgDataHist); 
    RooAbsReal* nllNeg = combinedPdf.createNLL(negDataHist); 
    RooAddition totalNLL("totalNLL", "Total NLL", RooArgList(*nllSig, *nllBkg, *nllNeg));

    const int max_iterations = 10;
    const double convergence_threshold = 0.01;

    nsig = initial_scal * Ntot;
    nsig_error = 0.15 * nsig;
    bool fit_successful = false;

    double current_params[11] = {
        sig_mean.getVal(), sig_sigmaL.getVal(), sig_alphaL.getVal(),
        sig_sigmaR.getVal(), sig_alphaR.getVal(), bkg_mean.getVal(),
        bkg_sigmaL.getVal(), bkg_alphaL.getVal(), bkg_sigmaR.getVal(),
        bkg_alphaR.getVal(), scal_sig.getVal()
    };

    for (int iter = 0; iter < max_iterations; iter++) {
        ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
        minimizer.SetMaxFunctionCalls(10000);
        minimizer.SetMaxIterations(1000);
        minimizer.SetTolerance(0.001);
        minimizer.SetPrintLevel(0);

        ROOT::Math::Functor f([&](const double* params) {
            sig_mean.setVal(params[0]);
            sig_sigmaL.setVal(params[1]);
            sig_alphaL.setVal(params[2]);
            sig_sigmaR.setVal(params[3]);
            sig_alphaR.setVal(params[4]);
            bkg_mean.setVal(params[5]);
            bkg_sigmaL.setVal(params[6]);
            bkg_alphaL.setVal(params[7]);
            bkg_sigmaR.setVal(params[8]);
            bkg_alphaR.setVal(params[9]);
            scal_sig.setVal(params[10]);
            return totalNLL.getVal();
        }, 11);

        minimizer.SetFunction(f);

        for (int i = 0; i < 11; i++) {
            minimizer.SetVariable(i, Form("par%d", i), current_params[i], 0.01);
        }

        minimizer.SetVariableLimits(0, 0.5, 1.1);
        minimizer.SetVariableLimits(1, 0.05, 1.5);
        minimizer.SetVariableLimits(2, 0.01, 8.0);
        minimizer.SetVariableLimits(3, 0.05, 1.5);
        minimizer.SetVariableLimits(4, 0.01, 8.0);
        minimizer.SetVariableLimits(5, 0.001, 0.5);
        minimizer.SetVariableLimits(6, 0.05, 1.5);
        minimizer.SetVariableLimits(7, 0.01, 8.0);
        minimizer.SetVariableLimits(8, 0.05, 1.5);
        minimizer.SetVariableLimits(9, 0.01, 8.0);
        minimizer.SetVariableLimits(10, 0.1, 1.0);

        bool success = minimizer.Minimize();

        if (success) {
            fit_successful = true;
            const double* xs = minimizer.X();
            const double* errors = minimizer.Errors();

            double scal_diff = std::fabs(xs[10] - prev_scal);
            prev_scal = xs[10];

            for (int i = 0; i < 11; i++) {
                current_params[i] = xs[i];
            }

            sig_mean.setVal(xs[0]);
            sig_sigmaL.setVal(xs[1]);
            sig_alphaL.setVal(xs[2]);
            sig_sigmaR.setVal(xs[3]);
            sig_alphaR.setVal(xs[4]);
            bkg_mean.setVal(xs[5]);
            bkg_sigmaL.setVal(xs[6]);
            bkg_alphaL.setVal(xs[7]);
            bkg_sigmaR.setVal(xs[8]);
            bkg_alphaR.setVal(xs[9]);
            scal_sig.setVal(xs[10]);

            nsig = xs[10] * Ntot;
            nsig_error = errors[10] * Ntot;

            if (scal_diff < convergence_threshold) {
                //std::cout << "Converged after " << (iter+1) << " iterations." << std::endl;
                break;
            }

            if (iter == max_iterations - 1) {
                //std::cout << "Reached max iterations without convergence." << std::endl;
            }
        } else {
            

            prev_scal *= 0.8;
            if (prev_scal < 0.1) prev_scal = 0.1;
            current_params[10] = prev_scal;

            if (!fit_successful) {
                nsig = -1;
                nsig_error = -1;
            }
        }
    }

    if (!fit_successful) {
        //std::cout << "All minimization attempts failed. Using default value." << std::endl;
    }

 

    delete nllSig;
    delete nllBkg;
    delete nllNeg;
}



void fit_start(int num_pr, int num_pi){
    // 创建 TStopwatch 对象
    TStopwatch timer;
    // 启动计时器
    timer.Start();
    TFile *f = new TFile(Form("/data03/tianye/pbar/root/both_vary/toymc_bin1_seq0001_npi%04d.root",num_pr),"read");
    TFile *f1 = new TFile(Form("/data03/tianye/pbar/root/both_vary/toymc_bin1_seq0001_npi%04d.root",num_pi),"read");
    TFile *fOutput = new TFile(Form("/data03/tianye/pbar/root/2d_result_check/binned_roocomb_npr%04d_npi%04d.root",num_pr,num_pi),"recreate");
    TTree *fitParamsTree = new TTree("fitParamsTree", "Fit Parameters Tree");
    std::vector<double> nsig(2);
    fitParamsTree->Branch("nsig", &nsig);
    std::cout<<"numpr="<<num_pr<<" numpi="<<num_pi<<std::endl;
    
    for(int j=0; j<10000; j++){
        
        TH1D *hsig = (TH1D *)f->Get(Form("data_pr30_%04d", j));
        TH1D *hbkg = (TH1D *)f1->Get(Form("data_pi30_%04d", j));
        TH1D *hneg = (TH1D *)f1->Get(Form("data_neg30_%04d",j));

        // 两个拟合使用相同的迭代优化
        RooCombFit(hsig, hbkg, hneg, nsig[0], nsig[1],  false);

        


        fitParamsTree->Fill();



  

        // 计算并显示剩余时间
        if(j%100==0){
            double elapsedTime = timer.RealTime();
            double progress = static_cast<double>(j + 1) / (10000);
            double estimatedTotalTime = elapsedTime / progress;
            double remainingTime = estimatedTotalTime - elapsedTime;
            std::cout << "Estimated remaining time: " << remainingTime << " seconds" << std::endl;
            std::cout << "j=" << j << std::endl;
            std::cout << "sig_yield1: " << nsig[0] << std::endl;
            std::cout << "sig_error1: " << nsig[1] << std::endl;
        }
        hsig->Delete();
        hbkg->Delete();
        hneg->Delete();
    }
    fOutput->cd();
    fitParamsTree->Write();
    fOutput->Close();
    f->Close();
}  

int main(int argc, char *argv[]) {
    fit_start(atoi(argv[1]), atoi(argv[2]));
    std::cout << "Finished!" << std::endl;
    return 0;
}