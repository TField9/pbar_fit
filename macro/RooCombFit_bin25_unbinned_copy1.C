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
void RooCombFit(TTree *sigtree, TTree *bkgtree, TTree *negtree, double &nsig,double &nsig_error) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooRealVar x("x", "x", -1.0, 3.0);
    double Ntot = negtree->GetEntries();
    double initial_scal = 1.0;  // 初始为1
    double prev_scal = initial_scal;
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

    RooRealVar scal_sig("scal_sig", "scal_sig", initial_scal, 0.1, 1.0);

    
    Double_t pr_mass2_value, pi_mass2_value, neg_mass2_value;
        sigtree->SetBranchAddress("pr_mass2", &pr_mass2_value);
        if (bkgtree->GetBranch("pi_mass2")) {
            bkgtree->SetBranchAddress("pi_mass2", &pi_mass2_value);
        }
        // else if (bkgtree->GetBranch("pi1_mass2")) {
        //     bkgtree->SetBranchAddress("pi1_mass2", &pi_mass2_value);
        // }
        // else if (bkgtree->GetBranch("pi2_mass2")) {
        //     bkgtree->SetBranchAddress("pi2_mass2", &pi_mass2_value);
        // }
        // else if (bkgtree->GetBranch("pi3_mass2")) {
        //     bkgtree->SetBranchAddress("pi3_mass2", &pi_mass2_value);
        // }
        // else if (bkgtree->GetBranch("pi4_mass2")) {
        //     bkgtree->SetBranchAddress("pi4_mass2", &pi_mass2_value);
        // }
         if(negtree->GetBranch("neg_mass2")) {
            negtree->SetBranchAddress("neg_mass2", &neg_mass2_value);
        }
        else if (negtree->GetBranch("neg1_mass2")) {
            negtree->SetBranchAddress("neg1_mass2", &neg_mass2_value);
        }
        RooDataSet sigData("sigData", "Signal Data", RooArgSet(x));
        RooDataSet bkgData("bkgData", "Background Data", RooArgSet(x));
        RooDataSet negData("negData", "Negative Data", RooArgSet(x));

        for (Long64_t i = 0; i < sigtree->GetEntries(); i++) {
            sigtree->GetEntry(i);
            x.setVal(pr_mass2_value);
            sigData.add(RooArgSet(x));
        }
        for (Long64_t i = 0; i < bkgtree->GetEntries(); i++) {
            bkgtree->GetEntry(i);
            x.setVal(pi_mass2_value);
            bkgData.add(RooArgSet(x));
        }
        for (Long64_t i = 0; i < negtree->GetEntries(); i++) {
            negtree->GetEntry(i);
            x.setVal(neg_mass2_value);
            negData.add(RooArgSet(x));
        }

    expgausexp sigPdf("sigPdf", "Signal PDF", x, sig_mean, sig_sigmaL, sig_alphaL, sig_sigmaR, sig_alphaR);
    expgausexp bkgPdf("bkgPdf", "Background PDF", x, bkg_mean, bkg_sigmaL, bkg_alphaL, bkg_sigmaR, bkg_alphaR);

    RooArgList pdfList(sigPdf, bkgPdf);
    RooArgList coefList(scal_sig);

    RooAddPdf combinedPdf("combinedPdf", "Combined PDF", pdfList, coefList);

    RooAbsReal* nllSig = sigPdf.createNLL(sigData);  // NLL for signal data
    //cout<<nllSig->getVal()<<endl;
    RooAbsReal* nllBkg = bkgPdf.createNLL(bkgData);  // NLL for background data
    //cout<<nllBkg->getVal()<<endl;
    RooAbsReal* nllNeg = combinedPdf.createNLL(negData);  // NLL for negative data
    //cout<<nllNeg->getVal()<<endl;
    // combine NLL
    RooAddition totalNLL("totalNLL", "Total NLL", RooArgList(*nllSig, *nllBkg, *nllNeg));

   const int max_iterations = 3;
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
        if(Ntot==400){
                minimizer.SetVariableLimits(10, 0.3, 0.7);
            }
            else{
                minimizer.SetVariableLimits(10, 0.1, 0.3);
            }
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
                std::cout << "Converged after " << (iter+1) << " iterations." << std::endl;
                break;
            }

            if (iter == max_iterations - 1) {
                std::cout << "Reached max iterations without convergence." << std::endl;
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
        std::cout << "All minimization attempts failed. Using default value." << std::endl;
    }


    
}

void fit_start(int seq, int npi, int low, int high){
    // 创建 TStopwatch 对象
    TStopwatch timer;
    // 启动计时器
    timer.Start();
    // TFile *f = new TFile(Form("/data03/tianye/pbar/root/template_bin25/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    TFile *f = new TFile(Form("/data03/tianye/pbar/root/10percent/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    TFile *fOutput = new TFile(Form("/data03/tianye/pbar/root/sys_result/unbinned_roocomb_seq%04d_npi%04d_%04d.root",seq,npi,low),"recreate");
    TTree *fitParamsTree = new TTree("fitParamsTree", "Fit Parameters Tree");
 
    std::vector<double> nsig(4);
    fitParamsTree->Branch("nsig", &nsig);



    for(int j=low;j<high;j++){
        TTree *sigtree=(TTree *)f->Get(Form("sigtree%04d",j));
        TTree *bkgtree=(TTree *)f->Get(Form("bkgtree%04d",j));
        // TTree *bkgtree1=(TTree *)f->Get(Form("bkgtree1_%04d",j));
        // TTree *bkgtree2=(TTree *)f->Get(Form("bkgtree2_%04d",j));
        // TTree *bkgtree3=(TTree *)f->Get(Form("bkgtree3_%04d",j));
        // TTree *bkgtree4=(TTree *)f->Get(Form("bkgtree4_%04d",j));
        TTree *negtree=(TTree *)f->Get(Form("negtree%04d",j));
        TTree *negtree1=(TTree *)f->Get(Form("negtree1_%04d",j));

        RooCombFit(sigtree, bkgtree, negtree, nsig[0], nsig[2]);
        RooCombFit(sigtree, bkgtree, negtree1, nsig[1], nsig[3]);
        // RooCombFit(sigtree, bkgtree1, negtree, nsig[1], nsig[6]);
        // RooCombFit(sigtree, bkgtree2, negtree, nsig[2], nsig[7]);
        // RooCombFit(sigtree, bkgtree3, negtree, nsig[3], nsig[8]);
        // RooCombFit(sigtree, bkgtree4, negtree, nsig[4], nsig[9]);

        fitParamsTree->Fill();

        std::cout <<"j="<< j << std::endl;
        std::cout <<"sig_yields"<< nsig[0] << "sig_yields"<<nsig[1]<<std::endl;
        //std::cout <<"sig_error"<< nsig[3] << std::endl;
    }
    fOutput->cd();
    fitParamsTree->Write();
    fOutput->Close();
    f->Close();
    timer.Stop();
    // 打印运行时间
    timer.Print();
}
    
int main(int argc, char *argv[]) {
    fit_start(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
    std::cout << "Finished!" << std::endl;
    return 0;

}
  