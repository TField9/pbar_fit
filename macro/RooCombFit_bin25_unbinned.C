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


void RooCombFit(TTree *sigtree, TTree *bkgtree, TTree *negtree,
                double &nsig, double &nsig_error,
                int max_iter = 10, double tol = 0.01)
{
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooRealVar x("x", "x", -1.0, 3.0);
    double Ntot = negtree->GetEntries();

    // 信号 shape 参数
    RooRealVar sig_mean("sig_mean", "sig_mean", 0.8, 0.5, 1.1);
    RooRealVar sig_sigmaL("sig_sigmaL", "sig_sigmaL", 1, 0.05, 1.1);
    RooRealVar sig_alphaL("sig_alphaL", "sig_alphaL", 1, 0.01, 5.0);
    RooRealVar sig_sigmaR("sig_sigmaR", "sig_sigmaR", 1, 0.05, 1.1);
    RooRealVar sig_alphaR("sig_alphaR", "sig_alphaR", 1, 0.01, 6.0);

    // 背景 shape 参数
    RooRealVar bkg_mean("bkg_mean", "bkg_mean", 0.1, 0.001, 0.35);
    RooRealVar bkg_sigmaL("bkg_sigmaL", "bkg_sigmaL", 1, 0.05, 1.1);
    RooRealVar bkg_alphaL("bkg_alphaL", "bkg_alphaL", 1, 0.01, 5.0);
    RooRealVar bkg_sigmaR("bkg_sigmaR", "bkg_sigmaR", 1, 0.05, 1.1);
    RooRealVar bkg_alphaR("bkg_alphaR", "bkg_alphaR", 1, 0.01, 5.0);

    // 比例参数
    RooRealVar *scal_sig;

    scal_sig = new RooRealVar("scal_sig", "scal_sig", 0.5, 0, 1);

    // 数据加载
    Double_t pr_val, pi_val, neg_val;
    if(sigtree->GetBranch("pr_mass2")){
        sigtree->SetBranchAddress("pr_mass2", &pr_val);
    }
    else{
        sigtree->SetBranchAddress("pr1_mass2", &pr_val);
    }



    if (bkgtree->GetBranch("pi_mass2")) {
        bkgtree->SetBranchAddress("pi_mass2",&pi_val);
    }
    else{
        bkgtree->SetBranchAddress("pi1_mass2",&pi_val);
    }
    if (negtree->GetBranch("neg_mass2")) negtree->SetBranchAddress("neg_mass2", &neg_val);
    else                             negtree->SetBranchAddress("neg1_mass2", &neg_val);

    RooDataSet sigData("sigData", "", RooArgSet(x));
    RooDataSet bkgData("bkgData", "", RooArgSet(x));
    RooDataSet negData("negData", "", RooArgSet(x));

    for (Long64_t i = 0; i < sigtree->GetEntries(); ++i) {
        sigtree->GetEntry(i);
        x.setVal(pr_val);  sigData.add(RooArgSet(x));
    }
    for (Long64_t i = 0; i < bkgtree->GetEntries(); ++i) {
        bkgtree->GetEntry(i);
        x.setVal(pi_val);  bkgData.add(RooArgSet(x));
    }
    for (Long64_t i = 0; i < negtree->GetEntries(); ++i) {
        negtree->GetEntry(i);
        x.setVal(neg_val); negData.add(RooArgSet(x));
    }

    // 构造 PDF
    expgausexp sigPdf("sigPdf", "", x, sig_mean, sig_sigmaL, sig_alphaL, sig_sigmaR, sig_alphaR);
    expgausexp bkgPdf("bkgPdf", "", x, bkg_mean, bkg_sigmaL, bkg_alphaL, bkg_sigmaR, bkg_alphaR);
    RooAddPdf combinedPdf("combinedPdf", "", RooArgList(sigPdf, bkgPdf), RooArgList(*scal_sig));

    // 分别创建各类 NLL
    RooAbsReal* nllSig = sigPdf.createNLL(sigData);
    RooAbsReal* nllBkg = bkgPdf.createNLL(bkgData);
    RooAbsReal* nllNeg = combinedPdf.createNLL(negData);
    RooAddition totalNLL("totalNLL", "", RooArgList(*nllSig, *nllBkg, *nllNeg));

    // Minuit2 配置
    ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
    minimizer.SetMaxFunctionCalls(1000000);
    minimizer.SetMaxIterations(100000);
    minimizer.SetTolerance(0.001);
//
    ROOT::Math::Functor func([&](const double* params) {
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
        scal_sig->setVal(params[10]);
        return totalNLL.getVal();
    }, 11);
    minimizer.SetFunction(func);

    // 设定变量初值与边界
    minimizer.SetVariable(0,  "sig_mean",  sig_mean.getVal(),  0.01);
    minimizer.SetVariable(1,  "sig_sigmaL",sig_sigmaL.getVal(),0.01);
    minimizer.SetVariable(2,  "sig_alphaL",sig_alphaL.getVal(),0.01);
    minimizer.SetVariable(3,  "sig_sigmaR",sig_sigmaR.getVal(),0.01);
    minimizer.SetVariable(4,  "sig_alphaR",sig_alphaR.getVal(),0.01);
    minimizer.SetVariable(5,  "bkg_mean",  bkg_mean.getVal(),  0.01);
    minimizer.SetVariable(6,  "bkg_sigmaL",bkg_sigmaL.getVal(),0.01);
    minimizer.SetVariable(7,  "bkg_alphaL",bkg_alphaL.getVal(),0.01);
    minimizer.SetVariable(8,  "bkg_sigmaR",bkg_sigmaR.getVal(),0.01);
    minimizer.SetVariable(9,  "bkg_alphaR",bkg_alphaR.getVal(),0.01);
    minimizer.SetVariable(10, "scal_sig",  scal_sig->getVal(),   0.01);

    // 边界
    minimizer.SetVariableLimits(0, 0.5, 1.1);
    minimizer.SetVariableLimits(1, 0.05,1.1);
    minimizer.SetVariableLimits(2, 0.01,5.0);
    minimizer.SetVariableLimits(3, 0.05,1.1);
    minimizer.SetVariableLimits(4, 0.01,6.0);
    minimizer.SetVariableLimits(5, 0.001,0.35);
    minimizer.SetVariableLimits(6, 0.05,1.1);
    minimizer.SetVariableLimits(7, 0.01,5.0);
    minimizer.SetVariableLimits(8, 0.05,1.1);
    minimizer.SetVariableLimits(9, 0.01,5.0);


    minimizer.SetVariableLimits(10, 0, 1);


    minimizer.Minimize();

    // 最终输出
    const double* xs   = minimizer.X();
    const double* errs = minimizer.Errors();
    nsig       = xs[10] * Ntot;
    nsig_error = errs[10] * Ntot;
    delete nllSig;
    delete nllBkg;
    delete nllNeg;
    delete scal_sig;
}


void fit_start(int npr,int npi){
    // 创建 TStopwatch 对象

    TFile *f = new TFile(Form("/data03/tianye/pbar/root/both_vary/toymc_bin1_seq0001_npi%04d.root",npr),"read");
    TFile *f1 = new TFile(Form("/data03/tianye/pbar/root/both_vary/toymc_bin1_seq0001_npi%04d.root",npi),"read");
    TFile *fOutput = new TFile(Form("/data03/tianye/pbar/root/2d_result/unbinned_roocomb_seq0001_npr%04d_npi%04d.root",npr,npi),"recreate");
    TTree *fitParamsTree = new TTree("fitParamsTree", "Fit Parameters Tree");
    std::cout<<"roocomb unbinned npr npi="<<npr<<npi<<std::endl;
    std::vector<double> nsig(2);
    fitParamsTree->Branch("nsig", &nsig);


    for(int j=0;j<10000;j++){

        TTree *tsigtree = (TTree *)f->Get(Form("sigtree_%04d", j));
        TTree *tbkgtree = (TTree *)f1->Get(Form("bkgtree_%04d", j));
        TTree *tnegtree = (TTree *)f1->Get(Form("negtree_%04d", j));

        RooCombFit(tsigtree, tbkgtree, tnegtree, nsig[0], nsig[1]);

        // 手动删除从文件读取的树对象，释放内存
        delete tsigtree;
        delete tbkgtree;
        delete tnegtree;




        fitParamsTree->Fill();
        if(j%100==0){
            std::cout <<"j="<< j << std::endl;
            std::cout <<"sig_yields"<< nsig[0] << std::endl;
        }
 

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
  