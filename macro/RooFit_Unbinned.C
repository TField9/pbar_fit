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
#include <fstream>
#include <chrono>
#include <thread>
#include <RooFitResult.h>

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

// 对单个 negtree 做组合拟合
void FitNeg(RooAddPdf &combinedPdf, RooRealVar &frac_sig,
            TTree *negtree, double &nsig, double &nsig_error) {
    RooRealVar x("x","x",-1.0,3.0);
    RooDataSet negData("negData","negData",RooArgSet(x));
    Double_t neg_val;
    if (negtree->GetBranch("neg_mass2")) {
        negtree->SetBranchAddress("neg_mass2", &neg_val);
    } else if (negtree->GetBranch("neg1_mass2")) {
        negtree->SetBranchAddress("neg1_mass2", &neg_val);
    }
    for (Long64_t i = 0, n = negtree->GetEntries(); i < n; ++i) {
        negtree->GetEntry(i);
        x.setVal(neg_val);
        negData.add(RooArgSet(x));
    }
    RooFitResult *fitResult = combinedPdf.fitTo(
        negData,
        RooFit::Save(true),
        RooFit::PrintLevel(-1),
        RooFit::Verbose(false),
        RooFit::Minimizer("Minuit2","migrad")
    );
    Int_t Ntot = negData.sumEntries();
    if (fitResult && fitResult->status() == 0) {
        nsig       = frac_sig.getVal() * Ntot;
        nsig_error = frac_sig.getError() * Ntot;
    } else {
        nsig = nsig_error = -1;
    }
}

void expgausexp_roofit(int npr,int npi) {
    TStopwatch timer; timer.Start();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    RooMsgService::instance().silentMode();

    TFile *f = new TFile(Form("/data03/tianye/pbar/root/both_vary/toymc_bin1_seq0001_npi%04d.root",npr),"read");
    TFile *f1 = new TFile(Form("/data03/tianye/pbar/root/both_vary/toymc_bin1_seq0001_npi%04d.root",npi),"read");

    TFile *fOut = TFile::Open(
      Form("/data03/tianye/pbar/root/2d_result/unbinned_roofit_seq0001_npr%04d_npi%04d.root",
           npr,npi),
      "RECREATE"
    );

    std::vector<double> nsig(2);
    TTree *fitParamsTree = new TTree("fitParamsTree","Fit Parameters");
    fitParamsTree->Branch("nsig",      &nsig);
    std::cout<<"roofit unbinned npr npi="<<npr<<npi<<std::endl;
    //fitParamsTree->Branch("nsig_error",&nsig_error);

    for (int j = 0; j < 10000; ++j) {
        // 第一套模板
        TTree *sigtree  = (TTree*)f->Get(Form("sigtree_%04d", j));
        TTree *bkgtree  = (TTree*)f1->Get(Form("bkgtree_%04d", j));
        RooRealVar x("x","x",-1.0,3.0);
        RooDataSet sigData("sigData","",RooArgSet(x)), bkgData("bkgData","",RooArgSet(x));
        Double_t v;
        sigtree->SetBranchAddress("pr_mass2",&v);
        for (Long64_t i=0,n=sigtree->GetEntries(); i<n; ++i) {
            sigtree->GetEntry(i);
            x.setVal(v); sigData.add(RooArgSet(x));
        }
        bkgtree->SetBranchAddress("pi_mass2",&v);
        for (Long64_t i=0,n=bkgtree->GetEntries(); i<n; ++i) {
            bkgtree->GetEntry(i);
            x.setVal(v); bkgData.add(RooArgSet(x));
        }
        RooRealVar sig_mean("sig_mean","",0.8,0.5,1.1),
                   sig_sigmaL("sig_sigmaL","",1,0.05,1.1),
                   sig_alphaL("sig_alphaL","",1,0.01,5.0),
                   sig_sigmaR("sig_sigmaR","",1,0.05,1.1),
                   sig_alphaR("sig_alphaR","",1,0.01,6.0);
        expgausexp f_sig("sigPdf","", x,
                         sig_mean, sig_sigmaL, sig_alphaL,
                         sig_sigmaR, sig_alphaR);
        f_sig.fitTo(sigData,
                    RooFit::Save(), RooFit::PrintLevel(-1),
                    RooFit::Verbose(false),
                    RooFit::Minimizer("Minuit2","migrad"));
        for (auto p:{&sig_mean,&sig_sigmaL,&sig_alphaL,&sig_sigmaR,&sig_alphaR})
            p->setConstant(true);

        RooRealVar bkg_mean("bkg_mean","",0.1,0.001,0.35),
                   bkg_sigmaL("bkg_sigmaL","",1,0.05,1.1),
                   bkg_alphaL("bkg_alphaL","",1,0.01,5.0),
                   bkg_sigmaR("bkg_sigmaR","",1,0.05,1.1),
                   bkg_alphaR("bkg_alphaR","",1,0.01,5.0);
        expgausexp f_bkg("bkgPdf","", x,
                         bkg_mean, bkg_sigmaL, bkg_alphaL,
                         bkg_sigmaR, bkg_alphaR);
        f_bkg.fitTo(bkgData,
                    RooFit::Save(), RooFit::PrintLevel(-1),
                    RooFit::Verbose(false),
                    RooFit::Minimizer("Minuit2","migrad"));
        for (auto p:{&bkg_mean,&bkg_sigmaL,&bkg_alphaL,&bkg_sigmaR,&bkg_alphaR})
            p->setConstant(true);

        RooRealVar frac_sig("frac_sig","",0.4,0.0,1.0);
        RooAddPdf combinedPdf("combinedPdf","",
                              RooArgList(f_sig,f_bkg),
                              RooArgList(frac_sig));

       
        // 对 negtree 与 negtree1 重用两套 PDF
        TTree *negtree  = (TTree*)f1->Get(Form("negtree_%04d", j));

        FitNeg(combinedPdf,  frac_sig,  negtree,  nsig[0], nsig[1]);


        fitParamsTree->Fill();
        delete sigtree;
        delete bkgtree;
        delete negtree;
         
        if(j%100==0){
            std::cout <<"j="<< j << std::endl;
            std::cout <<"sig_yields"<< nsig[0] << std::endl;
        }


    }

    fOut->cd();
    fitParamsTree->Write();
    fOut->Close();

}

int main(int argc, char* argv[]){
    expgausexp_roofit(atoi(argv[1]), atoi(argv[2]));
    std::cout<<"----done-----"<<std::endl;
    return 0;
}
