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


std::vector<double> RooFit_(TTree *sigtree, TTree *bkgtree, TTree *negtree){
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooRealVar x("x", "x", -1.0, 3.0);

    RooRealVar xs_0("xs_0", "xs_0",0.8, 0.5, 1.3);
    RooRealVar xs_1("xs_1", "sigmaL",1, 0.05, 1.1);
    RooRealVar xs_2("xs_2", "sigmaR",1, 0.05, 1.1);
    RooRealVar xs_3("xs_3", "alphaL",1, 0.01, 5.0);
    RooRealVar xs_4("xs_4", "nL",2,1.0001, 100);
    RooRealVar xs_5("xs_5", "alphaR",1, 0.01, 5.0);
    RooRealVar xs_6("xs_6", "nR",2,1.0001, 100);

    RooRealVar xb_0("xb_0", "xb_0",0.001, 0.001, 0.35);
    RooRealVar xb_1("xb_1", "sigmaL",1, 0.05, 1.1);
    RooRealVar xb_2("xb_2", "sigmaR",1, 0.05, 1.1);
    RooRealVar xb_3("xb_3", "alphaL",1, 0.01, 5.0);
    RooRealVar xb_4("xb_4", "nL",2,1.0001, 100);
    RooRealVar xb_5("xb_5", "alphaR",1, 0.01, 5.0);
    RooRealVar xb_6("xb_6", "nR",2,1.0001, 100);

    RooRealVar scal_sig("scal_sig", "scal_sig", 100, 0, 1000000);
    RooRealVar scal_bkg("scal_bkg", "scal_bkg", 100, 0, 1000000);

    Double_t pr_mass2_value, pi_mass2_value, neg_mass2_value;
        sigtree->SetBranchAddress("pr_mass2", &pr_mass2_value);
        bkgtree->SetBranchAddress("pi_mass2", &pi_mass2_value);
        negtree->SetBranchAddress("neg_mass2", &neg_mass2_value);

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

    RooCrystalBall cb1("cb1", "cb1", x, xs_0, xs_1,xs_2, xs_3, xs_4, xs_5, xs_6);
    RooCrystalBall cb2("cb2", "cb2", x, xb_0, xb_1,xb_2, xb_3, xb_4, xb_5, xb_6);

    xs_0.setConstant(kFALSE);
    xs_1.setConstant(kFALSE);
    xs_2.setConstant(kFALSE);
    xs_3.setConstant(kFALSE);
    xs_4.setConstant(kFALSE);
    xs_5.setConstant(kFALSE);
    xs_6.setConstant(kFALSE);

    xb_0.setConstant(kFALSE);
    xb_1.setConstant(kFALSE);
    xb_2.setConstant(kFALSE);
    xb_3.setConstant(kFALSE);
    xb_4.setConstant(kFALSE);
    xb_5.setConstant(kFALSE);
    xb_6.setConstant(kFALSE);

    cb1.fitTo(sigData, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Verbose(kFALSE));
    cb2.fitTo(bkgData, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Verbose(kFALSE));

    xs_0.setConstant(kTRUE);
    xs_1.setConstant(kTRUE);
    xs_2.setConstant(kTRUE);
    xs_3.setConstant(kTRUE);
    xs_4.setConstant(kTRUE);
    xs_5.setConstant(kTRUE);
    xs_6.setConstant(kTRUE);

    xb_0.setConstant(kTRUE);
    xb_1.setConstant(kTRUE);
    xb_2.setConstant(kTRUE);
    xb_3.setConstant(kTRUE);
    xb_4.setConstant(kTRUE);
    xb_5.setConstant(kTRUE);
    xb_6.setConstant(kTRUE);

    RooAddPdf model("model", "model", RooArgList(cb1, cb2), RooArgList(scal_sig, scal_bkg));

    // // Fit the combined PDF to the negative sample histogram data
    model.fitTo(negData,RooFit::Save(),RooFit::PrintLevel(-1000),RooFit::Verbose(kFALSE));
    // TCanvas* c = new TCanvas("c", "Fit Results", 1200, 400);
    // c->Divide(3, 1);

    // RooPlot* frame1 = x.frame();
    // sigData.plotOn(frame1);
    // cb1.plotOn(frame1);

    // RooPlot* frame2 = x.frame();
    // bkgData.plotOn(frame2);
    // cb2.plotOn(frame2);

    // RooPlot* frame3 = x.frame();
    // negData.plotOn(frame3);
    // model.plotOn(frame3);

    // c->cd(1);
    // frame1->Draw();

    // c->cd(2);
    // frame2->Draw();

    // c->cd(3);
    // frame3->Draw();

    // c->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/roofit_results.pdf");
    std::vector<double> params;
    params.push_back(scal_sig.getVal());
    params.push_back(scal_sig.getError());
    return params;
}


void expgausexp_roofit(int seq, int npi, int low, int high){
        TStopwatch timer;
    // 启动计时器
    timer.Start();
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
    TFile *f = new TFile(Form("/data03/tianye/pbar/root/template_bin25/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    TFile *fOutput = new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc_crys/unbinned_roofit_seq%04d_npi%04d_%04d.root",seq,npi,low),"recreate");

    std::vector<double> fitparams;
    TTree *fitParamsTree = new TTree("fitParamsTree", "Fit Parameters Tree");
    fitParamsTree->Branch("fitparams_unbinned", &fitparams);
    for(int j=low;j<high;j++){
    TTree *sigtree=(TTree *)f->Get(Form("sigtree%04d",j));
    TTree *bkgtree=(TTree *)f->Get(Form("bkgtree%04d",j));
    TTree *negtree=(TTree *)f->Get(Form("negtree%04d",j));

    fitparams=RooFit_(sigtree,bkgtree,negtree);

    fitParamsTree->Fill();


    fOutput->cd();

    std::cout<<"j = "<<j<<std::endl;
    std::cout<<"sig_yields = "<<fitparams[0]<<std::endl;
    std::cout<<"sig_error = "<<fitparams[1]<<std::endl;
    }
    fOutput->cd();
    fitParamsTree->Write();
    fitParamsTree->Delete();
    fOutput->Close();

    std::ofstream timeFile("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/time/RooFit_Unbinned_crys_time.txt", std::ios::app); // 以追加模式打开文件
    if (timeFile.is_open()) {
        timeFile << "RooFit_Unbinned_crys - Real time: " << timer.RealTime() << " seconds" << std::endl;
        timeFile << "RooFit_Unbinned_crys - CPU time: " << timer.CpuTime() << " seconds" << std::endl;
        timeFile.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
   
   
}

int main(int argc,char* argv[]){
    expgausexp_roofit(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
    return 0;
}