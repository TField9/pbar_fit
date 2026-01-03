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



std::vector<double> RooComb(TTree *sigtree, TTree *bkgtree, TTree *negtree){
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooRealVar x("x", "x", -1.0, 3.0);

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


    // 定义分箱数量和范围（例如将 x 分成 50 个 bin）
    x.setBins(30); 

    // 将数据集分箱
    RooDataHist sigDataBinned = RooDataHist("sigDataBinned", "Binned Signal Data", RooArgSet(x), sigData);
    RooDataHist bkgDataBinned = RooDataHist("bkgDataBinned", "Binned Background Data", RooArgSet(x), bkgData);
    RooDataHist negDataBinned = RooDataHist("negDataBinned", "Binned Negative Data", RooArgSet(x), negData);    


    expgausexp sigPdf("sigPdf", "Signal PDF", x, sig_mean, sig_sigmaL, sig_alphaL, sig_sigmaR, sig_alphaR);
    expgausexp bkgPdf("bkgPdf", "Background PDF", x, bkg_mean, bkg_sigmaL, bkg_alphaL, bkg_sigmaR, bkg_alphaR);

    RooArgList pdfList(sigPdf, bkgPdf);
    RooArgList coefList(*scal_sig);

    RooAddPdf combinedPdf("combinedPdf", "Combined PDF", pdfList, coefList);

    RooAbsReal* nllSig = sigPdf.createNLL(sigData);  // NLL for signal data
    //cout<<nllSig->getVal()<<endl;
    RooAbsReal* nllBkg = bkgPdf.createNLL(bkgData);  // NLL for background data
    //cout<<nllBkg->getVal()<<endl;
    RooAbsReal* nllNeg = combinedPdf.createNLL(negData);  // NLL for negative data
    //cout<<nllNeg->getVal()<<endl;
    // combine NLL
    RooAddition totalNLL("totalNLL", "Total NLL", RooArgList(*nllSig, *nllBkg, *nllNeg));

    // Use Minuit2 to minimize the NLL
    ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
    minimizer.SetMaxFunctionCalls(1000000);
    minimizer.SetMaxIterations(100000);
    minimizer.SetTolerance(0.001);

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
        scal_sig->setVal(params[10]);
        return totalNLL.getVal();
    }, 11);

    minimizer.SetFunction(f);

    minimizer.SetVariable(0, "sig_mean", sig_mean.getVal(), 0.01);
    minimizer.SetVariable(1, "sig_sigmaL", sig_sigmaL.getVal(), 0.01);
    minimizer.SetVariable(2, "sig_alphaL", sig_alphaL.getVal(), 0.01);
    minimizer.SetVariable(3, "sig_sigmaR", sig_sigmaR.getVal(), 0.01);
    minimizer.SetVariable(4, "sig_alphaR", sig_alphaR.getVal(), 0.01);
    minimizer.SetVariable(5, "bkg_mean", bkg_mean.getVal(), 0.01);
    minimizer.SetVariable(6, "bkg_sigmaL", bkg_sigmaL.getVal(), 0.01);
    minimizer.SetVariable(7, "bkg_alphaL", bkg_alphaL.getVal(), 0.01);
    minimizer.SetVariable(8, "bkg_sigmaR", bkg_sigmaR.getVal(), 0.01);
    minimizer.SetVariable(9, "bkg_alphaR", bkg_alphaR.getVal(), 0.01);
    minimizer.SetVariable(10, "scal_sig", scal_sig->getVal(), 0.01);
    minimizer.SetVariableLimits(0, 0.5, 1.1);
    minimizer.SetVariableLimits(1, 0.05, 1.1);
    minimizer.SetVariableLimits(2, 0.01, 5.0);
    minimizer.SetVariableLimits(3, 0.05, 1.1);
    minimizer.SetVariableLimits(4, 0.01, 6.0);
    minimizer.SetVariableLimits(5, 0.001, 0.35);
    minimizer.SetVariableLimits(6, 0.05, 1.1);
    minimizer.SetVariableLimits(7, 0.01, 5.0);
    minimizer.SetVariableLimits(8, 0.05, 1.1);
    minimizer.SetVariableLimits(9, 0.01, 5.0);
    minimizer.SetVariableLimits(10, 0.0, 1.0);

    minimizer.Minimize();;

    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(0);
    TCanvas *c = new TCanvas("c", "c", 1200, 400);
    c->Divide(3, 1);

    RooPlot* frame1 = x.frame(RooFit::Title("Signal Data"));
    sigData.plotOn(frame1);
    sigPdf.plotOn(frame1, RooFit::LineColor(kRed), RooFit::LineStyle(kDashed),RooFit::Name("sigPdf"));
    frame1->GetXaxis()->SetTitleSize(0.05);
    frame1->GetYaxis()->SetTitleSize(0.05);
    frame1->GetXaxis()->SetTitle("x");
    frame1->GetYaxis()->SetTitle("Count");
    frame1->GetXaxis()->SetLabelFont(43);
    frame1->GetXaxis()->SetLabelSize(20);
    frame1->GetXaxis()->SetTitleOffset(1.1);
    frame1->GetXaxis()->SetNdivisions(505);
    frame1->GetYaxis()->SetLabelFont(43);
    frame1->GetYaxis()->SetLabelSize(15);
    frame1->GetYaxis()->SetLabelOffset(0.01);
    frame1->GetYaxis()->SetNdivisions(505);

    RooPlot* frame2 = x.frame(RooFit::Title("Background Data"));
    bkgData.plotOn(frame2);
    bkgPdf.plotOn(frame2, RooFit::LineColor(kViolet), RooFit::LineStyle(kDashed),RooFit::Name("bkgPdf"));
    frame2->GetXaxis()->SetTitleSize(0.05);
    frame2->GetYaxis()->SetTitleSize(0.05);
    frame2->GetXaxis()->SetTitle("x");
    frame2->GetYaxis()->SetTitle("Count");
    frame2->GetXaxis()->SetLabelFont(43);
    frame2->GetXaxis()->SetLabelSize(20);
    frame2->GetXaxis()->SetTitleOffset(1.1);
    frame2->GetXaxis()->SetNdivisions(505);
    frame2->GetYaxis()->SetLabelFont(43);
    frame2->GetYaxis()->SetLabelSize(20);
    frame2->GetYaxis()->SetLabelOffset(0.02);
    frame2->GetYaxis()->SetNdivisions(505);

    RooPlot* frame3 = x.frame(RooFit::Title("Target Data"));
    negData.plotOn(frame3);
    combinedPdf.plotOn(frame3, RooFit::Name("combinedPdf"));
    combinedPdf.plotOn(frame3, RooFit::Components(sigPdf), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), RooFit::Name("sigComponent"));
    combinedPdf.plotOn(frame3, RooFit::Components(bkgPdf), RooFit::LineColor(kViolet), RooFit::LineStyle(kDashed), RooFit::Name("bkgComponent"));
    auto legend = new TLegend(0.85, 0.75, 0.95, 0.88);
    legend->AddEntry(frame3->findObject("sigComponent"), "Signal", "l");
    legend->AddEntry(frame3->findObject("bkgComponent"), "Background", "l");
    legend->AddEntry(frame3->findObject("combinedPdf"), "Combined", "l");
    legend->SetBorderSize(0);
    legend->Draw();
    frame3->GetXaxis()->SetTitleSize(0.05);
    frame3->GetYaxis()->SetTitleSize(0.05);
    frame3->GetXaxis()->SetTitle("x");
    frame3->GetYaxis()->SetTitle("Count");
    frame3->GetXaxis()->SetLabelFont(43);
    frame3->GetXaxis()->SetLabelSize(20);
    frame3->GetXaxis()->SetTitleOffset(1.1);
    frame3->GetXaxis()->SetNdivisions(505);
    frame3->GetYaxis()->SetLabelFont(43);
    frame3->GetYaxis()->SetLabelSize(20);
    frame3->GetYaxis()->SetLabelOffset(0.02);
    frame3->GetYaxis()->SetNdivisions(505);

    c->cd(1);
    TPad* pad1 = (TPad*)c->cd(1);
    pad1->SetLeftMargin(0.2);
    frame1->Draw();

    c->cd(2);
    TPad* pad2 = (TPad*)c->cd(2);
    pad2->SetLeftMargin(0.2);
    frame2->Draw();

    c->cd(3);
    TPad* pad3 = (TPad*)c->cd(3);
    pad3->SetLeftMargin(0.2);
    frame3->Draw();

    c->SaveAs("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/fit_component.pdf");
    std::vector<double> params;
    params.push_back(scal_sig->getVal());
    params.push_back(scal_sig->getError());
    return params;
}    

void fit_stack(){
TFile *f=new TFile("/data03/tianye/pbar/root/template_bin25/toymc_bin1_seq0001_npi0100.root","read");

TTree *sigtree = (TTree*)f->Get("sigtree9932");
TTree *bkgtree = (TTree*)f->Get("bkgtree9978");
TTree *negtree = (TTree*)f->Get("negtree9932");

std::vector<double> fit_results = RooComb(sigtree, bkgtree, negtree);




}