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
std::vector<double> RooCombFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg) {
    // RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    RooRealVar x("x", "x", -1.0, 3.0);
    int nbins = hsig->GetNbinsX();
    x.setBins(nbins);

    TH1D *hpr_clone = (TH1D*)hsig->Clone();
    TH1D *hpi_clone = (TH1D*)hbkg->Clone();
    TH1D *hneg_clone = (TH1D*)hneg->Clone();

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

    RooRealVar scal_sig("scal_sig", "scal_sig", 100, 10, 1000000);
    RooRealVar scal_bkg("scal_bkg", "scal_bkg", 100, 10, 1000000);

    RooDataHist sigDataHist("sigDataHist", "Signal Data", x, RooFit::Import(*hpr_clone));
    RooDataHist bkgDataHist("bkgDataHist", "Background Data", x, RooFit::Import(*hpi_clone));
    RooDataHist negDataHist("negDataHist", "Negative Data", x, RooFit::Import(*hneg_clone));

    RooCrystalBall cb1("cb1", "cb1", x, xs_0, xs_1,xs_2, xs_3, xs_4, xs_5, xs_6);
    RooCrystalBall cb2("cb2", "cb2", x, xb_0, xb_1,xb_2, xb_3, xb_4, xb_5, xb_6);
    RooArgList pdfList(cb1, cb2);
    RooArgList coefList(scal_sig,scal_bkg);

    RooAddPdf combinedPdf("combinedPdf", "Combined PDF", pdfList, coefList);
    RooAbsReal* nllSig = cb1.createNLL(sigDataHist);  // NLL for signal data
    RooAbsReal* nllBkg = cb2.createNLL(bkgDataHist);  // NLL for background data
    RooAbsReal* nllNeg = combinedPdf.createNLL(negDataHist);  // NLL for negative data

    // combine NLL
    RooAddition totalNLL("totalNLL", "Total NLL", RooArgList(*nllSig, *nllBkg, *nllNeg));

    // Use Minuit2 to minimize the NLL
    ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
    minimizer.SetMaxFunctionCalls(100000);
    minimizer.SetMaxIterations(10000);
    minimizer.SetTolerance(0.001);

    ROOT::Math::Functor f([&](const double* params) {
        xs_0.setVal(params[0]);
        xs_1.setVal(params[1]);
        xs_2.setVal(params[2]);
        xs_3.setVal(params[3]);
        xs_4.setVal(params[4]);
        xs_5.setVal(params[5]);
        xs_6.setVal(params[6]);
        xb_0.setVal(params[7]);
        xb_1.setVal(params[8]);
        xb_2.setVal(params[9]);
        xb_3.setVal(params[10]);
        xb_4.setVal(params[11]);
        xb_5.setVal(params[12]);
        xb_6.setVal(params[13]);
        scal_sig.setVal(params[14]);
        scal_bkg.setVal(params[15]);

        return totalNLL.getVal();
    }, 16);

    minimizer.SetFunction(f);

    minimizer.SetLimitedVariable(0, "xs_0", xs_0.getVal(), 0.01, 0.5, 1.3);
    minimizer.SetLimitedVariable(1, "xs_1", xs_1.getVal(), 0.01, 0.05, 1.1);
    minimizer.SetLimitedVariable(2, "xs_2", xs_2.getVal(), 0.01, 0.05, 1.1);
    minimizer.SetLimitedVariable(3, "xs_3", xs_3.getVal(), 0.01, 0.01, 5.0);
    minimizer.SetLimitedVariable(4, "xs_4", xs_4.getVal(), 0.01, 1.0001, 100);
    minimizer.SetLimitedVariable(5, "xs_5", xs_5.getVal(), 0.01, 0.01, 5.0);
    minimizer.SetLimitedVariable(6, "xs_6", xs_6.getVal(), 0.01, 1.0001, 100);
    minimizer.SetLimitedVariable(7, "xb_0", xb_0.getVal(), 0.01, 0.001, 0.35);
    minimizer.SetLimitedVariable(8, "xb_1", xb_1.getVal(), 0.01, 0.05, 1.1);
    minimizer.SetLimitedVariable(9, "xb_2", xb_2.getVal(), 0.01, 0.05, 1.1);
    minimizer.SetLimitedVariable(10, "xb_3", xb_3.getVal(), 0.01, 0.01, 5.0);
    minimizer.SetLimitedVariable(11, "xb_4", xb_4.getVal(), 0.01, 1.0001, 100);
    minimizer.SetLimitedVariable(12, "xb_5", xb_5.getVal(), 0.01, 0.01, 5.0);
    minimizer.SetLimitedVariable(13, "xb_6", xb_6.getVal(), 0.01, 1.0001, 100);
    minimizer.SetLimitedVariable(14, "scal_sig", scal_sig.getVal(), 0.01, 0.1, 1000000);
    minimizer.SetLimitedVariable(15, "scal_bkg", scal_bkg.getVal(), 0.01, 0.1, 1000000);

    minimizer.Minimize();

    const double* xs = minimizer.X();
    const double* errors = minimizer.Errors();
    for (int i = 0; i < 16; ++i) {
        std::cout << "Parameter " << minimizer.VariableName(i) << ": Value = " << xs[i] << ", Error = " << errors[i] << std::endl;
    }

    // TCanvas *c = new TCanvas("c", "Fit Results", 1200, 400);
    // c->Divide(3, 1);

    // RooPlot* frame1 = x.frame();
    // sigDataHist.plotOn(frame1);
    // cb1.plotOn(frame1);

    // RooPlot* frame2 = x.frame();
    // bkgDataHist.plotOn(frame2);
    // cb2.plotOn(frame2);

    // RooPlot* frame3 = x.frame();
    // negDataHist.plotOn(frame3);
    // combinedPdf.plotOn(frame3);

    // c->cd(1);
    // frame1->Draw();

    // c->cd(2);
    // frame2->Draw();

    // c->cd(3);
    // frame3->Draw();

    // c->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/roocomb_crys_results.pdf");
    std::vector<double> params;
    params.push_back(xs[14]);
    params.push_back(errors[14]);
    return params;

}

void fit_start(int seq, int npi, int low, int high){
            // 创建 TStopwatch 对象
    TStopwatch timer;
    // 启动计时器
    timer.Start();
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    TFile *f = new TFile(Form("/data03/tianye/pbar/root/template_bin25/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    TFile *fOutput = new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc_crys/binned_roocomb_seq%04d_npi%04d_%04d.root",seq,npi,low),"recreate");
    TTree *fitParamsTree = new TTree("fitParamsTree", "Fit Parameters Tree");
 
    std::vector<double> fitparams_bin25(2), fitparams_bin100(2);
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
    std::cout <<"sig_yields"<< fitparams_bin25[0] << std::endl;
    std::cout <<"sig_error"<< fitparams_bin25[1] << std::endl;
    }
    fOutput->cd();
    fitParamsTree->Write();
    fOutput->Close();
    f->Close();
    
    std::ofstream timeFile("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/time/RooCombFit_bin25_crys_time.txt", std::ios::app); // 以追加模式打开文件
    if (timeFile.is_open()) {
        timeFile << "RooCombFit_bin25_crys - Real time: " << timer.RealTime() << " seconds" << std::endl;
        timeFile << "RooCombFit_bin25_crys - CPU time: " << timer.CpuTime() << " seconds" << std::endl;
        timeFile.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
   
}
    
int main(int argc, char *argv[]) {
    fit_start(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
    std::cout << "Finished!" << std::endl;
    return 0;

}
  