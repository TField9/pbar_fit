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
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TRandom.h"
#include "TFile.h"
#include <cstdlib>
#include <string>
#include "TSystem.h"
#include "TH3F.h"
#include "TMath.h"
#include "RooCrystalBall.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "RooAddPdf.h"
#include <map>
#include <stdio.h>
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"




void setParameterWithDebug(TF1* func, int parIndex, double value, double low, double up) {
    if (value < low || value > up) {
        std::cerr << func->GetName()<<" "<<"Warning: Parameter " << parIndex << " value " << value << " is outside the range [" << low << ", " << up << "]. Setting to (" << low << " + " << up << ") / 2." << std::endl;
        value = (low + up) / 2;
    }
    func->SetParameter(parIndex, value);
    func->SetParLimits(parIndex, low, up);
}

double funcExpGausExp(double *x, double *p)
{
  constexpr double sqrtPiOver2 = 1.2533141373;
  constexpr double sqrt2 = 1.4142135624;

  const double &x0 = p[0];
  const double &sigmaL = p[1];
  const double &alphaL = p[2];
  const double &sigmaR = p[3];
  const double &alphaR = p[4];
  const double &xmin = p[5];
  const double &xmax = p[6];

  double t = (x[0] - x0) / (x[0] < x0 ? sigmaL : sigmaR);
  double v = 0;
  if (t < -alphaL)
  {
    double a = 0.5 * alphaL * alphaL;
    double b = alphaL * t;
    v = std::exp(a + b);
  }
  else if (t <= alphaR)
  {
    v = std::exp(-0.5 * t * t);
  }
  else
  {
    double a = 0.5 * alphaR * alphaR;
    double b = alphaR * (-t);
    v = std::exp(a + b);
  }

  double tmin = (xmin - x0) / (xmin < x0 ? sigmaL : sigmaR);
  double tmax = (xmax - x0) / (xmax < x0 ? sigmaL : sigmaR);
  double sum = 0;
  if (tmin < -alphaL)
  {
    double a = 0.5 * alphaL * alphaL;
    double lv = tmin;
    double uv = std::min(tmax, -alphaL);
    sum += (sigmaL / alphaL) * (std::exp(a + alphaL * uv) - std::exp(a + alphaL * lv));
  }
  if (tmax > alphaR)
  {
    double a = 0.5 * alphaR * alphaR;
    double lv = std::max(tmin, alphaR);
    double uv = tmax;
    sum += (sigmaR / alphaR) * (std::exp(a - alphaR * lv) - std::exp(a - alphaR * uv));
  }
  if (tmin < alphaR && tmax > -alphaL)
  {
    double sigmaMin = (tmin < double(0)) ? sigmaL : sigmaR;
    double sigmaMax = (tmax < double(0)) ? sigmaL : sigmaR;
    sum += sqrtPiOver2 * (sigmaMax * std::erf(std::min(tmax, alphaR) / sqrt2) - sigmaMin * std::erf(std::max(tmin, -alphaL) / sqrt2));
  }

  return (v / sum);
}

/**
 * @brief 自定义的 RooFit 概率密度函数 (PDF)，表示带有指数尾的高斯分布。
 * 
 * 该类继承自 RooAbsPdf，用于在 RooFit 中定义和使用带有指数尾的高斯分布。
 */
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


/**
 * @brief 使用 RooFit 对信号、背景和观测数据进行联合拟合。
 * 
 * @param hsig 指向信号直方图的指针。
 * @param hbkg 指向背景直方图的指针。
 * @param hneg 指向观测数据直方图的指针。
 * @return std::vector<double> 包含拟合参数和误差的向量。
 *         - 前 11 个元素为拟合参数。
 *         - 后 11 个元素为对应的误差。
 **/
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

    minimizer.Minimize();

    const double* xs = minimizer.X();
    const double* errors = minimizer.Errors();

    std::vector<double> params(xs, xs + 11);
    params.insert(params.end(), errors, errors + 11);
    return params;
}


/**
 * @brief 使用迭代加权最小二乘法 (IWLS) 对信号、背景和观测数据进行拟合。
 * 
 * @param hsig 指向信号直方图的指针。
 * @param hbkg 指向背景直方图的指针。
 * @param hneg 指向观测数据直方图的指针。
 * @param nsig 引用变量，用于存储拟合得到的信号事件数量。
 * @param nsig_error 引用变量，用于存储拟合得到的信号事件数量的误差。
 * @param iteration 拟合过程的迭代次数。
 * 
 **/
void IWLS_FunctionFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg,double &nsig, double &nsig_error,Int_t iteration)
{
  TF1 *func = new TF1("func", HistFitFunc, -1, 3, 2);
  func->SetNpx(400);
 
  TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone")+TString(hbkg->GetTitle()));
  TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone")+TString(hsig->GetTitle()));
  TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone")+TString(hneg->GetTitle()));

  hpi_clone ->Scale(1/hpi_clone->GetSumOfWeights()/hpi_clone->GetBinWidth(1));
  hpr_clone ->Scale(1/hpr_clone->GetSumOfWeights()/hpr_clone->GetBinWidth(1));
  hneg_clone->Scale(1/hneg_clone->GetSumOfWeights()/hneg_clone->GetBinWidth(1));

	double hneg_error[hneg_clone->GetNbinsX()];
	for (int i = 0; i < hneg_clone->GetNbinsX(); ++i) {
			hneg_error[i] = 1.0;
	}
   	
   double hpi_error[hpi_clone->GetNbinsX()];
	for (int i = 0; i < hpi_clone->GetNbinsX(); ++i) {
			hpi_error[i] = 1.0;
	}


    TF1* fit1=new TF1("fpr",funcExpGausExp,-1,3,7);
    TF1* fit2=new TF1("fpi",funcExpGausExp,-1,3,7);
    TF1* fit3=new TF1("fneg",combine_chi2,-1,3,15);

    setParameterWithDebug(fit1, 0, hpr_clone->GetMean(), 0.5, 1.1);
    setParameterWithDebug(fit1, 1, hpr_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit1, 2, 1, 0.01, 5.0);
    setParameterWithDebug(fit1, 3, hpr_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit1, 4, 1, 0.01, 6.0);
    fit1->FixParameter(5, -1);
    fit1->FixParameter(6, 3);

    setParameterWithDebug(fit2, 0, 0.001, 0.001, 0.35);
    setParameterWithDebug(fit2, 1, hpi_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit2, 2, 1, 0.01, 5.0);
    setParameterWithDebug(fit2, 3, hpi_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit2, 4, 1, 0.01, 5.0);
    fit2->FixParameter(5, -1);
    fit2->FixParameter(6, 3);

    hpr_clone->Fit(fit1,"q0", "R");


   for(int iter=0;iter<iteration;iter++){
   
      for(int i = 1; i <= hpi_clone->GetNbinsX(); ++i) {
      hpi_clone->SetBinError(i,hpi_error[i-1]);
      }
     //std::cout<<"pi iteration: "<<iter<<std::endl;
      hpi_clone->Fit(fit2, "q0", "R");
      
      for(int i = 1; i <= hpi_clone->GetNbinsX(); ++i) {
      hpi_error[i-1]=TMath::Sqrt(fit2->Eval(hpi_clone->GetBinCenter(i)));
      }
   }

    double comb_params[15];
    for (int i = 0; i < 7; ++i) {
        comb_params[i] = fit1->GetParameter(i);
    }
    for (int i = 0; i < 7; ++i) {
        comb_params[i + 7] = fit2->GetParameter(i);
    }
    comb_params[14] = 0.1;
    fit3->SetParameters(comb_params);
    for (int i = 0; i < 14; ++i) {
      fit3->FixParameter(i, comb_params[i]);
    }

    fit3->SetParLimits(14,0.,1.);
    fit3->SetNpx(150);


   for(int iter=0;iter<iteration;iter++){
      if(iter!=0){
         for(int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
         hneg_clone->SetBinError(i,hneg_error[i-1]);
         }
      }

      hneg_clone->Fit(fit3, "q0", "R");

      for(int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
      hneg_error[i-1]=TMath::Sqrt(fit3->Eval(hneg_clone->GetBinCenter(i)));
      }
   }


  nsig = fit3->GetParameter(14);
  nsig_error=fit3->GetParError(14);


  delete fit1;
  delete fit2;
  delete fit3;
}

class QMinimizer : public ROOT::Math::IBaseFunctionMultiDim {
public:
    QMinimizer(TH1* signal, TH1* background, TH1* observed)
        : sig(dynamic_cast<TH1*>(signal->Clone("signal_clone"))),
          bkg(dynamic_cast<TH1*>(background->Clone("background_clone"))),
          data(observed),
          y1_opt(0), y2_opt(0), y1_err(0), y2_err(0) {}

    ~QMinimizer() override {
        delete sig;
        delete bkg;
    }

    // 克隆函数
    QMinimizer* Clone() const override {
        return new QMinimizer(sig, bkg, data);
    }

    // 维度数量：2
    unsigned int NDim() const override {
        return 2;
    }

    // 核心评估函数
    double DoEval(const double* p) const override {
        double y1 = p[0];
        double y2 = p[1];
        return calculateQ(y1, y2);
    }

    // 优化接口：运行优化并返回结果
    void Minimize() {
        ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimizer->SetFunction(*this);
        minimizer->SetMaxFunctionCalls(1000000);
        minimizer->SetMaxIterations(100000);
        minimizer->SetTolerance(0.001);

        // 设置参数初始值及边界
        minimizer->SetLimitedVariable(0, "y1", 150, 0.01, 0.0, 300.); // 参数 y1
        minimizer->SetLimitedVariable(1, "y2", 250, 0.01, 0.0, 300.); // 参数 y2

        // 执行最小化
        minimizer->Minimize();

        // 提取优化结果
        const double* results = minimizer->X();
        const double* errors = minimizer->Errors();

        y1_opt = results[0];
        y2_opt = results[1];
        y1_err = errors[0];
        y2_err = errors[1];

        delete minimizer;
    }

    // 获取结果
    double GetY1() const { return y1_opt/400.; }
    double GetY2() const { return y2_opt/400.; }
    double GetY1Error() const { return y1_err/400.; }
    double GetY2Error() const { return y2_err/400.; }

private:
    TH1* sig;
    TH1* bkg;
    TH1* data;
    double y1_opt, y2_opt; // 最优参数
    double y1_err, y2_err; // 参数误差

    // Q 计算逻辑
    double calculateQ(double y1, double y2) const {
           double Q = 0;
        for (int i = 0; i < sig->GetNbinsX(); i++) {
            double n = data->GetBinContent(i + 1);
            double Vn = data->GetBinContent(i + 1);
            double t = 1;
            double mu0 = y1 * sig->GetBinContent(i + 1) / sig->GetSumOfWeights() +
                         y2 * bkg->GetBinContent(i + 1) / bkg->GetSumOfWeights();

            if(mu0!=0){

                double Vmu = y1 * y1 * sig->GetBinContent(i + 1) / (sig->GetSumOfWeights()*sig->GetSumOfWeights()) +
                            y2 * y2 * bkg->GetBinContent(i + 1) / (bkg->GetSumOfWeights()*bkg->GetSumOfWeights());
                double s = (Vmu != 0) ? mu0 / Vmu : 0;
                double beta=(t * n + s * mu0) / (t * mu0 + s * mu0);
                double cash1 = cash(t * n, beta * t * mu0);
                double cash2 = cash(s * mu0, beta * s * mu0);
                Q += cash1 + cash2;
            }
            else{
                Q+=0;
            }
        } 
        return Q;
    }

    // Cash 统计量计算
    double cash(double observed, double expected) const {
        if (expected <= 0) return 0;
        if (observed <= 0) return 0;
        return expected - observed + observed * log(observed / expected);
    }
};


/** paper method
 * @brief 对信号、背景和观测数据进行 Q 优化。
 * 
 * @param signal 指向信号直方图的指针。
 * @param background 指向背景直方图的指针。
 * @param data 指向观测数据直方图的指针。
 * @return std::vector<double> 包含优化结果和误差的向量。
 *         - par[0]: 优化得到的 y1 值（signal frac）。
 *         - par[1]: 优化得到的 y2 值（bkg    frac）。
 *         - par[2]: y1 的误差。
 *         - par[3]: y2 的误差。
 **/
std::vector <double> RunQOptimization(TH1* signal, TH1* background, TH1* data) {
    QMinimizer qMinimizer(signal, background, data);
    qMinimizer.Minimize();
    std::vector<double> par(4);
    par[0] = qMinimizer.GetY1();
    par[1] = qMinimizer.GetY2();
    par[2] = qMinimizer.GetY1Error();
    par[3] = qMinimizer.GetY2Error();
    return par;
}

/**
 * @brief 使用 TFractionFitter 对信号、背景和观测数据进行拟合。
 * 
 * @param hsig 指向信号直方图的指针。
 * @param hbkg 指向背景直方图的指针。
 * @param hneg 指向观测数据直方图的指针。
 * @param nsig 引用变量，用于存储拟合得到的信号事件数量。
 * @param nsig_error 引用变量，用于存储拟合得到的信号事件数量的误差。
 * 
 **/
void TFracFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double &nsig,double &nsig_error)
{

   TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone")+TString(hbkg->GetTitle()));
   TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone")+TString(hsig->GetTitle()));
   TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone")+TString(hneg->GetTitle()));
   //   hpi_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.22))/hpi_clone->GetBinContent(hpi_clone->FindBin(0.22)));
   //   hpr_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.88))/hpr_clone->GetBinContent(hpr_clone->FindBin(0.88)));
   hpi_clone ->Scale(1/hpi_clone->GetSumOfWeights()/hpi_clone->GetBinWidth(1));
   hpr_clone ->Scale(1/hpr_clone->GetSumOfWeights()/hpr_clone->GetBinWidth(1));
   hneg_clone->Scale(1/hneg_clone->GetSumOfWeights()/hneg_clone->GetBinWidth(1));
   TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
   mc->Add(hpi_clone);
   mc->Add(hpr_clone);
   TFractionFitter* fit = new TFractionFitter(hneg_clone, mc); // initialise
   fit->Constrain(0, 0.0, 1.0);               // constrain fraction 0 to be between 0 and 1
   fit->Constrain(1, 0.0, 1.0);               // constrain fraction 1 to be between 0 and 1
   fit->SetRangeX(1, hneg->GetNbinsX());      // use all bins in the fit
   Int_t status = fit->Fit();                 // perform the fit
   std::cout << "fit status: " << status << std::endl;

   if (status == 0) {                         // check on fit status
      //TH1F* result = (TH1F*) fit->GetPlot();
      //hneg->Draw("Ep");
      //result->Draw("same");
      fit->GetResult(1,nsig,nsig_error);
      // nsig=nsig*hneg->GetEntries();
   }



   delete fit;
}