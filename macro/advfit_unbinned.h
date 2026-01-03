/* Usage:
bool MLfit = true;
bool fixsum = true;
TF_Fitter *myFit = new TF_Fitter(hsig, hbkg, hdata, MLfit, fixsum);
double nsig = myFit->nsig; 
double nbkg = myFit->nbkg; 
*/
#include <Math/IFunction.h>
#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include "TF1.h"
#include "TROOT.h"
#include <TTree.h>
#include "TFractionFitter.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "TError.h"
#include "TFractionFitter.h"
#include "TMinuit.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "Fit/Fitter.h"



TH1D *fHistPi, *fHistPr;

double FitGaussian(TH1D *h, double &m, double &s)
{
  TF1 *fgaus = new TF1("fgaus","gaus", -3, 3);
  fgaus->SetParameters(h->GetMaximum(), 1, 0.05);
  h->Fit(fgaus,"qmnl", "",0.8, 1.2);
  m = fgaus->GetParameter(1);
  s = fgaus->GetParameter(2);

  h->Fit(fgaus, "qmnl", "", m - 1 * s, m + 1 * s);
  m = fgaus->GetParameter(1);
  s = fgaus->GetParameter(2);
  h->Fit(fgaus, "qmnl", "", m - 1.5 * s, m + 1.5 * s);
  m = fgaus->GetParameter(1);
  s = fgaus->GetParameter(2);

  h->SetTitle(Form("Bias=%.1f%%, Width=%.1f%%, #sqrt{Bias^{2}+Width^{2}} = %.1f %%",
                       (m-1)*100,
                       s*100, 
                       TMath::Sqrt( (m-1)*(m-1)+s*s)*100)
                       );
  return TMath::Sqrt((m - 1) * (m - 1) + s * s) * 100;
}

void setParameterWithDebug(TF1* func, int parIndex, double value, double low, double up) {
    if (value < low || value > up) {
        std::cerr << func->GetName()<<" "<<"Warning: Parameter " << parIndex << " value " << value << " is outside the range [" << low << ", " << up << "]. Setting to (" << low << " + " << up << ") / 2." << std::endl;
        value = (low + up) / 2;
    }
    func->SetParameter(parIndex, value);
    func->SetParLimits(parIndex, low, up);
}

// Double_t HistFitFunc(double *x, double *p)
// {
//   Double_t xx = x[0];
//   Double_t p0 = p[0];
//   Double_t p1 = p[1];
//   Int_t ib = fHistPi ? fHistPi->FindBin(xx) : 0;
//   Double_t y0 = fHistPi?fHistPi->GetBinContent(ib):0;
//   Double_t y1 = fHistPr?fHistPr->GetBinContent(ib):0;
//   return p0*y0+p1*y1;
// }
// double fHistFit_nbkg = 0;
// void HistFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double &nsig)
// {
//   TF1 *func = new TF1("func", HistFitFunc, -2.55, 3.55, 2);
//   func->SetNpx(400);
 
//   TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone")+TString(hbkg->GetTitle()));
//   TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone")+TString(hsig->GetTitle()));
//   hpi_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.22))/hpi_clone->GetBinContent(hpi_clone->FindBin(0.22)));
//   hpr_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.88))/hpr_clone->GetBinContent(hpr_clone->FindBin(0.88)));
//   double p[2] = {1,1};
//   func->SetParameters(p);
//   func->SetParLimits(0, 0, 1e3);
//   func->SetParLimits(1, 0, 1e3);

//   fHistPi = hpi_clone;
//   fHistPr = hpr_clone;

//   hneg->Fit(func,"qL0", "");
//   func->GetParameters(p);
//   nsig = hpr_clone->Integral()*p[1];
//   fHistFit_nbkg = hpi_clone->Integral()*p[0];

//   delete func;
//   // delete hpi_clone;
//   // delete hpr_clone;

// }


Double_t HistFitFunc(double *x, double *p)
{
  Double_t xx = x[0];
  Double_t p0 = p[0];
  Double_t p1 = 1-p[0];
  Int_t ib = fHistPi ? fHistPi->FindBin(xx) : 0;
  Double_t y0 = fHistPi?fHistPi->GetBinContent(ib):0;
  Double_t y1 = fHistPr?fHistPr->GetBinContent(ib):0;
  return p0*y0+p1*y1;
}
double fHistFit_nbkg = 0;
void CHistFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double &nsig)
{
  TF1 *func = new TF1("func", HistFitFunc, -1, 3, 2);
  func->SetNpx(400);
 
  TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone")+TString(hbkg->GetTitle()));
  TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone")+TString(hsig->GetTitle()));
  TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone")+TString(hneg->GetTitle()));
//   hpi_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.22))/hpi_clone->GetBinContent(hpi_clone->FindBin(0.22)));
//   hpr_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.88))/hpr_clone->GetBinContent(hpr_clone->FindBin(0.88)));
  hpi_clone->Scale(1/hpi_clone->GetSumOfWeights()/hpi_clone->GetBinWidth(1));
  hpr_clone->Scale(1/hpr_clone->GetSumOfWeights()/hpr_clone->GetBinWidth(1));
  hneg_clone->Scale(1/hneg_clone->GetSumOfWeights()/hneg_clone->GetBinWidth(1));
  double p[1] = {0.5};
  func->SetParameters(p);
  func->SetParLimits(0, 0, 1);


  fHistPi = hpi_clone;
  fHistPr = hpr_clone;

  hneg_clone->Fit(func,"q0", "");
  func->GetParameters(p);
  nsig = hneg->GetEntries()*(1-p[0]);
  //fHistFit_nbkg = hneg->GetEntries()*p[0];

  delete func;
  // delete hpi_clone;
  // delete hpr_clone;

}

void LHistFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double &nsig)
{
  TF1 *func = new TF1("func", HistFitFunc, -1, 3, 2);
  func->SetNpx(400);
 
  TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone")+TString(hbkg->GetTitle()));
  TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone")+TString(hsig->GetTitle()));
  TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone")+TString(hneg->GetTitle()));
//   hpi_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.22))/hpi_clone->GetBinContent(hpi_clone->FindBin(0.22)));
//   hpr_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.88))/hpr_clone->GetBinContent(hpr_clone->FindBin(0.88)));
  hpi_clone->Scale(1/hpi_clone->GetSumOfWeights()/hpi_clone->GetBinWidth(1));
  hpr_clone->Scale(1/hpr_clone->GetSumOfWeights()/hpr_clone->GetBinWidth(1));
  hneg_clone->Scale(1/hneg_clone->GetSumOfWeights()/hneg_clone->GetBinWidth(1));
  double p[1] = {0.5};
  func->SetParameters(p);
  func->SetParLimits(0, 0, 1);


  fHistPi = hpi_clone;
  fHistPr = hpr_clone;

  hneg_clone->Fit(func,"qL0", "");
  func->GetParameters(p);
  nsig = hneg->GetEntries()*(1-p[0]);
  //fHistFit_nbkg = hneg->GetEntries()*p[0];

  delete func;
  // delete hpi_clone;
  // delete hpr_clone;

}

// void IWLS_HistFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double * hneg_error,double &nsig)
// {
//   TF1 *func = new TF1("func", HistFitFunc, -1, 3, 2);
//   func->SetNpx(400);
 
//   TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone")+TString(hbkg->GetTitle()));
//   TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone")+TString(hsig->GetTitle()));
//   TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone")+TString(hneg->GetTitle()));
// //   hpi_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.22))/hpi_clone->GetBinContent(hpi_clone->FindBin(0.22)));
// //   hpr_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.88))/hpr_clone->GetBinContent(hpr_clone->FindBin(0.88)));
//   hpi_clone->Scale(1/hpi_clone->GetSumOfWeights()/hpi_clone->GetBinWidth(1));
//   hpr_clone->Scale(1/hpr_clone->GetSumOfWeights()/hpr_clone->GetBinWidth(1));
//   hneg_clone->Scale(1/hneg_clone->GetSumOfWeights()/hneg_clone->GetBinWidth(1));

//    // for (int i = 1; i <= hpi_clone->GetNbinsX(); ++i) {
//    //    hpi_clone->SetBinError(i, hbkg_error[i-1]);
//    // }

//    // for (int i = 1; i <= hpr_clone->GetNbinsX(); ++i) {
//    //    hpr_clone->SetBinError(i, hsig_error[i-1]);
//    // }

//    for (int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
//       hneg_clone->SetBinError(i, hneg_error[i-1]);
//    }
//   double p[1] = {0.5};
//   func->SetParameters(p);
//   func->SetParLimits(0, 0, 1);


//   fHistPi = hpi_clone;
//   fHistPr = hpr_clone;

//   hneg_clone->Fit(func,"q0", "");
  
//   for(int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
//    hneg_error[i-1]=TMath::Sqrt(func->Eval(hneg_clone->GetBinCenter(i)));
//   }

//   func->GetParameters(p);
//   nsig = hneg->GetEntries()*(1-p[0]);
//   //fHistFit_nbkg = hneg->GetEntries()*p[0];

//   delete func;
//   // delete hpi_clone;
//   // delete hpr_clone;

// }

void IWLS_HistFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg,double &nsig,double &nsig_error, Int_t iteration)
{
  TF1 *func = new TF1("func", HistFitFunc, -1, 3, 2);
  func->SetNpx(400);
  int Ntot= hneg->GetEntries();
  TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone")+TString(hbkg->GetTitle()));
  TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone")+TString(hsig->GetTitle()));
  TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone")+TString(hneg->GetTitle()));
//   hpi_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.22))/hpi_clone->GetBinContent(hpi_clone->FindBin(0.22)));
//   hpr_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.88))/hpr_clone->GetBinContent(hpr_clone->FindBin(0.88)));
  hpi_clone ->Scale(1/hpi_clone->GetSumOfWeights());
  hpr_clone ->Scale(1/hpr_clone->GetSumOfWeights());
  hneg_clone->Scale(1/hneg_clone->GetSumOfWeights());

	double hneg_error[hneg_clone->GetNbinsX()];
	for (int i = 0; i < hneg_clone->GetNbinsX(); ++i) {
			hneg_error[i] = 1.0;
	}


  double p[1] = {0.4};
  func->SetParLimits(0, 0, 1);

  double m=0;
  fHistPi = hpi_clone;
  fHistPr = hpr_clone;
//error iteration
	for(int iter=0;iter<iteration;iter++){
		
		for(int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
		hneg_clone->SetBinError(i,hneg_error[i-1]);
		}

		hneg_clone->Fit(func,"q0", "");
		
		for(int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
         m=func->GetParameter(0)*hpi_clone->GetBinContent(i)+
            (1-func->GetParameter(0))*hpr_clone->GetBinContent(i);
		   hneg_error[i-1]=TMath::Sqrt(m*(1-m)/Ntot);
		}
	}
  func->GetParameters(p);
  nsig = hneg->GetEntries()*(1-p[0]);
   nsig_error = Ntot*func->GetParError(0);
  //fHistFit_nbkg = hneg->GetEntries()*p[0];

  delete func;
  // delete hpi_clone;
  // delete hpr_clone;

}

TF1 *fpi, *fpr;
Double_t FuncFitFunc(double *x,double *p){
    Double_t xx = x[0];
    Double_t p0 = p[0];
    Double_t p1 = p[1];
    Double_t y0 = fpi?fpi->Eval(xx):0;
    Double_t y1 = fpr?fpr->Eval(xx):0;
    return p0*y0+p1*y1;
}
void FuncFit(TF1 fsig, TF1 fbkg, TH1D hneg, double &nsig){
    TF1 *func = new TF1("func", FuncFitFunc, -2.55, 3.55, 2);
    func->SetNpx(400);
    fpi = &fbkg;
    fpr = &fsig;


    double p[2] = {
      hneg.GetBinContent(hneg.FindBin(0.22))/fbkg.GetParameter(0),
      hneg.GetBinContent(hneg.FindBin(0.88))/fsig.GetParameter(0)};

    func->SetParameters(p);
    func->SetParLimits(0, 0, 10);
    func->SetParLimits(1, 0, 10);

    hneg.Fit(func,"qLN0", "");
    func->GetParameters(p);
    TH1D *hpr = (TH1D*)hneg.Clone("hneg");
    hpr->Reset();
    for(int j=1;j<=hpr->GetNbinsX();j++)
      hpr->SetBinContent(j, fsig.Eval(hpr->GetBinCenter(j)));
    nsig = hpr->Integral()*p[1];
    delete hpr;
    delete func;
}

double calculate_norm(const double *const pars, const double &xmin, const double &xmax) 
{
   double sqrtPiOver2 = double(1.2533141373);
   double sqrt2 = double(1.4142135624);

   const double &x0 = pars[0];
   const double &sigmaL = pars[1];
   const double &alphaL = pars[2];
   const double &nL = pars[3];
   const double &sigmaR = pars[4];
   const double &alphaR = pars[5];
   const double &nR = pars[6];

   double tmin = (xmin - x0) / (xmin < x0 ? sigmaL : sigmaR);
   double tmax = (xmax - x0) / (xmax < x0 ? sigmaL : sigmaR);

   double sum = double(0);
   if (tmin < -alphaL)
   {
      double a = TMath::Power(nL / alphaL, nL) * TMath::Exp(-0.5 * alphaL * alphaL);
      double b = nL / alphaL - alphaL;
      double lv = tmin;
      double uv = std::min(tmax, -alphaL);
      sum += a * sigmaL / (1.0 - nL) * (1.0 / (TMath::Power(b - lv, nL - 1.0)) - 1.0 / (TMath::Power(b - uv, nL - 1.0)));
   }
   if (tmax > alphaR)
   {
      double a = TMath::Power(nR / alphaR, nR) * TMath::Exp(-0.5 * alphaR * alphaR);
      double b = nR / alphaR - alphaR;
      double lv = std::max(tmin, alphaR);
      double uv = tmax;
      sum += a * sigmaR / (1.0 - nR) * (1.0 / (TMath::Power(b + uv, nR - 1.0)) - 1.0 / (TMath::Power(b + lv, nR - 1.0)));
   }
   if (tmin < alphaR && tmax > -alphaL)
   {
      double sigmaMin = (tmin < double(0)) ? sigmaL : sigmaR;
      double sigmaMax = (tmax < double(0)) ? sigmaL : sigmaR;
      sum += sqrtPiOver2 * (sigmaMax * TMath::Erf(std::min(tmax, alphaR) / sqrt2) - sigmaMin * TMath::Erf(std::max(tmin, -alphaL) / sqrt2));
   }
   return sum;
}

double CrystalBall(double* x, double* p) {
  double xx = x[0];
  double x0 = p[0];
  double sigmaL = std::abs(p[1]);
  double alphaL = std::abs(p[2]);
  double nL     = std::abs(p[3]);
  double sigmaR = std::abs(p[4]);
  double alphaR = std::abs(p[5]);
  double nR     = std::abs(p[6]);
  double norm   = std::abs(p[7]);
  double xmin   = p[8];
  double xmax   = p[9];

  double t = (xx - x0) / (xx < x0 ? sigmaL : sigmaR);

   //FIXME: Try to Normalize it so that n0 = total events, Integral(v) = 1
  double nrm = calculate_norm(p, xmin, xmax);
  double v = 0.0;
  if (t < -alphaL) {
    double a = std::pow(nL / alphaL, nL) * std::exp(-0.5 * alphaL * alphaL);
    double b = nL / alphaL - alphaL;
    v = a / std::pow(b - t, nL);
  } else if (t <= alphaR) {
    v = std::exp(-0.5 * t * t); 
  } else {
    double a = std::pow(nR / alphaR, nR) * std::exp(-0.5 * alphaR * alphaR);
    double b = nR / alphaR - alphaR;
    v = a / std::pow(b + t, nR);
  }   
//   return  v/nrm;
   return  norm*v;
}
TF1 CB_Fit(TH1D* hist, Double_t *p){
   gROOT->cd();
   TF1 cb(Form("cb")+TString(hist->GetTitle()), CrystalBall, p[7], p[8], 9);
   // cb.SetParameters(p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
   // cb.SetParLimits(0, p[0] - p[1]*3, p[0] + p[1]*3);                        // x0
   // cb.SetParLimits(1, p[1]*0.2, p[1]*5);                                     // sigL
   // cb.SetParLimits(2, 0.1, 10);                                       // alpha L
   // cb.SetParLimits(3, 1.0001, 100);                                         // nL
   // cb.SetParLimits(4, p[1]*0.2, p[1]*5);                                     // sig R
   // cb.SetParLimits(5, 0.1, 10);                                       // alpha R
   // cb.SetParLimits(6, 1.0001, 100);                                         // nR
   cb.SetParameters(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]);
   cb.SetNpx(150);
   TH1D hfit = *(TH1D*)hist->Clone( TString(hist->GetTitle())+"_CBFit");
   hfit.Scale(1/hfit.GetSumOfWeights()/hfit.GetBinWidth(1));
   hfit.Fit(&cb, "LMNq", "");
   return cb;
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
  //const double &norm = p[5];
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

double combine_chi2(double *x, double *p) {
    TF1 *fFit1 = new TF1("fFpr", funcExpGausExp, -1, 3, 7);
    TF1 *fFit2 = new TF1("fFpi", funcExpGausExp, -1, 3, 7);
   //  fFit1->SetParameters(p);
   //  fFit2->SetParameters(&p[7]);
    double fsig = fFit1->EvalPar(x, p);
    double fbkg = fFit2->EvalPar(x, &p[7]);
    return p[14]*fsig + (1-p[14])*fbkg;
}


void IWLS_FunctionFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg,double &nsig, double &nsig_error,Int_t iteration)
{
  TF1 *func = new TF1("func", HistFitFunc, -1, 3, 2);
  func->SetNpx(400);
 
  TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone")+TString(hbkg->GetTitle()));
  TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone")+TString(hsig->GetTitle()));
  TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone")+TString(hneg->GetTitle()));
//   hpi_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.22))/hpi_clone->GetBinContent(hpi_clone->FindBin(0.22)));
//   hpr_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.88))/hpr_clone->GetBinContent(hpr_clone->FindBin(0.88)));
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
   //std::cout<<"---pi done---"<<std::endl;
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
    //fit3->SetParError(14, 0.01);



   for(int iter=0;iter<iteration;iter++){
      if(iter!=0){
         for(int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
         hneg_clone->SetBinError(i,hneg_error[i-1]);
         }
      }
      //std::cout<<"---fit3---"<<std::endl;
      hneg_clone->Fit(fit3, "q0", "R");
      //std::cout<<"---fit3 d---"<<std::endl;
      //std::cout<<fit3->GetParameter(14)<<std::endl;
      for(int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
      hneg_error[i-1]=TMath::Sqrt(fit3->Eval(hneg_clone->GetBinCenter(i)));
      }
   }


  nsig = fit3->GetParameter(14);
  nsig_error=fit3->GetParError(14);
  //fHistFit_nbkg = hneg->GetEntries()*p[0];

  delete fit1;
  delete fit2;
  delete fit3;
  // delete hpi_clone;
  // delete hpr_clone;

}

void IWLS_FunctionFit_Minuit2(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double &nsig, double &nsig_error, Int_t iteration) {
    // 创建 TF1 对象
    TF1* fit1 = new TF1("fit1", funcExpGausExp, -1, 3, 7);
    TF1* fit2 = new TF1("fit2", funcExpGausExp, -1, 3, 7);
    TF1* fit3 = new TF1("fit3", combine_chi2, -1, 3, 15);

    // 设置初始参数
    fit1->SetParameters(0.8, 1.0, 1.0, 1.0, 1.0, -1.0, 3.0);
    fit2->SetParameters(0.001, 1.0, 1.0, 1.0, 1.0, -1.0, 3.0);

    // 克隆直方图
    TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone") + TString(hbkg->GetTitle()));
    TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone") + TString(hsig->GetTitle()));
    TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone") + TString(hneg->GetTitle()));

    hpi_clone->Scale(1 / hpi_clone->GetSumOfWeights() / hpi_clone->GetBinWidth(1));
    hpr_clone->Scale(1 / hpr_clone->GetSumOfWeights() / hpr_clone->GetBinWidth(1));
    hneg_clone->Scale(1 / hneg_clone->GetSumOfWeights() / hneg_clone->GetBinWidth(1));

    double hneg_error[hneg_clone->GetNbinsX()];
    for (int i = 0; i < hneg_clone->GetNbinsX(); ++i) {
        hneg_error[i] = 1.0;
    }

    double hpi_error[hpi_clone->GetNbinsX()];
    for (int i = 0; i < hpi_clone->GetNbinsX(); ++i) {
        hpi_error[i] = 1.0;
    }

    // 直接拟合 hpr_clone
    hpr_clone->Fit(fit1, "q0", "R");

    // 定义目标函数2
    auto chi2_2 = [&](const double* params) {
        fit2->SetParameters(params);
        double chi2 = 0.0;
      for (int i = 1; i <= hpi_clone->GetNbinsX(); ++i) {
          double x = hpi_clone->GetBinCenter(i);
          double y = hpi_clone->GetBinContent(i);
          double y_fit = fit2->Eval(x);
          double error = hpi_error[i - 1];
          if (error != 0) {
             chi2 += std::pow((y - y_fit) / error, 2);
          } else {
             chi2 += 0;
          }
      }

        return chi2;
    };

    // 创建 Minuit2 最小化器
    ROOT::Minuit2::Minuit2Minimizer minimizer2(ROOT::Minuit2::kMigrad);
    minimizer2.SetMaxFunctionCalls(1000000);
    minimizer2.SetMaxIterations(100000);
    minimizer2.SetTolerance(0.001);

    // 设置目标函数
    ROOT::Math::Functor f2(chi2_2, 7);
    minimizer2.SetFunction(f2);

   // 设置初始参数和范围
   minimizer2.SetVariable(0, "param0", 0.001, 0.01);
   minimizer2.SetVariableLimits(0, 0.001, 0.35);
   minimizer2.SetVariable(1, "param1", hpi_clone->GetRMS(), 0.01);
   minimizer2.SetVariableLimits(1, 0.05, 1.1);
   minimizer2.SetVariable(2, "param2", 1, 0.01);
   minimizer2.SetVariableLimits(2, 0.01, 5.0);
   minimizer2.SetVariable(3, "param3", hpi_clone->GetRMS(), 0.01);
   minimizer2.SetVariableLimits(3, 0.05, 1.1);
   minimizer2.SetVariable(4, "param4", 1, 0.01);
   minimizer2.SetVariableLimits(4, 0.01, 5.0);
   minimizer2.SetFixedVariable(5, "param5", -1);
   minimizer2.SetFixedVariable(6, "param6", 3);

    // 执行迭代最小化
    for (int iter = 0; iter < iteration; ++iter) {
        if (iter != 0) {
            for (int i = 1; i <= hpi_clone->GetNbinsX(); ++i) {
                hpi_clone->SetBinError(i, hpi_error[i - 1]);
            }
        }

        minimizer2.Minimize();

        for (int i = 1; i <= hpi_clone->GetNbinsX(); ++i) {
            hpi_error[i - 1] = TMath::Sqrt(fit2->Eval(hpi_clone->GetBinCenter(i)));
        }
    }

    // 获取拟合结果
    const double* xs2 = minimizer2.X();
    const double* errors2 = minimizer2.Errors();

    // 打印拟合结果
    for (int i = 0; i < 7; ++i) {
        std::cout << "Parameter " << i << " (fit2): " << xs2[i] << " ± " << errors2[i] << std::endl;
    }

    // 更新组合参数
    double comb_params[15];
    for (int i = 0; i < 7; ++i) {
        comb_params[i] = fit1->GetParameter(i);
        comb_params[i + 7] = xs2[i];
    }
    comb_params[14] = 0.1;
   //  fit3->SetParameters(comb_params);
   //  for (int i = 0; i < 14; ++i) {
   //      fit3->FixParameter(i, comb_params[i]);
   //  }

   //  fit3->SetParLimits(14, 0., 1.);
   //  fit3->SetNpx(150);

    // 定义目标函数3
    auto chi2_3 = [&](const double* params) {
        fit3->SetParameters(params);
        double chi2 = 0.0;
        for (int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
            double x = hneg_clone->GetBinCenter(i);
            double y = hneg_clone->GetBinContent(i);
            double y_fit = fit3->Eval(x);
            double error = hneg_error[i - 1];
            if (error != 0) {
               chi2 += std::pow((y - y_fit) / error, 2);
            } else {
               chi2 += 0;
            }
        }
        return chi2;
    };

    // 创建 Minuit2 最小化器
    ROOT::Minuit2::Minuit2Minimizer minimizer3(ROOT::Minuit2::kMigrad);
    minimizer3.SetMaxFunctionCalls(1000000);
    minimizer3.SetMaxIterations(100000);
    minimizer3.SetTolerance(0.001);

    // 设置目标函数
    ROOT::Math::Functor f3(chi2_3, 15);
    minimizer3.SetFunction(f3);

    // 设置初始参数和范围
   minimizer3.SetFixedVariable(0, "param0", comb_params[0]);
   minimizer3.SetFixedVariable(1, "param1", comb_params[1]);
   minimizer3.SetFixedVariable(2, "param2", comb_params[2]);
   minimizer3.SetFixedVariable(3, "param3", comb_params[3]);
   minimizer3.SetFixedVariable(4, "param4", comb_params[4]);
   minimizer3.SetFixedVariable(5, "param5", -1);
   minimizer3.SetFixedVariable(6, "param6", 3);

   minimizer3.SetFixedVariable(7, "param7", comb_params[7]);
   minimizer3.SetFixedVariable(8, "param8", comb_params[8]);
   minimizer3.SetFixedVariable(9, "param9", comb_params[9]);
   minimizer3.SetFixedVariable(10, "param10", comb_params[10]);
   minimizer3.SetFixedVariable(11, "param11", comb_params[11]);
   minimizer3.SetFixedVariable(12, "param12", -1);
   minimizer3.SetFixedVariable(13, "param13", 3);

   minimizer3.SetVariable(14, "param14", comb_params[14], 0.01);
   minimizer3.SetVariableLimits(14, 0.0, 1.0);

    // 执行迭代最小化
    for (int iter = 0; iter < iteration; ++iter) {
        if (iter != 0) {
            for (int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
                hneg_clone->SetBinError(i, hneg_error[i - 1]);
            }
        }

        minimizer3.Minimize();

        for (int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
            hneg_error[i - 1] = TMath::Sqrt(fit3->Eval(hneg_clone->GetBinCenter(i)));
        }
    }

    // 获取拟合结果
    const double* xs3 = minimizer3.X();
    const double* errors3 = minimizer3.Errors();

    nsig = xs3[14];
   if (errors3[14] == 0) {
      nsig_error = -1;
   } else {
      nsig_error = errors3[14];
   }

    delete fit1;
    delete fit2;
    delete fit3;
}


double combine_crys(double *x, double *p) {
    TF1 *fFit1=new TF1("fit1",CrystalBall,-1,3,10);
    TF1 *fFit2=new TF1("fit2",CrystalBall,-1,3,10);
   //  fFit1->SetParameters(p);
   //  fFit2->SetParameters(&p[7]);
    double fsig = fFit1->EvalPar(x, p);
    double fbkg = fFit2->EvalPar(x, &p[10]);
    return p[20]*fsig + p[21]*fbkg;
}

void IWLS_CrysFunctionFit_Minuit2(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double &nsig, double &nsig_error, Int_t iteration) {
    // 创建 TF1 对象
    TF1 *fit1=new TF1("fit1",CrystalBall,-1,3,10);
    TF1 *fit2=new TF1("fit2",CrystalBall,-1,3,10);
    TF1* fit3=new TF1("fneg",combine_crys,-1,3,22);




    // 克隆直方图
    TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone") + TString(hbkg->GetTitle()));
    TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone") + TString(hsig->GetTitle()));
    TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone") + TString(hneg->GetTitle()));

   //  hpi_clone->Scale(1 / hpi_clone->GetSumOfWeights() / hpi_clone->GetBinWidth(1));
   //  hpr_clone->Scale(1 / hpr_clone->GetSumOfWeights() / hpr_clone->GetBinWidth(1));
   //  hneg_clone->Scale(1 / hneg_clone->GetSumOfWeights() / hneg_clone->GetBinWidth(1));

   setParameterWithDebug(fit1, 0, hpr_clone->GetMean(), 0.5, 1.3);
   setParameterWithDebug(fit1, 1, hpr_clone->GetRMS(), 0.05, 1.1);
   setParameterWithDebug(fit1, 2, 1, 0.01, 5.0);
   setParameterWithDebug(fit1, 3, 2,1.0001, 100);
   setParameterWithDebug(fit1, 4, hpr_clone->GetRMS(), 0.05, 1.1);
   setParameterWithDebug(fit1, 5, 1, 0.01, 6.0);
   setParameterWithDebug(fit1, 6, 2,1.0001, 100);
   setParameterWithDebug(fit1, 7, 100, 0, 1000000);//norm
   setParameterWithDebug(fit1, 8, -1, 0, -1); // xmin
   setParameterWithDebug(fit1, 9, 3, 3, 3);   // xmax

   setParameterWithDebug(fit2, 0, 0.001, 0.001, 0.35);
   setParameterWithDebug(fit2, 1, hpi_clone->GetRMS(), 0.05, 1.1);
   setParameterWithDebug(fit2, 2, 1, 0.01, 5.0);
   setParameterWithDebug(fit2, 3, 2,1.0001, 100);
   setParameterWithDebug(fit2, 4, hpi_clone->GetRMS(), 0.05, 1.1);
   setParameterWithDebug(fit2, 5, 1, 0.01, 6.0);
   setParameterWithDebug(fit2, 6, 2,1.0001, 100);
   setParameterWithDebug(fit2, 7, 100, 0, 1000000); //norm
   setParameterWithDebug(fit2, 8, -1, -1, -1); // xmin
   setParameterWithDebug(fit2, 9, 3, 3, 3);   // xmax

    double hneg_error[hneg_clone->GetNbinsX()];
    for (int i = 0; i < hneg_clone->GetNbinsX(); ++i) {
        hneg_error[i] = 1.0;
    }

    double hpi_error[hpi_clone->GetNbinsX()];
    for (int i = 0; i < hpi_clone->GetNbinsX(); ++i) {
        hpi_error[i] = 1.0;
    }

    // 直接拟合 hpr_clone
    hpr_clone->Fit(fit1, "q0", "R");

    // 定义目标函数2
    auto chi2_2 = [&](const double* params) {
        fit2->SetParameters(params);
        double chi2 = 0.0;
      for (int i = 1; i <= hpi_clone->GetNbinsX(); ++i) {
          double x = hpi_clone->GetBinCenter(i);
          double y = hpi_clone->GetBinContent(i);
          double y_fit = fit2->Eval(x);
          double error = hpi_error[i - 1];
          if (error != 0) {
             chi2 += std::pow((y - y_fit) / error, 2);
          } else {
             chi2 += 0;
          }
      }

        return chi2;
    };

    // 创建 Minuit2 最小化器
    ROOT::Minuit2::Minuit2Minimizer minimizer2(ROOT::Minuit2::kMigrad);
    minimizer2.SetMaxFunctionCalls(100000);
    minimizer2.SetMaxIterations(1000);
    minimizer2.SetTolerance(0.01);

    // 设置目标函数
    ROOT::Math::Functor f2(chi2_2, 10);
    minimizer2.SetFunction(f2);

   // 设置初始参数和范围
   minimizer2.SetVariable(0, "param0", 0.001, 0.01);
   minimizer2.SetVariableLimits(0, 0.001, 0.35);
   minimizer2.SetVariable(1, "param1", hpi_clone->GetRMS(), 0.01);
   minimizer2.SetVariableLimits(1, 0.05, 1.1);
   minimizer2.SetVariable(2, "param2", 1, 0.01);
   minimizer2.SetVariableLimits(2, 0.01, 5.0);
   minimizer2.SetVariable(3, "param3", 2, 0.01);
   minimizer2.SetVariableLimits(3, 1.0001, 100);
   minimizer2.SetVariable(4, "param4", hpi_clone->GetRMS(), 0.01);
   minimizer2.SetVariableLimits(4, 0.05, 1.1);
   minimizer2.SetVariable(5, "param5", 1, 0.01);
   minimizer2.SetVariableLimits(5, 0.01, 6.0);
   minimizer2.SetVariable(6, "param6", 2, 0.01);
   minimizer2.SetVariableLimits(6, 1.0001, 100);
   minimizer2.SetVariable(7, "param7", 10, 0.011);
   minimizer2.SetVariableLimits(7, 0, 1000000);
   minimizer2.SetFixedVariable(8, "param8", -1);
   minimizer2.SetFixedVariable(9, "param9", 3);

    // 执行迭代最小化
    for (int iter = 0; iter < iteration; ++iter) {
        if (iter != 0) {
            for (int i = 1; i <= hpi_clone->GetNbinsX(); ++i) {
                hpi_clone->SetBinError(i, hpi_error[i - 1]);
            }
        }

        minimizer2.Minimize();

        for (int i = 1; i <= hpi_clone->GetNbinsX(); ++i) {
            hpi_error[i - 1] = TMath::Sqrt(fit2->Eval(hpi_clone->GetBinCenter(i)));
        }
    }

    // 获取拟合结果
    const double* xs2 = minimizer2.X();
    const double* errors2 = minimizer2.Errors();

   //  // 打印拟合结果
   //  for (int i = 0; i < 10; ++i) {
   //      std::cout << "Parameter " << i << " (fit2): " << xs2[i] << " ± " << errors2[i] << std::endl;
   //  }

    // 更新组合参数
    double comb_params[22];
    for (int i = 0; i < 10; ++i) {
        comb_params[i] = fit1->GetParameter(i);
        comb_params[i + 10] = xs2[i];
    }
    comb_params[20] = 0.001;
    comb_params[21] = 1;
   //  fit3->SetParameters(comb_params);
   //  for (int i = 0; i < 14; ++i) {
   //      fit3->FixParameter(i, comb_params[i]);
   //  }

   //  fit3->SetParLimits(14, 0., 1.);
   //  fit3->SetNpx(150);

    // 定义目标函数3
    auto chi2_3 = [&](const double* params) {
        fit3->SetParameters(params);
        double chi2 = 0.0;
        for (int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
            double x = hneg_clone->GetBinCenter(i);
            double y = hneg_clone->GetBinContent(i);
            double y_fit = fit3->Eval(x);
            double error = hneg_error[i - 1];
            if (error != 0) {
               chi2 += std::pow((y - y_fit) / error, 2);
            } else {
               chi2 += 0;
            }
        }
        return chi2;
    };

    // 创建 Minuit2 最小化器
    ROOT::Minuit2::Minuit2Minimizer minimizer3(ROOT::Minuit2::kMigrad);
    minimizer3.SetMaxFunctionCalls(1000000);
    minimizer3.SetMaxIterations(100000);
    minimizer3.SetTolerance(0.001);

    // 设置目标函数
    ROOT::Math::Functor f3(chi2_3, 22);
    minimizer3.SetFunction(f3);

    // 设置初始参数和范围
   minimizer3.SetFixedVariable(0, "param0", comb_params[0]);
   minimizer3.SetFixedVariable(1, "param1", comb_params[1]);
   minimizer3.SetFixedVariable(2, "param2", comb_params[2]);
   minimizer3.SetFixedVariable(3, "param3", comb_params[3]);
   minimizer3.SetFixedVariable(4, "param4", comb_params[4]);
   minimizer3.SetFixedVariable(5, "param5", comb_params[5]);
   minimizer3.SetFixedVariable(6, "param6", comb_params[6]);
   minimizer3.SetFixedVariable(7, "param7", comb_params[7]);
   minimizer3.SetFixedVariable(8, "param8", comb_params[8]);
   minimizer3.SetFixedVariable(9, "param9", comb_params[9]);
   minimizer3.SetFixedVariable(10, "param10", comb_params[10]);
   minimizer3.SetFixedVariable(11, "param11", comb_params[11]);
   minimizer3.SetFixedVariable(12, "param12", comb_params[12]);
   minimizer3.SetFixedVariable(13, "param13", comb_params[13]);
   minimizer3.SetFixedVariable(14, "param14", comb_params[14]);
   minimizer3.SetFixedVariable(15, "param15", comb_params[15]);
   minimizer3.SetFixedVariable(16, "param16", comb_params[16]);
   minimizer3.SetFixedVariable(17, "param17", comb_params[17]);
   minimizer3.SetFixedVariable(18, "param18", comb_params[18]);
   minimizer3.SetFixedVariable(19, "param19", comb_params[19]);
   minimizer3.SetVariable(20, "param20", comb_params[20], 0.01);
   minimizer3.SetVariableLimits(20, 0.0, 100000);
   minimizer3.SetVariable(21, "param21", comb_params[21], 0.01);
   minimizer3.SetVariableLimits(21, 0.0, 100000);



    // 执行迭代最小化
    for (int iter = 0; iter < iteration; ++iter) {
        if (iter != 0) {
            for (int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
                hneg_clone->SetBinError(i, hneg_error[i - 1]);
            }
        }

        minimizer3.Minimize();

        for (int i = 1; i <= hneg_clone->GetNbinsX(); ++i) {
            hneg_error[i - 1] = TMath::Sqrt(fit3->Eval(hneg_clone->GetBinCenter(i)));
        }
    }

    // 获取拟合结果
    const double* xs3 = minimizer3.X();
    const double* errors3 = minimizer3.Errors();

    nsig = xs3[20]*hpr_clone->GetEntries();
   //  std::cout<<"scale_sig="<<fit3->GetParameter(20)<<std::endl;
   //  std::cout<<"scale_sig="<<xs3[20]<<std::endl;
   if (errors3[20] == 0) {
      nsig_error = -1;
   } else {
      nsig_error = errors3[20]*hpr_clone->GetEntries();
   }
   // // Draw fit3 and hneg_clone
   // TCanvas *c = new TCanvas("c", "Fit Results", 1800, 600);
   // c->Divide(3, 1);
   
   // c->cd(1);
   // hneg_clone->Draw();
   // fit3->SetLineColor(kRed);
   // fit3->Draw("same");
   
   // c->cd(2);
   // hpi_clone->Draw();
   // fit2->SetLineColor(kBlue);
   // fit2->Draw("same");

   // c->cd(3);
   // hpr_clone->Draw();
   // fit1->SetLineColor(kGreen);
   // fit1->Draw("same");
   
   // c->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/fit_results.pdf");

    delete fit1;
    delete fit2;
    delete fit3;
}


// void TFracFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double &nsig,double &nsig_error)
// {
//    double Ntot= hneg->GetEntries();
//    TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone")+TString(hbkg->GetTitle()));
//    TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone")+TString(hsig->GetTitle()));
//    TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone")+TString(hneg->GetTitle()));
//    //   hpi_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.22))/hpi_clone->GetBinContent(hpi_clone->FindBin(0.22)));
//    //   hpr_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.88))/hpr_clone->GetBinContent(hpr_clone->FindBin(0.88)));
//    hpi_clone ->Scale(1/hpi_clone->GetSumOfWeights()/hpi_clone->GetBinWidth(1));
//    hpr_clone ->Scale(1/hpr_clone->GetSumOfWeights()/hpr_clone->GetBinWidth(1));
//    hneg_clone->Scale(1/hneg_clone->GetSumOfWeights()/hneg_clone->GetBinWidth(1));
//    TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
//    mc->Add(hpi_clone);
//    mc->Add(hpr_clone);
//    TFractionFitter* fit = new TFractionFitter(hneg_clone, mc); // initialise
//    fit->Constrain(0, 0.0, 1.0);               // constrain fraction 0 to be between 0 and 1
//    fit->Constrain(1, 0.0, 1.0); 
   
//    fit->SetRangeX(1, hneg->GetNbinsX());      // use all bins in the fit
//    Int_t status = fit->Fit();                 // perform the fit
//    std::cout << "fit status: " << status << std::endl;

//    if (status == 0) {                         // check on fit status
//       //TH1F* result = (TH1F*) fit->GetPlot();
//       //hneg->Draw("Ep");
//       //result->Draw("same");
//       fit->GetResult(1,nsig,nsig_error);
//       // nsig=nsig*hneg->GetEntries();
//    }



//    delete fit;
// }

// void TFracFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg,
//                         double &nsig, double &nsig_error,
//                         double tol = 0.01, int max_iter = 10)
// {  
//     RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
//     RooMsgService::instance().setSilentMode(true);
//     gErrorIgnoreLevel = kFatal;
//     const double Ntot = hneg->GetEntries();

//     TH1D *hpi_clone  = (TH1D*)hbkg->Clone(Form("clone") + TString(hbkg->GetTitle()));
//     TH1D *hpr_clone  = (TH1D*)hsig->Clone(Form("clone") + TString(hsig->GetTitle()));
//     TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone") + TString(hneg->GetTitle()));

//     hpi_clone ->Scale(1.0 / hpi_clone->GetSumOfWeights() / hpi_clone->GetBinWidth(1));
//     hpr_clone ->Scale(1.0 / hpr_clone->GetSumOfWeights() / hpr_clone->GetBinWidth(1));
//     hneg_clone->Scale(1.0 / hneg_clone->GetSumOfWeights() / hneg_clone->GetBinWidth(1));

//     double prev_nsig = -1.0;
//     nsig = 0.0;
//     nsig_error = 0.0;

//     for (int iter = 0; iter < max_iter; ++iter) {
//         TObjArray *mc = new TObjArray(2);
//         mc->Add(hpi_clone);
//         mc->Add(hpr_clone);

//         TFractionFitter *fit = new TFractionFitter(hneg_clone, mc);
//         fit->Constrain(0, 0.0, 1.0);
//         fit->Constrain(1, 0.0, 1.0);
//         fit->SetRangeX(1, hneg->GetNbinsX());

//         // ——3. TMinuit 静默 ——  
//         if (gMinuit) {
//             gMinuit->SetPrintLevel(-1);
//             #ifdef TMINUIT_SETNOWARNINGS
//             gMinuit->SetNoWarnings();
//             #endif
//         }

//         TFitResultPtr r = fit->Fit();
//         int status = r.Get() ? r.Get()->Status() : 1;
//         std::cout << "Iter " << iter + 1 << ", status = " << status << std::endl;

//         if (status != 0) {
//             delete fit;
//             delete mc;
//             std::cerr << "Fit failed with status: " << status << std::endl;
//             break;
//         }

//         double frac = 0.0;
//         double frac_err = 0.0;
//         fit->GetResult(1, frac, frac_err);  // Signal 为第1项

//         nsig = frac ;
//         nsig_error = frac_err ;

//         if (prev_nsig >= 0 && std::abs(nsig - prev_nsig) < tol ) {
//             delete fit;
//             delete mc;
//             std::cout << "Converged after " << iter + 1 << " iterations." << std::endl;
//             break;  // 收敛
//         }

//         prev_nsig = nsig;

//         delete fit;
//         delete mc;
//     }
// }

void TFracFit(TH1D *hsig, TH1D *hbkg, TH1D *hneg, double &nsig,double &nsig_error)
{
  int Ntot=hneg->GetEntries();
   TH1D *hpi_clone = (TH1D*)hbkg->Clone(Form("clone")+TString(hbkg->GetTitle()));
   TH1D *hpr_clone = (TH1D*)hsig->Clone(Form("clone")+TString(hsig->GetTitle()));
   TH1D *hneg_clone = (TH1D*)hneg->Clone(Form("clone")+TString(hneg->GetTitle()));
   //   hpi_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.22))/hpi_clone->GetBinContent(hpi_clone->FindBin(0.22)));
   //   hpr_clone->Scale(hneg->GetBinContent(hneg->FindBin(0.88))/hpr_clone->GetBinContent(hpr_clone->FindBin(0.88)));
   hpi_clone ->Scale(1/hpi_clone->GetSumOfWeights()/hpi_clone->GetBinWidth(1));
   hpr_clone ->Scale(1/hpr_clone->GetSumOfWeights()/hpr_clone->GetBinWidth(1));
   hneg_clone->Scale(1/hneg_clone->GetSumOfWeights()/hneg_clone->GetBinWidth(1));
   TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
   mc->Add(hpr_clone);
   mc->Add(hpi_clone);
   TFractionFitter* fit = new TFractionFitter(hneg_clone, mc); // initialise
   fit->Constrain(0, 0.0, 1.0);               // constrain fraction 1 to be between 0 and 1
   fit->Constrain(1, 0.0, 1.0);               // constrain fraction 1 to be between 0 and 1
   fit->SetRangeX(1, hneg->GetNbinsX());      // use all bins in the fit


   ROOT::Fit::Fitter* fitter = fit->GetFitter();
   // 获取当前参数配置
   ROOT::Fit::FitConfig& config = fitter->Config();
    // 直接设置步长
    const int nParams = 2; // 根据您的模型参数数量修改
    for (int i = 0; i < nParams; i++) {
        fitter->Config().ParSettings(i).SetStepSize(0.001);
    }
       if(Ntot==400){
      // 设置第一个参数（索引0）的初始值
    config.ParSettings(0).SetValue(0.5);   // 设置初始值为0.8
    // 设置第二个参数（索引1）的初始值
    config.ParSettings(1).SetValue(0.5);   // 设置初始值为0.2 
  }
  else{
        // 设置第一个参数（索引0）的初始值
    config.ParSettings(0).SetValue(0.2);   // 设置初始值为0.8
    // 设置第二个参数（索引1）的初始值
    config.ParSettings(1).SetValue(0.8);   // 设置初始值为0.2
  }
    // 创建并配置最小化器选项
    ROOT::Math::MinimizerOptions minOptions;
    
    // 设置容差（精度） - 使用正确的 API
    minOptions.SetTolerance(1e-6);
    
    // 设置最大迭代次数 - 使用正确的 API
    minOptions.SetMaxIterations(2000);
    
    // 设置最大函数调用次数
    minOptions.SetMaxFunctionCalls(1000);
    
    // 应用最小化器选项
    fitter->Config().SetMinimizerOptions(minOptions);
    
    // 设置最小化器类型和算法
    fitter->Config().SetMinimizer("Minuit2", "Migrad");


   Int_t status = fit->Fit();                 // perform the fit
   //std::cout << "fit status: " << status << std::endl;

   if (status == 0) {                         // check on fit status
      //TH1F* result = (TH1F*) fit->GetPlot();
      //hneg->Draw("Ep");
      //result->Draw("same");
      fit->GetResult(0,nsig,nsig_error);
      nsig*=hneg->GetEntries();
      nsig_error*=hneg->GetEntries();
   }



   delete fit;
}


inline double sq(double x){return x*x;}

struct TF_FCN: public  ROOT::Math::IBaseFunctionMultiDim {
   // TH1 *ghsig, *ghbkg, *ghdata;
   TTree* gsigTree, *gbkgTree, *gnegTree;
   //int sbin, ebin;
   bool MLfit, fixtot, fit2D;
   mutable int nbinfit;
   double SigTot, BkgTot, NegTot;
   double ntot;
   TF1 *fFuncPi;
   TF1 *fFuncPr;
   //mutable double nsig, nbkg, nbkg2;
   mutable double Chi2;
   TF_FCN( TTree* sigTree, TTree* bkgTree, TTree* negTree, bool maxlikelihood_ = true, bool fixtot_=false ):
      // ghsig(dynamic_cast<TH1*>(hsig->Clone("hsig_clone"))),
      // ghbkg(dynamic_cast<TH1*>(hbkg->Clone("hbkg_clone"))),
      // ghdata(hdata),
      gsigTree(dynamic_cast<TTree*>(sigTree->CloneTree())), gbkgTree(dynamic_cast<TTree*>(bkgTree->CloneTree())),
      gnegTree(dynamic_cast<TTree*>(negTree->CloneTree())),
      //sbin(startbin), ebin(endbin),
      MLfit(maxlikelihood_), fixtot(fixtot_), fit2D(2),
      SigTot( sigTree->GetEntries() ), BkgTot( bkgTree->GetEntries() ), NegTot( negTree->GetEntries() )

   {
      //if(norm1>0){ if(ghsig->GetSumw2N()==0)ghsig->Sumw2();ghsig->Scale(1/norm1);}
      //if(norm2>0){ if(ghbkg->GetSumw2N()==0)ghbkg->Sumw2();ghbkg->Scale(1/norm2);}

      fFuncPi = new TF1("fFuncPi", funcExpGausExp, -1, 3, 7);
      fFuncPr = new TF1("fFuncPr", funcExpGausExp, -1, 3, 7);
   }
   TF_FCN* Clone()const override{ return new TF_FCN(gsigTree, gbkgTree, gnegTree, MLfit, fixtot); }
   virtual ~TF_FCN(){ delete gsigTree; delete gbkgTree; delete gnegTree;delete fFuncPi; delete fFuncPr;}
   //int ndf(){ return nbinfit-1 -1 - !fixtot - bool(ghbkg2);}
   //virtual double Up()const{ return MLfit ? 0.5 : 1.; }
   virtual unsigned int NDim()const override{ return 12; }  // nPars
   virtual double DoEval(double const* p)const override
   {
      double scl_sig, scl_bkg;
      if (not fixtot)
      {
         scl_sig = p[0], scl_bkg = p[1];
      }
      else
      {
         // scl_sig = p[0], scl_bkg = (NegTot - scl_sig*SigTot)/BkgTot;
         scl_sig = p[0], scl_bkg = 1-scl_sig;
      }
      if (MLfit){
         // {std::cout<<scl_sig<<"wowowowowowowow"<<scl_bkg<<std::endl;
         // {std::cout<<"xmin"<<p[7];
         return getNLL(scl_sig, scl_bkg, &p[2]);}
      else
         // return getChi2(scl_sig, scl_bkg, &p[2]);
         return 1e99;
   }
   // double getChi2( double nsig, double nbkg, double const *funcpar)const  {
      
   //    double chi2 = 0, nll = 0;
   //    nbinfit = 0;
   //    for(int i=0;i<8;i++)fFuncPi->SetParameter(i, funcpar[i]);
   //    for(int i=0;i<8;i++)fFuncPr->SetParameter(i, funcpar[i+8]);
   //    //double pi_tot = fFuncPi->Integral(-2.55,3);
   //    //double pr_tot = fFuncPr->Integral(-2.55,3);
   //    for( int i=1; i<ghdata->GetNbinsX(); ++i)
   //    {
   //       double bkg = ghbkg->GetBinContent(i); double e_bkg = ghbkg->GetBinError(i); e_bkg = e_bkg>0?e_bkg:1;
   //       double sig = ghsig->GetBinContent(i); double e_sig = ghsig->GetBinError(i); e_sig = e_sig>0?e_sig:1;
   //       double x = ghdata->GetBinCenter(i);

   //       double y_bkg = fFuncPi->Eval(x);
   //       double y_sig = fFuncPr->Eval(x);
   //       if(TMath::IsNaN(y_bkg)||TMath::IsNaN(y_sig)) continue;
   //       double data = ghdata->GetBinContent(i); double e = ghdata->GetBinError(i); e = e>0?e:1;
   //       double model = nsig * y_sig + nbkg * y_bkg;
         
   //       double chi2_bkg = sq ( (y_bkg-bkg)/e_bkg);
   //       double chi2_sig = sq ( (y_sig-sig)/e_sig);
   //       chi2 += chi2_bkg + chi2_sig;
   //       chi2 += sq( (model-data)/e );
   //       if(bkg>0) nbinfit++;
   //       if(sig>0) nbinfit++;
   //       if(data>0) nbinfit++;
   //       //if(0<x && x<0.10) cout << log(p_data) << ' ' << log(p_bkg) << ' '<<log(p_sig) << endl;
   //    }
   //    //if(not TMath::IsNaN(nll) ) cout << nsig << ' ' << ntot-nsig << ' ' << ' ' << nll << endl;
   //    Chi2 = chi2;
   //    return chi2;
   // }
   double getNLL(  double scl_sig, double scl_bkg, double const *funcpar) const  {
      double chi2 = 0, nll = 0;
      nbinfit = 0;
      for(int i=0;i<5;i++)fFuncPi->SetParameter(i, funcpar[i]);
      for(int i=0;i<5;i++)fFuncPr->SetParameter(i, funcpar[i+5]);
      fFuncPi->FixParameter(5, -1);
      fFuncPr->FixParameter(5, -1);
      fFuncPi->FixParameter(6, 3);
      fFuncPr->FixParameter(6, 3);
      //double pi_tot = fFuncPi->Integral(-2.55,3);
      //double pr_tot = fFuncPr->Integral(-2.55,3);
      // std::cout<<"xmin"<<funcpar[5];
      if( scl_sig<=0 || scl_bkg <=0) return 1e100;
 
      Double_t sig_mass2;
      gsigTree->SetBranchAddress("pr_mass2", &sig_mass2);
      Double_t bkg_mass2;
      gbkgTree->SetBranchAddress("pi_mass2", &bkg_mass2);
      Double_t neg_mass2;
      gnegTree->SetBranchAddress("neg_mass2", &neg_mass2);

      Double_t p_sig = 0, p_bkg = 0, p_neg = 0;
      for(int i=0;i<gsigTree->GetEntries();i++){
         gsigTree->GetEntry(i);
         p_sig+=(fFuncPr->Eval(sig_mass2)>0?std::log(fFuncPr->Eval(sig_mass2)):0);
         //cout<<"sig_mass="<<sig_mass2<<endl;
      }
      for(int i=0;i<gbkgTree->GetEntries();i++){
         gbkgTree->GetEntry(i);
         p_bkg+=(fFuncPi->Eval(bkg_mass2)>0?std::log(fFuncPi->Eval(bkg_mass2)):0);
      }
      for(int i=0;i<gnegTree->GetEntries();i++){
         gnegTree->GetEntry(i);
         p_neg+=(scl_sig*fFuncPr->Eval(neg_mass2)+scl_bkg*fFuncPi->Eval(neg_mass2)>0?std::log(scl_sig*fFuncPr->Eval(neg_mass2)+scl_bkg*fFuncPi->Eval(neg_mass2)):0);
      }
      nll=-2.0*(p_sig+p_bkg+p_neg);
      // std::cout<<fFuncPr->Eval(sig_mass2)<<std::endl;
      // std::cout<<fFuncPi->Eval(bkg_mass2)<<std::endl;
      // cout<<nll<<endl;
      return nll;
   }
};


struct TF_Fitter{
   ROOT::Math::Minimizer *fitter;
   TF_FCN *fcn;
   Float_t nsig, nbkg;
   Float_t chi2, ndf;
   Double_t pars[12];
   Int_t ret;
   TH1D *fHistSig, *fHistBkg;
   TF_Fitter(TTree* sigTree, TTree* bkgTree,TTree* negTree, bool MLFit, bool fixsum, bool Print=false):
   nsig(0.),nbkg(0.),chi2(0.),ndf(0.),ret(0)
   {

      TH1D *hsig = new TH1D("hsig", "hsig", 100, -1, 3);
      TH1D *hbkg = new TH1D("hbkg", "hbkg", 100, -1, 3);
      TH1D *hneg = new TH1D("hneg", "hneg", 100, -1, 3);
      sigTree->Draw("pr_mass2>>hsig", "", "goff");
      bkgTree->Draw("pi_mass2>>hbkg", "", "goff");
      negTree->Draw("neg_mass2>>hneg", "", "goff");
      
      if(!sigTree || !bkgTree || !negTree) return;
      memset(pars,0,sizeof(pars));

      fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
      fcn = new TF_FCN(sigTree,bkgTree,negTree, MLFit, fixsum);
      //if(fcn) cout<<"fcn is ok"<<endl;
      fitter->SetFunction(*fcn);
      // fitter->SetLimitedVariable(0, "Nsig", 0.5, 1e-8, 0., 100);
      // fitter->SetLimitedVariable(1, "Nbkg", 0.5, 1e-8, 0., 100);
      fitter->SetLimitedVariable(0, "Nsig", 0.8, 1e-8, 0., 1.);
      fitter->SetLimitedVariable(1, "Nbkg", 0.2, 1e-8, 0., 1.);
      //fixsum ? fitter->FixVariable(1) : 0;

      TF1 *func_bkg = new TF1("bkg_template",funcExpGausExp,-1,3,7);
      func_bkg->SetParNames("x0","sigmaL","alphaL","sigmaR","alphaR","xmin","xmax");
      func_bkg->SetParameters(hbkg->GetMean(),hbkg->GetRMS(),1,hbkg->GetRMS(),1,-1,3);
      //func_bkg->SetParameters(0.8,0.3,1,1,0.30,1,1);
      func_bkg->SetParLimits(0,0.001,0.35);//x0
      func_bkg->SetParLimits(1,0.005,1.1);//sigL
      func_bkg->SetParLimits(2,0.1,5.0);//alphaL
      func_bkg->SetParLimits(3,0.005,1.1);//sigR
      func_bkg->SetParLimits(4,0.01,5.0);//alphaR
      //func_bkg->SetParLimits(5,0,h->GetMaximum()*2);//norm
      func_bkg->FixParameter(5,-1);
      func_bkg->FixParameter(6,3);
      func_bkg->SetNpx(150);
      hbkg->Scale(1/hbkg->GetSumOfWeights()/hbkg->GetBinWidth(1));
      hbkg->Fit(func_bkg,"LMNQ","R");

      TF1 *func_sig = new TF1("sig_template",funcExpGausExp,-1,3,7);
      func_sig->SetParNames("x0","sigmaL","alphaL","sigmaR","alphaR","xmin","xmax");
      func_sig->SetParameters(hsig->GetMean(),hsig->GetRMS(),1,hsig->GetRMS(),1,-1,3);
      //func_sig->SetParameters(0.8,0.3,1,1,0.30,1,1);
      func_sig->SetParLimits(0,0.5,1.1);//x0
      func_sig->SetParLimits(1,0.005,1.1);//sigL
      func_sig->SetParLimits(2,0.1,5.0);//alphaL
      func_sig->SetParLimits(3,0.005,1.1);//sigR
      func_sig->SetParLimits(4,0.01,5.0);//alphaR
      //func_sig->SetParLimits(5,0,h->GetMaximum()*2);//norm
      func_sig->FixParameter(5,-1);
      func_sig->FixParameter(6,3);
      func_sig->SetNpx(150);
      hsig->Scale(1/hsig->GetSumOfWeights()/hsig->GetBinWidth(1));
      hsig->Fit(func_sig,"LMNQ","R");




      SetParameter(func_bkg, func_sig);
      //cout<<"kkkkkkkkk"<<endl;
      //cout<<"before"<<cb_bkg.GetParameter(0)<<endl;
      // TCanvas *c = new TCanvas("c", "c", 800, 600);
      // c->Divide(2,1);
      // c->cd(1);
      // cb_sig.Draw();
      // c->cd(2);
      // cb_bkg.Draw();
      //fitter->SetMaxFunctionCalls(10000);
      ret = fitter->Minimize();
      if(ret==0 && Print){
         std::cout <<"TF_Fitter failed" << std::endl;
         fitter->PrintResults();
      }
      //cout<<"after"<<fitter->X()[2]<<endl;
      //fitter->PrintResults();
      // c->cd(1);
      // cb_sig.Draw("same");
      // c->cd(2);
      // cb_bkg.Draw("same");
      //c->SaveAs("/data03/tianye/pbar/pic/unbinned_comb.pdf");
      nsig = fitter->X()[0]*negTree->GetEntries();
      nbkg = negTree->GetEntries()-nsig;
      //chi2 = fcn->Chi2;
      //ndf  = fcn->nbinfit-16-!fixsum;
      memcpy(pars, fitter->X(), sizeof(Double_t)*12);
      pars[0] = nsig;
      pars[1] = nbkg;
      // cout<<pars[0]<<endl;
      // cout<<"x[0]="<<fitter->X()[0]<<endl;
      // cout<<"x[1]="<<fitter->X()[1]<<endl;
      hsig->Delete();
      hbkg->Delete();
      hneg->Delete();
   }
   Double_t Nsig()
   {
      return fHistSig?fHistSig->Integral():-1;
      TH1D *hsig2 = (TH1D*)fHistSig->Clone("hsig2");
      hsig2->Reset();
      double pars_sig[8];
      for(int j=0;j<8;j++) 
         pars_sig[j] = fitter->X()[j+10];
      TF1 *cb_sig = new TF1("cb_sig2", CrystalBall, -2.55, 3, 8);
      cb_sig->SetParameters(pars_sig);
      for(int j=1;j<hsig2->GetNbinsX();j++)
         hsig2->SetBinContent(j, cb_sig->Eval(hsig2->GetBinCenter(j)));
      double nsig = hsig2->Integral();
      delete hsig2;
      return nsig;
   }
   void SetParameter(TF1 *func_bkg, TF1 *func_sig)
   {

      fitter->SetLimitedVariable(2, "x0_bkg", func_bkg->GetParameter(0), 0.001, 0.001,0.35);
      fitter->SetLimitedVariable(3, "sL_bkg", func_bkg->GetParameter(1), 0.001, 0.005,1.1);
      fitter->SetLimitedVariable(4, "aL_bkg", func_bkg->GetParameter(2), 0.001, 0.1,5.0);
      fitter->SetLimitedVariable(5, "sR_bkg", func_bkg->GetParameter(3), 0.001, 0.005,1.1);
      fitter->SetLimitedVariable(6, "aR_bkg", func_bkg->GetParameter(4), 0.001, 0.01,5.0);
      // fitter->FixVariable(7);
      // fitter->FixVariable(8);

      fitter->SetLimitedVariable(7,  "x0_sig", func_sig->GetParameter(0), 0.001, 0.5,1.1);
      fitter->SetLimitedVariable(8,  "sL_sig", func_sig->GetParameter(1), 0.001, 0.005,1.1);
      fitter->SetLimitedVariable(9,  "aL_sig", func_sig->GetParameter(2), 0.001, 0.1,5.0);
      fitter->SetLimitedVariable(10, "sR_sig", func_sig->GetParameter(3), 0.001,  0.005,1.1);
      fitter->SetLimitedVariable(11, "aR_sig", func_sig->GetParameter(4), 0.001,  0.01,5.0);
      // fitter->FixVariable(14);
      // fitter->FixVariable(15);

   }
   virtual ~TF_Fitter(){
      delete fcn;
   };

};
