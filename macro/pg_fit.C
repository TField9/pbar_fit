#include <vector>
#include <cmath>
#include <iostream>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TString.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "poisson_gamma.h"  // provided interface

void generalized_pg_mixture(int k, double *alphas, double *betas, size_t w_size, double *result)
{
    int i=0,j=0;

    double first_var_vec[w_size];

    double deltas[k+2], sum_terms[k+1];
    deltas[0]=1.0;

    double prefac=1.0;
    double running_vec[w_size];

    for (i=0; i < w_size;i++)
    {   
        first_var_vec[i]=1.0/(1.0+betas[i]);
        prefac*=pow((1.0/(1.0+1.0/betas[i])), alphas[i]);
        running_vec[i]=1.0;
    }

    if(k>0)
    {
        for(i=1; i<k+1; i++)
        {
            sum_terms[i]=0.0;

            for(j=0; j<w_size;j++)
            {
                running_vec[j]*=first_var_vec[j];
                sum_terms[i]+=alphas[j]*running_vec[j];
            }

            deltas[i]=0.0;
            for(j=1;j<=i; j++)
            {
                deltas[i]+=sum_terms[j]*deltas[i-j];
            }
            deltas[i]/=(double)(i);

        }
        
    }
    *result=prefac*deltas[k];
    
}
/* eq. 96 - generalized mixture */
void generalized_pg_mixture_marginalized_combined(int k,double *new_alphas,  double *betas,  double *gammas, double *alphas_2, double *betas_2, size_t w_size,size_t w_size_2, double *result)
{

    int i=0,j=0;

    double E=0, o_o_m_c=0;

    double first_var_vec[w_size];
    double sec_var_vec[w_size];
    double old_var_vec[w_size_2];

    double deltas[k+2], sum_terms[k+1];
    deltas[0]=1.0;

    double prefac=1.0;

    double running_first_vec[w_size];
    double running_second_vec[w_size];

    double running_old_vec[w_size_2];

    for (i=0; i < w_size;i++)
    {   
        E=1.0/(1.0+betas[i]);
        o_o_m_c=(1.0+betas[i]/(1.0+gammas[i]*(1.0+betas[i])));

        prefac*=pow( (1.0/(1.0+1.0/gammas[i]))*o_o_m_c, new_alphas[i]);

        first_var_vec[i]=E*o_o_m_c;
        sec_var_vec[i]=E;
    
        running_first_vec[i]=1.0;
        running_second_vec[i]=1.0;
        
    }

    for(i=0; i < w_size_2;i++)
    {
        prefac*=pow( 1.0/(1.0+1.0/betas_2[i]), alphas_2[i]);
        old_var_vec[i]=1.0/(1.0+betas_2[i]);
        running_old_vec[i]=1.0;
    }


    if(k>0)
    {
        for(i=1; i<k+1; i++)
        {
            sum_terms[i]=0.0;

            for(j=0; j<w_size;j++)
            {
                running_first_vec[j]*=first_var_vec[j];
                running_second_vec[j]*=sec_var_vec[j];
                
                sum_terms[i]+=new_alphas[j]*running_first_vec[j] - new_alphas[j]*running_second_vec[j];
                //printf("obp, %f , obp ma , %f, firstvarvec: %f, secvarvec: %f\n", o_p_b, o_p_b_minus_a, first_var_vec[j], sec_var_vec[j]);
            }
            for(j=0; j<w_size_2;j++)
            {
                running_old_vec[j]*=old_var_vec[j];
                sum_terms[i]+=alphas_2[j]*running_old_vec[j];
            }
          
            deltas[i]=0.0;
            for(j=1;j<=i; j++)
            {
                deltas[i]+=sum_terms[j]*deltas[i-j];
                //printf("SUMMING IN DELTAS ... , %f, %f\n", sum_terms[j], deltas[i-j]);
            }
            deltas[i]/=(double)(i);

            
        }
        
    }
    *result=prefac*deltas[k];
    
}
class PG_Likelihood {
public:
  PG_Likelihood(const TH1D* sig, const TH1D* bkg, const TH1D* data)
    : nbins(sig->GetNbinsX()),
      Ssum(sig->Integral()), Bsum(bkg->Integral()),
      Ntot(data->Integral()) {
    nSigMC.resize(nbins);
    nBkgMC.resize(nbins);
    nobs .resize(nbins);
    for (int i = 1; i <= nbins; ++i) {
      nSigMC[i-1] = static_cast<int>(sig->GetBinContent(i) + 0.5);
      nBkgMC[i-1] = static_cast<int>(bkg->GetBinContent(i) + 0.5);
      nobs   [i-1] = static_cast<int>(data->GetBinContent(i) + 0.5);
    }
  }

  // operator: x[0]=Nsig
  double operator()(const double *x) const {
    double Nsig = x[0];
    if (Nsig < 0 || Nsig > Ntot) return 1e6;
    double Nbkg = Ntot - Nsig;
    double w_sig = (Ssum>0 ? Nsig/Ssum : 0.0);
    double w_bkg = (Bsum>0 ? Nbkg/Bsum : 0.0);
    double beta_sig = (w_sig>0 ? (1.0 - w_sig)/w_sig : 1e9);
    double beta_bkg = (w_bkg>0 ? (1.0 - w_bkg)/w_bkg : 1e9);

    double nll = 0.0;
    for (int i = 0; i < nbins; ++i) {
      int k  = nobs[i];
      int ns = nSigMC[i];
      int nb = nBkgMC[i];
      int total = ns + nb;
      // 构造统一权重数组：等权情形
      std::vector<double> alphas(total, 1.0);
      std::vector<double> betas(total);
      for (int j = 0; j < ns; ++j) betas[j] = beta_sig;
      for (int j = ns; j < total; ++j) betas[j] = beta_bkg;
      double Li = 0;
      generalized_pg_mixture(k,
                             alphas.data(), betas.data(), total,
                             &Li);
      if (Li <= 0) return 1e6;
      nll -= std::log(Li);
    }
    return nll;
  }

private:
  int nbins;
  double Ssum, Bsum, Ntot;
  std::vector<int> nSigMC, nBkgMC, nobs;
};

// Fit 函数
void Likelihood(TH1D *hsig, TH1D *hbkg, TH1D *hdata,
                double &nsig, double &nsig_error) {
  PG_Likelihood pg(hsig, hbkg, hdata);
  ROOT::Math::Functor functor(
    [&](const double *x){ return pg(x); }, 1);
  ROOT::Minuit2::Minuit2Minimizer minim(ROOT::Minuit2::kMigrad);
  minim.SetFunction(functor);
  minim.SetMaxFunctionCalls(50000);
  minim.SetPrecision(1e-6);
  double init = hdata->Integral() * 0.2;
  minim.SetLimitedVariable(0, "Nsig", init, init*0.01, 0.0, hdata->Integral());
  minim.Minimize();
  nsig       = minim.X()[0];
  nsig_error = std::sqrt(minim.CovMatrix(0,0));
}

// 主循环
void fit_start(int seq, int npi, int low, int high) {
  TStopwatch timer; timer.Start();
  TFile *fin = TFile::Open(Form("/data03/tianye/pbar/root/10percent/toymc_bin1_seq%04d_npi%04d.root", seq, npi));
  TFile *fout = TFile::Open(Form("/data03/tianye/pbar/root/sys_result/binned_pg_seq%04d_npi%04d_%04d.root", seq, npi, low), "RECREATE");
  TTree *tree = new TTree("fitParamsTree", "Fit Params");
  std::vector<double> vec(8);
  tree->Branch("nsig", &vec);
  for (int j = low; j <= high; ++j) {
    TH1D *hs0 = (TH1D*)fin->Get(Form("data_pr30_%04d", j));
    TH1D *hb0 = (TH1D*)fin->Get(Form("data_pi30_%04d", j));
    TH1D *hd0 = (TH1D*)fin->Get(Form("data_neg30_%04d", j));
    TH1D *hs1 = (TH1D*)fin->Get(Form("data1_pr30_%04d", j));
    TH1D *hb1 = (TH1D*)fin->Get(Form("data1_pi30_%04d", j));
    TH1D *hd1 = (TH1D*)fin->Get(Form("data1_neg30_%04d", j));
    Likelihood(hs0, hb0, hd0, vec[0], vec[4]);
    Likelihood(hs0, hb0, hd1, vec[1], vec[5]);
    Likelihood(hs1, hb1, hd0, vec[2], vec[6]);
    Likelihood(hs1, hb1, hd1, vec[3], vec[7]);
    tree->Fill();
    if (j % 10 == 0) {
      double e = timer.RealTime();
      double p = double(j - low + 1) / double(high - low + 1);
      std::cout<< "sig="<< vec[0] << " sig1=" << vec[1] << " sig2=" << vec[2] << " sig3=" << vec[3] << "\n"; 
      std::cout<<"hs1="<<hs1->Integral()<<" hb1="<<hb1->Integral()<<" hd1="<<hd1->Integral()<<"\n";
      std::cout << "j=" << j << " rem=" << (e / p - e) << "s\n";
    }
  }
  fout->cd(); tree->Write(); fout->Close(); fin->Close();
}

int main(int argc, char** argv) {
  if (argc != 5) { std::cerr << "Usage: ... seq npi low high\n"; return 1; }
  fit_start(std::atoi(argv[1]), std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]));
  return 0;
}
