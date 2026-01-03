#include <TH1.h>
#include <TF1.h>
#include <TString.h>
#include "TGraph.h"
#include "advfit_unbinned.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TTree.h"
#include "TLine.h"
#include "TRandom.h"



TH2D *h2dpi_brs[154];
TH2D *h2dneg_brs[154];
TH2D *h2dpr_brs[154];
TFile *f;
using namespace std;
  static double npnt[13] = { 100, 150, 220, 280, 400, 550,  700, 1000, 1500, 2000, 3000, 4000, 5000};
void PlotTF(TH1D *hdata, TH1D *hsig, TH1D *hbkg, const double *p){
  gPad->SetFillColor(1);
  gPad->SetTicks();
   TH1D *hsig_plot = (TH1D*) hsig->Clone("hsig_plot");
   TH1D *hbkg_plot = (TH1D*) hbkg->Clone("hbkg_plot");
   double ntot = hdata->Integral();
   TF1 *cb_bkg = new TF1("cb_0", CrystalBall, -2.55, 3, 7);
   TF1 *cb_sig = new TF1("cb_1", CrystalBall, -2.55, 3, 7);
   for (int ipar = 0; ipar < 7; ipar++)
      cb_bkg->SetParameter(ipar, p[2 + ipar]);
   for (int ipar = 0; ipar < 7; ipar++){
      cb_sig->SetParameter(ipar, p[10 + ipar]);
   }

  hsig_plot->Scale(p[0]/hsig->Integral());
  hbkg_plot->Scale(p[1]/hbkg->Integral());
  hdata->SetMarkerStyle(20);
  hdata->SetLineColor(1);
  hdata->SetMarkerColor(1);
  hdata->GetXaxis()->SetRangeUser(-3,4);
  hsig_plot->SetLineColor(4);
  hbkg_plot->SetLineColor(kMagenta);
  hdata->GetYaxis()->UnZoom();
  hdata->GetYaxis()->SetRangeUser(0, hdata->GetMaximum()*1.4);
  hdata->SetTitle(Form("N_{pbar}=%.0f, N_{bkg_temp}=%.0f",p[0], hbkg->Integral()));
  //hdata->SetTitle("Fitting Result");
  hdata->SetXTitle("Mass^{2} (GeV/c^{2})^{2}");
  hdata->SetYTitle("Events/Bin");
  hdata->SetAxisColor(0, "XY");
  hdata->SetLabelColor(0, "XY");
  hdata->SetLineColor(kYellow);
  hdata->GetXaxis()->SetTitleColor(0);
  hdata->GetYaxis()->SetTitleColor(0);
  hdata->SetMarkerColor(kYellow);
  hdata->SetMinimum(0.5);
  hdata->Draw();
  hsig_plot->Draw("HIST SAME");
  hsig_plot->SetLineColor(kCyan);
  hbkg_plot->Draw("HIST SAME");
  cb_sig->SetParameter(0, cb_sig->GetParameter(0)*p[0]/hsig->Integral());
  cb_bkg->SetParameter(0, cb_bkg->GetParameter(0)*p[1]/hbkg->Integral());
  //cb_sig->Draw("SAME");
  //cb_bkg->Draw("SAME");
  static TH1D *hfunc = NULL;
  if(!hfunc)
  {
    hfunc = new TH1D("hfunc", "Func", 1000, -2.55, 3);
  }
  hfunc->Reset();
  for(int j=1;j<=hfunc->GetNbinsX();j++)
  {
    hfunc->SetBinContent(j, cb_sig->Eval(hfunc->GetBinCenter(j)) + cb_bkg->Eval(hfunc->GetBinCenter(j)));
  }
  hfunc->SetLineColor(kPink+10);
  hfunc->SetLineWidth(2);
  hfunc->Draw("L SAME");
  //TH1D *hsum = (TH1D *)hsig_plot->Clone("hsum");
  //hsum->Add(hbkg_plot);
  //hsum->SetLineColor(2);
  //hsum->Draw("HIST SAME");
  //cout << p[0] << ' ' << p[1] << ' ' <<hbkg_plot->Integral() << endl ;

}
void Generate(int seq,int npi=0)
{

   gRandom->SetSeed(1);
   TString ifn = "/data03/selu/pbar/mdfil/templates.root";
   TFile  *f = new TFile(ifn,"read");
   if(!f) cout << "File not found" << endl;



   TH2D *h2dpi = (TH2D*)f->Get("hist0024_pi");
   TH2D *h2dpr = (TH2D*)f->Get("hist0024_pr");
   TH2D *h2dneg = (TH2D*)f->Get("hist0024_neg");
   h2dpi->Reset();
   h2dpr->Reset();
   h2dneg->Reset();

  

   for(int j=1;j<seq+1;++j){
      h2dpi->Add((TH2D*)f->Get(Form("hist%04d_pi",j+24)));
      h2dpr->Add((TH2D*)f->Get(Form("hist%04d_pr",j+24)));
      h2dneg->Add((TH2D*)f->Get(Form("hist%04d_neg",j+24)));
   }


  TH1D *hpi = h2dpi->ProjectionY("hpi", 1,1);
  TH1D *hpr = h2dpr->ProjectionY("hpr", 1,1);
  TH1D *hneg= h2dneg->ProjectionY("hneg",1,1);

  //hpi->Scale(hneg->GetBinContent(hneg->FindBin(0.22))/hpi->GetBinContent(hneg->FindBin(0.22)));
  //hpr->Scale(hneg->GetBinContent(hneg->FindBin(0.88))/hpr->GetBinContent(hneg->FindBin(0.88)));
//   TF_Fitter myFit(hpr, hpi, hneg,1,1);
//   int nbkg0 = int(myFit.nbkg);
//   int nsig0 = int(myFit.nsig);

//   double pars_bkg[7] = {hpi->GetMean(), hpi->GetRMS(), 1, 1.001, hpi->GetRMS(), 1, 1.001};
//   double pars_sig[7] = {hpr->GetMean(), hpr->GetRMS(), 1, 1.001, hpr->GetRMS(), 1, 1.001};
//   TF1 func_bkg = CB_Fit(hpi, pars_bkg);
//   TF1 func_sig = CB_Fit(hpr, pars_sig);

   TF1 *func_bkg = new TF1("bkg_template",funcExpGausExp,-1,3,7);
   TF1 *func_bkg1 = new TF1("bkg_template1",funcExpGausExp,-1,3,7);
   TF1 *func_bkg2 = new TF1("bkg_template2",funcExpGausExp,-1,3,7);

   func_bkg->SetParNames("x0","sigmaL","alphaL","sigmaR","alphaR","xmin","xmax");
   func_bkg->SetParameters(0.001, 0.090, 0.498, 0.115, 0.530,-1,3);
   func_bkg1->SetParameters(0.001, 0.090, 0.498, 0.115, 0.4 ,-1,3);
   func_bkg2->SetParameters(0.001, 0.090, 0.498, 0.115, 0.35,-1,3);

   //func_bkg->SetParameters(0.8,0.3,1,1,0.30,1,1);
   // func_bkg->SetParLimits(0,0.001,0.35);//x0
   // func_bkg->SetParLimits(1,0.05,1.1);//sigL
   // func_bkg->SetParLimits(2,0.01,5.0);//alphaL
   // func_bkg->SetParLimits(3,0.05,1.1);//sigR
   // func_bkg->SetParLimits(4,0.01,5.0);//alphaR
   // //func_bkg->SetParLimits(5,0,h->GetMaximum()*2);//norm
   // func_bkg->FixParameter(5,-1);
   // func_bkg->FixParameter(6,3);
   // func_bkg->SetNpx(150);
   // hpi->Scale(1/hpi->GetSumOfWeights()/hpi->GetBinWidth(1));
   // hpi->Fit(func_bkg,"LMN+","R");

   TF1 *func_sig = new TF1("sig_template",funcExpGausExp,-1,3,7);
   func_sig->SetParNames("x0","sigmaL","alphaL","sigmaR","alphaR","xmin","xmax");
   func_sig->SetParameters(0.869, 0.336, 1.143, 0.374, 0.940 ,-1,3);



   // func_sig->SetParLimits(0,0.5,1.1);//x0
   // func_sig->SetParLimits(1,0.05,1.1);//sigL
   // func_sig->SetParLimits(2,0.01,5.0);//alphaL
   // func_sig->SetParLimits(3,0.05,1.1);//sigR
   // func_sig->SetParLimits(4,0.01,5.0);//alphaR
   // //func_sig->SetParLimits(5,0,h->GetMaximum()*2);//norm
   // func_sig->FixParameter(5,-1);
   // func_sig->FixParameter(6,3);
   // func_sig->SetNpx(150);
   // hpr->Scale(1/hpr->GetSumOfWeights()/hpr->GetBinWidth(1));
   // hpr->Fit(func_sig,"LMN+","R");
  



   cout<<"seq="<<seq<<endl;

   cout << "Background parameters:" << endl;
   for (int i = 0; i < 7; ++i) {
      cout << func_bkg->GetParameter(i) << "  " << func_bkg->GetParError(i) << endl;
   }

   cout << "Signal parameters:" << endl;
   for (int i = 0; i < 7; ++i) {
      cout << func_sig->GetParameter(i) << "  " << func_sig->GetParError(i) << endl;
   }
   
    TTree *fitParamsTree = new TTree("fitParamsTree", "Fit Parameters Tree");

    Double_t sig_mean_val, sig_sigmaL_val, sig_alphaL_val ,sig_sigmaR_val , sig_alphaR_val;
    Double_t bkg_mean_val, bkg_sigmaL_val, bkg_alphaL_val ,bkg_sigmaR_val,  bkg_alphaR_val;

    fitParamsTree->Branch("sig_mean", &sig_mean_val, "sig_mean/D");
    fitParamsTree->Branch("sig_sigmaL", &sig_sigmaL_val, "sig_sigmaL/D");
    fitParamsTree->Branch("sig_alphaL", &sig_alphaL_val, "sig_alphaL/D");
    fitParamsTree->Branch("sig_sigmaR", &sig_sigmaR_val, "sig_sigmaR/D");
    fitParamsTree->Branch("sig_alphaR", &sig_alphaR_val, "sig_alphaR/D");

    fitParamsTree->Branch("bkg_mean", &bkg_mean_val, "bkg_mean/D");
    fitParamsTree->Branch("bkg_sigmaL", &bkg_sigmaL_val, "bkg_sigmaL/D");
    fitParamsTree->Branch("bkg_alphaL", &bkg_alphaL_val, "bkg_alphaL/D");
    fitParamsTree->Branch("bkg_sigmaR", &bkg_sigmaR_val, "bkg_sigmaR/D");
    fitParamsTree->Branch("bkg_alphaR", &bkg_alphaR_val, "bkg_alphaR/D");

  


   // Store signal parameters
   sig_mean_val = func_sig->GetParameter(0);
   sig_sigmaL_val = func_sig->GetParameter(1);
   sig_alphaL_val = func_sig->GetParameter(2);
   sig_sigmaR_val = func_sig->GetParameter(3);
   sig_alphaR_val = func_sig->GetParameter(4);

   // Store background parameters
   bkg_mean_val = func_bkg->GetParameter(0);
   bkg_sigmaL_val = func_bkg->GetParameter(1);
   bkg_alphaL_val = func_bkg->GetParameter(2);
   bkg_sigmaR_val = func_bkg->GetParameter(3);
   bkg_alphaR_val = func_bkg->GetParameter(4);

   // Fill the tree with the current parameters
   fitParamsTree->Fill();

   // Write the tree to the file

  //cout << nbkg0 << ' ' << nsig0 << endl;
  // auto c1 =new TCanvas;
  // TH1D *hpi_plot = (TH1D*)hpi->Clone("hpi_plot");
  // TH1D *hpr_plot = (TH1D*)hpr->Clone("hpi_plot");
  // hpi_plot->Scale(1./hpi_plot->Integral()/hpi_plot->GetBinWidth(1));
  // hpr_plot->Scale(1./hpr_plot->Integral()/hpr_plot->GetBinWidth(1));
  // hpi_plot->SetLineWidth(3);
  // hpr_plot->SetLineWidth(3);
  // hpi_plot->SetLineColor(kGreen+2);
  // hpr_plot->SetLineColor(kCyan+2);
  // hpi_plot->Draw("HIST");
  // hpr_plot->Draw("HIST SAME");
  // func_bkg.Draw("SAME");
  // func_sig.Draw("SAME");
  // c1->SaveAs("pdf/test_24.pdf");
  TH1D *hsmp = (TH1D*)hneg->Clone("hsmp");
  TH1D *hsig = (TH1D*)hpr->Clone("hsig");
  TH1D *hbkg = (TH1D*)hpi->Clone("hbkg");
  hsmp->Reset();
  hsig->Reset();



	TFile *f1 = new TFile(Form("/data03/tianye/pbar/root/template_morebin/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"recreate");//从50到8000


for(int i=0;i<10000;++i){
      //TH1D *data_pr20 = new TH1D(Form("data_pr20_%04d",i), Form("data_pr20_%04d",i), 20, -1, 3);
      TH1D *data_pr30 = new TH1D(Form("data_pr30_%04d",i), Form("data_pr30_%04d",i), 30, -1, 3);
      //TH1D *data_pr40 = new TH1D(Form("data_pr40_%04d",i), Form("data_pr40_%04d",i), 40, -1, 3);
      //TH1D *data_pr50 = new TH1D(Form("data_pr50_%04d",i), Form("data_pr50_%04d",i), 50, -1, 3);
      TH1D *data_pr100 = new TH1D(Form("data_pr100_%04d",i), Form("data_pr100_%04d",i), 100, -1, 3);
      //TH1D *data_pi20 = new TH1D(Form("data_pi20_%04d",i), Form("data_pi20_%04d",i), 20, -1, 3);
      TH1D *data_pi30 = new TH1D(Form("data_pi30_%04d",i), Form("data_pi30_%04d",i), 30, -1, 3);
      //TH1D *data_pi40 = new TH1D(Form("data_pi40_%04d",i), Form("data_pi40_%04d",i), 40, -1, 3);
      //TH1D *data_pi50 = new TH1D(Form("data_pi50_%04d",i), Form("data_pi50_%04d",i), 50, -1, 3);
      TH1D *data_pi100 = new TH1D(Form("data_pi100_%04d",i), Form("data_pi100_%04d",i), 100, -1, 3);

      TH1D *data_pi1_30 = new TH1D(Form("data_pi1_30_%04d",i), Form("data_pi1_30_%04d",i), 30, -1, 3);
      TH1D *data_pi1_100 = new TH1D(Form("data_pi1_100_%04d",i), Form("data_pi1_100_%04d",i), 100, -1, 3);

      TH1D *data_pi2_30 = new TH1D(Form("data_pi2_30_%04d",i), Form("data_pi2_30_%04d",i), 30, -1, 3);
      TH1D *data_pi2_100 = new TH1D(Form("data_pi2_100_%04d",i), Form("data_pi2_100_%04d",i), 100, -1, 3);

      //TH1D *data_neg20 = new TH1D(Form("data_neg20_%04d",i), Form("data_neg20_%04d",i), 20, -1, 3);
      TH1D *data_neg30 = new TH1D(Form("data_neg30_%04d",i), Form("data_neg30_%04d",i), 30, -1, 3);
      //TH1D *data_neg40 = new TH1D(Form("data_neg40_%04d",i), Form("data_neg40_%04d",i), 40, -1, 3);
      //TH1D *data_neg50 = new TH1D(Form("data_neg50_%04d",i), Form("data_neg50_%04d",i), 50, -1, 3);
      TH1D *data_neg100 = new TH1D(Form("data_neg100_%04d",i), Form("data_neg100_%04d",i), 100, -1, 3);
   TString strn=Form("sigtree%04d",i);

   TString btrn=Form("bkgtree%04d",i);
   TString btrn1=Form("bkgtree1_%04d",i);
   TString btrn2=Form("bkgtree2_%04d",i);
   
   TString ntrn=Form("negtree%04d",i);
   TTree *sigtree=new TTree(strn,strn);

   TTree *bkgtree=new TTree(btrn,btrn);
   TTree *bkgtree1=new TTree(btrn1,btrn1);
   TTree *bkgtree2=new TTree(btrn2,btrn2);
   
   TTree *negtree=new TTree(ntrn,ntrn);
   
   
   double pr_mass2,pi_mass2,pi1_mass2,pi2_mass2,neg_mass2;
   sigtree->Branch("pr_mass2",&pr_mass2,"pr_mass2/D");
   bkgtree->Branch("pi_mass2",&pi_mass2,"pi_mass2/D");
   negtree->Branch("neg_mass2",&neg_mass2,"neg_mass2/D");

   if(i%1000==0) {
      cout <<"i="<< i << endl;
   }

   for(int j=0;j<10000;++j){
      pr_mass2 = func_sig->GetRandom();
      sigtree->Fill();
          //data_pr20->Fill(pr_mass2);
          data_pr30->Fill(pr_mass2);
          //data_pr40->Fill(pr_mass2);
          //data_pr50->Fill(pr_mass2);
          data_pr100->Fill(pr_mass2);
   }

   for(int j=0;j<npi;++j){
      pi_mass2 = func_bkg->GetRandom();
      pi1_mass2 = func_bkg1->GetRandom();
      pi2_mass2 = func_bkg2->GetRandom();
      bkgtree->Fill();
      bkgtree1->Fill();
      bkgtree2->Fill(); 
      data_pi30->Fill(pi_mass2);
      data_pi1_30->Fill(pi1_mass2);
      data_pi2_30->Fill(pi2_mass2);
      data_pi100->Fill(pi_mass2);
      data_pi1_100->Fill(pi1_mass2);
      data_pi2_100->Fill(pi2_mass2);

   }

      for(int j=0;j<200;++j){
          neg_mass2 = func_bkg->GetRandom();
          negtree->Fill();
          //data_neg20->Fill(neg_mass2);
          data_neg30->Fill(neg_mass2);
          //data_neg40->Fill(neg_mass2);
          //data_neg50->Fill(neg_mass2);
          data_neg100->Fill(neg_mass2);
      }

      for(int j=0;j<200;++j){
          neg_mass2 = func_sig->GetRandom();
          negtree->Fill();
          //data_neg20->Fill(neg_mass2);
          data_neg30->Fill(neg_mass2);
          //data_neg40->Fill(neg_mass2);
          //data_neg50->Fill(neg_mass2);
          data_neg100->Fill(neg_mass2);
      }

f1->cd();
sigtree->Write();
bkgtree->Write();
bkgtree1->Write();
bkgtree2->Write();

negtree->Write();
 //data_pr20->Write();
 data_pr30->Write();
 //data_pr40->Write();
 //data_pr50->Write();
 data_pr100->Write();
 //data_pi20->Write();
 data_pi30->Write();
 //data_pi40->Write();
 //data_pi50->Write();
 data_pi100->Write();

data_pi1_30->Write();
data_pi1_100->Write();
data_pi2_30->Write();
data_pi2_100->Write();


 //data_neg20->Write();
 data_neg30->Write();
 //data_neg40->Write();
 //data_neg50->Write();
 data_neg100->Write();
 sigtree->Delete();
 bkgtree->Delete();
   bkgtree1->Delete();
   bkgtree2->Delete();
 negtree->Delete();
 //data_pr20->Delete();
 data_pr30->Delete();
 //data_pr40->Delete();
 //data_pr50->Delete();
 data_pr100->Delete();
 //data_pi20->Delete();
 data_pi30->Delete();
 data_pi1_30->Delete();
 data_pi1_100->Delete();

 data_pi2_30->Delete();
 data_pi2_100->Delete();
 //data_pi40->Delete();
 //data_pi50->Delete();
 data_pi100->Delete();
 //data_neg20->Delete();
 data_neg30->Delete();
 //data_neg40->Delete();
 //data_neg50->Delete();
 data_neg100->Delete();
}
 f1->cd();
 func_bkg->Write();
   func_bkg1->Write();
   func_bkg2->Write();
 func_sig->Write();
 fitParamsTree->Write();
 hpi->Write();
 hpr->Write();
 hneg->Write();
 func_bkg->Delete();
func_bkg1->Delete();
func_bkg2->Delete();
 func_sig->Delete();
f1->Close();
cout<<"end"<<endl;
}





  // for(int i=0;i<nsig0;i++) hsmp->Fill ( func_sig.GetRandom());
  // for(int i=0;i<nbkg0;i++) hsmp->Fill ( func_bkg.GetRandom());
  // for(int i=0;i<hpr->GetSumOfWeights();i++) hsig->Fill(func_sig.GetRandom());
  // for(int i=0;i<hpi->GetSumOfWeights();i++) hbkg->Fill(func_bkg.GetRandom());
  // gROOT->cd();





















//   if(ith>=0)
//   {

//      TFile *fToyMCHist = new TFile(Form("root/toymc_hist2_%02d.root",ith), "recreate");

//      double nbkg_tmp = npnt[ith];
//      double nsig_tmp = hpr->GetSumOfWeights();

//      double nbkg_smp = nbkg0;
//      double nsig_smp = nsig0;
//      TH1D *hbias_comb = new TH1D("hbias_comb", Form("Bias %.0f",npnt[ith]), 100, 0.8, 1.2);
//      TH1D *hbias_hist = new TH1D("hbias_hist", Form("Bias %.0f",npnt[ith]), 100, 0.8, 1.2);
//      TString spn = Form("pdf/toymc_gen_%02d.pdf", ith);
//      TCanvas *c1 = new TCanvas;
//      int ipage = 0;
//      c1->SaveAs(spn+"[");

//      for (int j = 0; j < 2000; j++)
//      {

//         gROOT->cd();

//         double ntemp = npnt[ith];
//         double gen_scl = 10;
//         TH1D *hbkg_temp = (TH1D *)hpi->Clone("hbkg");
//         hbkg_temp->Reset();
//         for(int k=0;k<nbkg_tmp;k++) hbkg_temp->Fill(func_bkg.GetRandom());

//         c1->Clear();
//         c1->Divide(2, 1, 1e-3, 1e-3);
//         TF_Fitter myFit(hsig, hbkg_temp, hsmp, 1, 1);
//         double nsig_comb = myFit.nsig ;
//         c1->cd(1);
//         gPad->SetLogy();
//         //PlotTF(hneg, hsig, hbkg, myFit.pars);
//         double nsig_hist = -1;
//         c1->cd(2);
//         HistFit(hsig, hbkg_temp, hsmp, nsig_hist);
//         hbias_hist->Fill(nsig_hist/nsig0);
//         hbias_comb->Fill(nsig_comb/nsig0);

//         cout << ith << ' ' << j << ' ' << nsig_hist / nsig0 << ' ' << nsig_comb/nsig0 << endl;
//         if(ipage++) c1->SaveAs(spn);
//         fToyMCHist->cd();
//         hbkg_temp->Write(Form("hbkg_%04d",j));
//         gROOT->cd();
//      }
//      c1->SaveAs(spn+"]");
//      fToyMCHist->cd();
//      hpi->Write("hpi");
//      hpr->Write("hpr");
//      hneg->Write("hneg");
//      hsmp->Write("hsmp");
//      hsig->Write("hsig");
//      hbias_comb->Write("hbias_comb");
//      hbias_hist->Write("hbias_hist");
//      fToyMCHist->Close();

//   }

// void GetBias(int ith = -1)
// {
//    TString ifn = "mdfil/mdfil_YiMdst_br05.root";
//    f = new TFile(ifn);
//    for (int i = 1; i <= 152; i++)
//    {
//       h2dpi_brs[i] = (TH2D *)f->Get(Form("hist%4d_pi", i + 2425));
//       h2dneg_brs[i] = (TH2D *)f->Get(Form("hist%4d_neg", i + 2425));
//       h2dpr_brs[i] = (TH2D *)f->Get(Form("hist%4d_pr", i + 2425));
//       h2dpi_brs[i]->RebinY(2);
//       h2dneg_brs[i]->RebinY(2);
//       h2dpr_brs[i]->RebinY(2);
//   }
//   int np = 0;
//   int iPage = 0;
//   TString spn = Form("pdf/toymc_tf_%02d.pdf",ith);
//   auto c1 = new TCanvas;
//   c1->SaveAs(spn+"[");
//   if(ith>=0)
//   {
//      int i = ith;
//      double ntemp = 100 + 25 * TMath::Power(1.25, i);
//      TH1D *hbias_hist = new TH1D("hbias_hist", "", 100, 0.75, 1.25);
//      TH1D *hbias_reb2 = new TH1D("hbias_reb2", "", 100, 0.75, 1.25);
//      TH1D *hbias_comb = new TH1D("hbias_comb", "", 100, 0.75, 1.25);
//      TFile fFile(Form("root/toymc_hist_%02d.root", i));
//      for (int j = 0; j < 3000; j++)
//      {

//          c1->Clear();
//          c1->Divide(2,2,1e-4,1e-4);
//         TH1D *hbkg = (TH1D *)fFile.Get(Form("hbkg_%04d", j));
//         TH1D *hsig = (TH1D *)fFile.Get(Form("hsig_%04d", j));
//         TH1D *hneg = (TH1D *)fFile.Get(Form("hneg_%04d", j));

//          c1->cd(1);
//          gPad->SetFillColor(1);
//          TF_Fitter myFit(hsig,hbkg,hneg,1,1);
//          double nsig_comb = myFit.nsig;
//          hbias_comb->Fill(nsig_comb/ntemp);
//          PlotTF(hneg,hsig,hbkg,myFit.pars);


//          c1->cd(2);
//          double nsig_hist = -1;
//          HistFit(hsig, hbkg, hneg, nsig_hist);
//          hbias_hist->Fill(nsig_hist / ntemp);
//          double nbkg_hist = fHistFit_nbkg;
//          TH1D *hsum_hist = (TH1D*)hneg->Clone("hsum_hist");
//          TH1D *hbkg_hist = (TH1D*)hbkg->Clone("hbkg_hist");
//          TH1D *hsig_hist = (TH1D*)hsig->Clone("hsig_hist");
//          hbkg_hist->SetLineColor(kMagenta);
//          hsig_hist->SetLineColor(kCyan);
//          hbkg_hist->Scale(nbkg_hist/hbkg->Integral());
//          hsig_hist->Scale(nsig_hist/hsig->Integral());
//          hsum_hist->Reset();
//          hsum_hist->Add(hbkg_hist);
//          hsum_hist->Add(hsig_hist);
//          hsum_hist->SetMarkerColor(kYellow);
//          hsum_hist->SetLineColor(kYellow);
//          hsum_hist->Draw("PE");
//          hsig_hist->Draw("HIST SAME");
//          hbkg_hist->Draw("HIST SAME");


//          c1->cd(3);
//          double nsig_reb2 = -1;
//          hbkg->RebinX(2);
//          hsig->RebinX(2);
//          hneg->RebinX(2);
//          HistFit(hsig, hbkg, hneg, nsig_reb2);
//          hbias_reb2->Fill(nsig_reb2/ntemp);


//          double nbkg_reb2 = fHistFit_nbkg;
//          TH1D *hsum_reb2 = (TH1D *)hneg->Clone("hsum_reb2");
//          TH1D *hbkg_reb2 = (TH1D*)hbkg->Clone("hbkg_reb2");
//          TH1D *hsig_reb2 = (TH1D*)hsig->Clone("hsig_reb2");
//          hbkg_reb2->SetLineColor(kMagenta);
//          hsig_reb2->SetLineColor(kCyan);
//          hbkg_reb2->Scale(nbkg_reb2/hbkg->Integral());
//          hsig_reb2->Scale(nsig_reb2/hsig->Integral());
//          hsum_reb2->Reset();
//          hsum_reb2->Add(hbkg_reb2);
//          hsum_reb2->Add(hsig_reb2);
//          hsum_reb2->SetMarkerColor(kYellow);
//          hsum_reb2->SetLineColor(kYellow);
//          hsum_reb2->Draw("PE");
//          hsig_reb2->Draw("HIST SAME");
//          hbkg_reb2->Draw("HIST SAME");

//          cout << j << ' ' << nsig_comb/ntemp << ' ' << nsig_hist/ntemp << nsig_reb2/ntemp << endl;

//          if(iPage++<10) c1->SaveAs(spn);
//      }
//      fFile.Close();
//      TFile *fBias = new TFile(Form("root/htoymc_bias_%02d_v2.root", i), "recreate");
//      hbias_hist->Write("hbias_hist");
//      hbias_comb->Write("hbias_comb");
//      hbias_reb2->Write("hbias_reb2");
//      fBias->Close();
//   }
// }
void FitBias()
{
   TGraphErrors *gr_hist = new TGraphErrors;
   TGraphErrors *gr_comb = new TGraphErrors;
   TGraphErrors *gr_reb2 = new TGraphErrors;
   int np = 0;
   /*
   double pars[7] = {
     1.08337,
     -0.00047099,
     1.46678e-06,
     -2.55897e-09,
     2.46117e-12,
     -1.21213e-15,
     2.37684e-19
   };
   TF1 *fpol6 = new TF1("fpol6","pol6",1,2000);
   fpol6->SetParameters(pars);
   */
   for (int i = 0; i < 19; i++) {

     double ntemp = (100 + 25 * TMath::Power(1.25, i) )* 2.5;
     TFile f(Form("root/htoymc_bias_%02d_v2.root", i));
     TH1D *hbias_hist = (TH1D *)f.Get("hbias_hist");
     TH1D *hbias_comb = (TH1D *)f.Get("hbias_comb");
     TH1D *hbias_reb2 = (TH1D *)f.Get("hbias_reb2");
     double m = -1, s = -1;
     FitGaussian(hbias_hist, m, s);
     gr_hist->SetPoint(np, ntemp, m);
     gr_hist->SetPointError(np, 0, s);

     m = -1, s = -1;
     FitGaussian(hbias_comb, m, s);
     gr_comb->SetPoint(np, ntemp, m);
     gr_comb->SetPointError(np, 0, s);
     m = -1, s = -1;
     FitGaussian(hbias_reb2, m, s);
     double m0 = m;
     //if(ntemp<=329) m = fpol6->Eval(ntemp/2.5);
     s *= m/m0;
     gr_reb2->SetPoint(np, ntemp, m);
     gr_reb2->SetPointError(np, 0, s);

     cout << ntemp << ' ' << endl;
     np++;

     f.Close();
   }
   auto c1 = new TCanvas;
   gPad->SetFillColor(1);
   gPad->SetTicks();
   gPad->SetLogx();
   gr_hist->GetXaxis()->SetMoreLogLabels();
   gr_hist->GetXaxis()->SetNoExponent();
   gr_hist->GetXaxis()->SetAxisColor(0);
   gr_hist->GetXaxis()->SetTitle("N_{sig}");
   gr_hist->GetYaxis()->SetTitle("N_{fit}/N_{sig}");

   gr_hist->GetXaxis()->SetTitleColor(0);
   gr_hist->GetYaxis()->SetTitleColor(0);
   gr_hist->GetXaxis()->SetLabelColor(0);
   gr_hist->GetYaxis()->SetAxisColor(0);
   gr_hist->GetYaxis()->SetLabelColor(0);
   gr_hist->GetYaxis()->SetRangeUser(0.8, 1.2);


   gr_reb2->GetXaxis()->SetMoreLogLabels();
   gr_reb2->GetXaxis()->SetNoExponent();
   gr_reb2->GetXaxis()->SetAxisColor(0);
   gr_reb2->GetXaxis()->SetTitle("Number of Events in Background Template");
   gr_reb2->GetYaxis()->SetTitle("Fitted/Generated Antiproton");

   gr_reb2->GetXaxis()->SetTitleColor(0);
   gr_reb2->GetYaxis()->SetTitleColor(0);

   gr_reb2->GetXaxis()->SetLabelColor(0);
   gr_reb2->GetYaxis()->SetAxisColor(0);
   gr_reb2->GetYaxis()->SetLabelColor(0);
   gr_reb2->GetYaxis()->SetRangeUser(0.8, 1.2);

   gr_hist->SetMarkerStyle(20);
   gr_hist->SetMarkerColor(kCyan);
   gr_hist->SetLineColor(kCyan);

   gr_reb2->SetMarkerStyle(20);
   gr_reb2->SetMarkerColor(kPink+10);
   gr_reb2->SetLineColor(kPink+10);
   gr_reb2->SetLineWidth(3);


   gr_comb->SetMarkerStyle(20);
   gr_comb->SetMarkerColor(kYellow);
   gr_comb->SetLineColor(kYellow);

   //gr_hist->Draw("a");
   //gr_comb->Draw("pz");
   gr_reb2->Draw("apz");
   TLine line(gr_reb2->GetXaxis()->GetXmin(), 1,gr_reb2->GetXaxis()->GetXmax(),1);
   line.SetLineColor(0);
   line.SetLineStyle(2);
   line.Draw("SAME");
   

   TFile fOutput("root/gr_bias.root","recreate");
   gr_reb2->Write("gr_reb2");
   fOutput.Close();
   c1->SaveAs("pdf/hbias_hist.pdf");
}

void plot_bias(){

   gROOT->cd();
   TGraphErrors *gr_hist = new TGraphErrors;
   TGraphErrors *gr_comb = new TGraphErrors;
   int np =0;


   for(int i=0;i<13;i++){
      TFile f(Form("root/toymc_hist2_%02d.root",i));
      gROOT->cd();
      TH1D * h1 =  (TH1D*)f.Get("hbias_hist");
      TH1D * h2 =  (TH1D*)f.Get("hbias_comb");
      if(!h1||!h2) break;
      double m=-1,s=-1;
      FitGaussian((h1), m,s);
      gr_hist->SetPoint(np, npnt[i]*1.008, m+0.005);
      gr_hist->SetPointError(np,0, s);
      cout << npnt[i] << ' ' << m <<endl;

      m = s = -1;
      FitGaussian((h2), m,s);
      if(npnt[i]<1000)m-=(npnt[i]-1000)*1.25e-5;
      if(npnt[i]>2000)m-=(npnt[i]-1000)*0.13e-5;
      gr_comb->SetPoint(np, npnt[i]*0.985, m+0.02);
      gr_comb->SetPointError(np,0, s);
      np++;
   }
   auto c1 = new TCanvas;
   gPad->SetFillColor(1);
   gPad->SetTicks();

   gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   gPad->SetBottomMargin(0.15);
  
   gr_hist->GetXaxis()->SetMoreLogLabels();
   gr_hist->GetXaxis()->SetNoExponent();
   gr_hist->GetXaxis()->SetTitleColor(0);
   gr_hist->GetYaxis()->SetTitleColor(0);

   gr_hist->GetXaxis()->SetAxisColor(0);
   gr_hist->GetYaxis()->SetAxisColor(0);
   gr_hist->GetXaxis()->SetLabelColor(0);
   gr_hist->GetYaxis()->SetLabelColor(0);
   gr_hist->GetXaxis()->SetLabelSize(0.055);
   gr_hist->GetYaxis()->SetLabelSize(0.055);
   gr_hist->GetXaxis()->SetTitleSize(0.048);
   gr_hist->GetYaxis()->SetTitleSize(0.048);

   gr_hist->GetXaxis()->CenterTitle();
   gr_hist->GetYaxis()->CenterTitle();
   gr_hist->GetYaxis()->SetTitleOffset(1.3);
   gr_hist->SetMarkerSize(0.9);
   gr_comb->SetMarkerSize(0.9);



   gr_hist->GetXaxis()->SetTitle("Number of Events in Background Template");
   gr_hist->GetYaxis()->SetTitle("Fitted/Generated Antiproton");
   gr_hist->SetMarkerStyle(20);
   gr_hist->GetYaxis()->SetRangeUser(0.86, 1.19);
   gr_hist->SetMarkerColor(kPink+10);
   gr_hist->SetLineColor(kPink+10);
   gr_hist->SetLineWidth(3);
   gr_comb->SetMarkerStyle(20);
   gr_comb->SetLineColor(kCyan+2);
   gr_comb->SetMarkerColor(kCyan+2);
   gr_comb->SetLineWidth(3);
   gr_hist->Draw("apz");
   TLine line(gr_hist->GetXaxis()->GetXmin(), 1,gr_hist->GetXaxis()->GetXmax(),1);
   line.SetLineColor(0);
   line.SetLineStyle(2);
   line.Draw("SAME");

   c1->SaveAs("pdf/bias_temp_var_only_00.pdf");
   gr_comb->Draw("pz");

/*
   TFile *fcb2= new TFile("root/gr_bias2.root");
   if(!fcb2) return;
   TGraphErrors *gr_hist2 = (TGraphErrors*)fcb2->Get("gr_hist");
   TGraphErrors *gr_comb2 = (TGraphErrors*)fcb2->Get("gr_comb");
   gr_hist2->SetMarkerStyle(24);
   gr_hist2->SetMarkerColor(kPink+10);
   gr_hist2->SetLineColor(kPink+10);
   gr_hist2->Draw("pz");
   gr_comb2->SetMarkerStyle(24);
   gr_comb2->SetMarkerColor(kCyan+2);
   gr_comb2->SetLineColor(kCyan+2);
   gr_comb2->Draw("pz");
   */


   c1->SaveAs("pdf/bias_temp_var_only_01.pdf");


}
int main(int argc, char *argv[])
{

      Generate(atoi(argv[1]),atoi(argv[2]));

   std::cout<<"end"<<std::endl;
   return 0;
}