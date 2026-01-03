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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include <vector>

void oned_bias(){
    // 配置
    std::vector<int> size = {10,15,20,25,30,40,50,75,100,200,400};
    const char* names[] = {
        "RooComb","TFF","IWLS","DA","JSC","asy","URFit","URComb" ,
        "KP0.7","KP0.9","KP1.0","KP1.1","KP1.3"

    };

    int base_methods = 8; // 
    int nsize=size.size();
    const std::vector<double> widths = {0.7, 0.9, 1.0, 1.1, 1.3};
    int key_pdf_widths = widths.size();
    int met = base_methods + key_pdf_widths;
    TH1D* h_fit[met];
    TF1 *f_fit[met];
    TGraphErrors *g_bias[met][nsize], *g_sigma[met][nsize];




    for(int i=0;i<met;++i){
        for(int j=0;j<nsize;j++){
            TString bias;
            TString sigma;
            if(i<base_methods){
                bias = Form("g_%s_bias_Signal_%d", names[i], size[j]);
                sigma = Form("g_%s_sigma_Signal_%d", names[i], size[j]);
            }
            else{
                bias= Form("g_KP%02d_bias_Signal_%d", int(widths[i - base_methods] * 10), size[j]);
                sigma= Form("g_KP%02d_sigma_Signal_%d", int(widths[i - base_methods] * 10), size[j]);
            }



            g_bias[i][j]=new TGraphErrors();
            g_bias[i][j]->SetName(bias);
            g_bias[i][j]->SetTitle(Form("%s;N_{B};Bias(%%)", bias.Data()));

            g_sigma[i][j] = new TGraphErrors();
            g_sigma[i][j]->SetName(sigma);
            g_sigma[i][j]->SetTitle(Form("%s;N_{B};Sigma(%%)", sigma.Data()));
        }
    }

    TFile *f =new TFile("/data03/tianye/pbar/root/bias_2d/vary_signals_vs_bkg.root","recreate");


    for(int npr=0;npr<nsize;npr++){
        for(int npi=0;npi<nsize;npi++){

            for(int i=0;i<met;++i){
                if(i<base_methods){
                    h_fit[i]=new TH1D(Form("h_%s_fit_%d_%d", names[i],size[npr],size[npi]), Form("h_%s_fit_%d_%d", names[i],npr,npi), 400, 0, 400);
                    f_fit[i] = new TF1(Form("f_%s_fit_%d_%d", names[i],size[npr],size[npi]), "gaus", 0, 400);
                }
                else if(i>=base_methods && i<met){
                    h_fit[i]=new TH1D(Form("h_KP%02d_fit_%d_%d", int(widths[i - base_methods] * 10),size[npr],size[npi]), Form("h_KP%02d_fit_%d_%d", int(widths[i - base_methods] * 10),npr,npi), 400, 0, 400);
                    f_fit[i] = new TF1(Form("f_KP%02d_fit_%d_%d", int(widths[i - base_methods] * 10),size[npr], size[npi]), "gaus", 0, 400);
                }
            }


            TFile *f1=new TFile(Form("/data03/tianye/pbar/root/2d_result/binned_roocomb_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            TFile *f2=new TFile(Form("/data03/tianye/pbar/root/2d_result/binned_TFF_seq0001_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            TFile *f3=new TFile(Form("/data03/tianye/pbar/root/2d_result/binned_IWLS_seq0001_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            TFile *f4=new TFile(Form("/data03/tianye/pbar/root/2d_result/binned_OM_seq0001_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            TFile *f5=new TFile(Form("/data03/tianye/pbar/root/2d_result/unbinned_roofit_seq0001_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            TFile *f6=new TFile(Form("/data03/tianye/pbar/root/2d_result/unbinned_roocomb_seq0001_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            TFile *f7=new TFile(Form("/data03/tianye/pbar/root/2d_result/unbinned_rookeyspdf_diffwidth_seq0001_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            

            

                    // 取树
            TTree *t_roocomb = (TTree*)f1   ->Get("fitParamsTree");
            TTree *t_TFF     = (TTree*)f2    ->Get("fitParamsTree");
            TTree *t_IWLS    = (TTree*)f3   ->Get("fitParamsTree");
            TTree *t_OM      = (TTree*)f4    ->Get("fitParamsTree");
            TTree *t_unRF    = (TTree*)f5   ->Get("fitParamsTree");
            TTree *t_unRC    = (TTree*)f6   ->Get("fitParamsTree");
            TTree *t_keypdf  = (TTree*)f7   ->Get("fitParamsTree");
            

                        // 参数容器
            std::vector<double>* p_roocomb = new std::vector<double>(2);
            std::vector<double>* p_TFF     = new std::vector<double>(2);
            std::vector<double>* p_IWLS    = new std::vector<double>(2);
            Double_t pDA[2]={0}, pJSC[2]={0}, pasy[2]={0};
            std::vector<double>* p_unRF    = new std::vector<double>(2);
            std::vector<double>* p_unRC    = new std::vector<double>(2);
            std::vector<std::vector<double>>* p_keypdf = new std::vector<std::vector<double>>(widths.size());


            // 绑定
            t_roocomb->SetBranchAddress("nsig", &p_roocomb);
            t_TFF    ->SetBranchAddress("nsig", &p_TFF);
            t_IWLS   ->SetBranchAddress("nsig", &p_IWLS);
            t_OM->SetBranchAddress("Ysig_DA",   pDA);
            t_OM->SetBranchAddress("Ysig_JSC",  pJSC);
            t_OM->SetBranchAddress("Ysig_asy",  pasy);
            t_unRF   ->SetBranchAddress("nsig", &p_unRF);
            t_unRC   ->SetBranchAddress("nsig", &p_unRC);
            t_keypdf ->SetBranchAddress("nsig", &p_keypdf);


            //////////must reset//////////

            for(int i = 0; i < met; ++i){
                h_fit[i]->Reset();
            }


            for(int j=0;j<10000;++j){
                t_roocomb->GetEntry(j);
                t_TFF    ->GetEntry(j);
                t_IWLS   ->GetEntry(j);
                t_OM     ->GetEntry(j);
                t_unRF   ->GetEntry(j);
                t_unRC   ->GetEntry(j);
                t_keypdf ->GetEntry(j);

                // 填充
                h_fit[0]->Fill((*p_roocomb)[0]);
                h_fit[1]->Fill((*p_TFF)[0]);
                h_fit[2]->Fill((*p_IWLS)[0]);
                h_fit[3]->Fill(pDA[0]);
                h_fit[4]->Fill(pJSC[0]);
                h_fit[5]->Fill(pasy[0]);
                h_fit[6]->Fill((*p_unRF)[0]);
                h_fit[7]->Fill((*p_unRC)[0]);

                for(int w=0; w<key_pdf_widths; ++w){
                    h_fit[base_methods + w]->Fill((*p_keypdf)[w][0]);
                }
                

            }

            for(int i = 0; i < met; ++i) {
                f_fit[i]->SetParameters(h_fit[i]->GetMean(), h_fit[i]->GetRMS());
                h_fit[i]->Fit(f_fit[i], "Q0");
                g_bias[i][npr]->SetPoint(npi, size[npi], (f_fit[i]->GetParameter(1)/200-1)*100);
                g_bias[i][npr]->SetPointError(npi, size[npi], (f_fit[i]->GetParError(1)/200-1)*100);
                g_sigma[i][npr]->SetPoint(npi, size[npi], f_fit[i]->GetParameter(2)/200*100);
                g_sigma[i][npr]->SetPointError(npi, size[npi],f_fit[i]->GetParError(2)/200*100);
            }
            f->cd();
            
            
            for(int i=0;i<met;++i){
                h_fit[i]->Write();
                f_fit[i]->Write();
            }

            


        }
        cout<<"npr="<<npr<<endl;
        
    }
    f->cd();

    for(int i=0;i<met;++i){
        for(int j=0;j<nsize;j++){
            g_bias[i][j]->Write();
            g_sigma[i][j]->Write();

        }
    }
    f->Close();


}