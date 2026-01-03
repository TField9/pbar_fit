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

void table(){
    // 配置
    std::vector<int> size = {10,20,50,100};
    const char* names[] = {
      "RooComb","DA","URComb", "URFit" ,"KP10"};
       TFile *f =new TFile("/data03/tianye/pbar/root/bias_2d/vary_signals_vs_bkg.root","read");
         int nsize=size.size();
         static std::ofstream outfile("/data03/tianye/pbar/syscode/pdf/means.txt");
        // outfile<<"binned bias"<<endl;
        //  for(int i=0;i<nsize;++i){
        //     outfile <<size[i]<<" &";
        //     for(int j=0;j<nsize;++j){
        //         TF1* f1 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[0],size[j],size[i]));
        //         TF1* f2 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[1],size[j],size[i]));
        //         //TF1* f3 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[2],size[j],size[i]));
        //         // 打开文件用于追加写入
                
        //         // outfile << Form("Nb=%d, Ns=%d: URComb mean = %.2f,  URFit mean = %.2f,  KP10 mean = %.2f, URComb sigma=%.2f, URFit sigma=%.2f, KP10 sigma=%.2f\n",
        //         //                 size[i], size[j],
        //         //                 f1 ? (f1->GetParameter(1)/200-1)*100 : 0,
        //         //                 f2 ? (f2->GetParameter(1)/200-1)*100 : 0,
        //         //                 f3 ? (f3->GetParameter(1)/200-1)*100 : 0,
        //         //                 f1 ? (f1->GetParameter(2)/200*100) : 0,
        //         //                 f2 ? f2->GetParameter(2)/2 : 0,
        //         //                 f3 ? f3->GetParameter(2)/2 : 0);
        //         outfile << Form("%.2f & %.2f &",
        //                         f1 ? (f1->GetParameter(1)/200-1)*100 : 0,
        //                         f2 ? (f2->GetParameter(1)/200-1)*100 : 0);
        //     }

        //     outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
        //  }

        //  outfile<<"unbinned bias"<<endl;
        //  for(int i=0;i<nsize;++i){
        //     outfile <<size[i]<<" &";
        //     for(int j=0;j<nsize;++j){
        //         TF1* f1 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[2],size[j],size[i]));
        //         TF1* f2 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[3],size[j],size[i]));
        //         //TF1* f3 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[4],size[j],size[i]));
        //         // 打开文件用于追加写入
                
        //         // outfile << Form("Nb=%d, Ns=%d: URComb mean = %.2f,  URFit mean = %.2f,  KP10 mean = %.2f, URComb sigma=%.2f, URFit sigma=%.2f, KP10 sigma=%.2f\n",
        //         //                 size[i], size[j],
        //         //                 f1 ? (f1->GetParameter(1)/200-1)*100 : 0,
        //         //                 f2 ? (f2->GetParameter(1)/200-1)*100 : 0,
        //         //                 f3 ? (f3->GetParameter(1)/200-1)*100 : 0,
        //         //                 f1 ? (f1->GetParameter(2)/200*100) : 0,
        //         //                 f2 ? f2->GetParameter(2)/2 : 0,
        //         //                 f3 ? f3->GetParameter(2)/2 : 0);
        //         outfile << Form("%.2f & %.2f &",
        //                         f1 ? (f1->GetParameter(1)/200-1)*100 : 0,
        //                         f2 ? (f2->GetParameter(1)/200-1)*100 : 0);
        //     }

        //     outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
        //  }
        //   outfile<<"binned sigma"<<endl;
        // for(int i=0;i<nsize;++i){
        //     outfile <<size[i]<<" &";
        //     for(int j=0;j<nsize;++j){
        //         TF1* f1 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[0],size[j],size[i]));
        //         TF1* f2 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[1],size[j],size[i]));
                
        //         // 打开文件用于追加写入
                
        //         // outfile << Form("Nb=%d, Ns=%d: URComb mean = %.2f,  URFit mean = %.2f,  KP10 mean = %.2f, URComb sigma=%.2f, URFit sigma=%.2f, KP10 sigma=%.2f\n",
        //         //                 size[i], size[j],
        //         //                 f1 ? (f1->GetParameter(1)/200-1)*100 : 0,
        //         //                 f2 ? (f2->GetParameter(1)/200-1)*100 : 0,
        //         //                 f3 ? (f3->GetParameter(1)/200-1)*100 : 0,
        //         //                 f1 ? (f1->GetParameter(2)/200*100) : 0,
        //         //                 f2 ? f2->GetParameter(2)/2 : 0,
        //         //                 f3 ? f3->GetParameter(2)/2 : 0);
        //         outfile << Form("%.2f & %.2f &",
        //                         f1 ? f1->GetParameter(2)/2 : 0,
        //                         f2 ? f2->GetParameter(2)/2 : 0);
        //     }
            
        //     outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
        //  }
        //   outfile<<"unbinned sigma"<<endl;
        // for(int i=0;i<nsize;++i){
        //     outfile <<size[i]<<" &";
        //     for(int j=0;j<nsize;++j){
        //         TF1* f1 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[2],size[j],size[i]));
        //         TF1* f2 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[3],size[j],size[i]));
        //         //TF1* f3 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[4],size[j],size[i]));
        //         // 打开文件用于追加写入
                
        //         // outfile << Form("Nb=%d, Ns=%d: URComb mean = %.2f,  URFit mean = %.2f,  KP10 mean = %.2f, URComb sigma=%.2f, URFit sigma=%.2f, KP10 sigma=%.2f\n",
        //         //                 size[i], size[j],
        //         //                 f1 ? (f1->GetParameter(1)/200-1)*100 : 0,
        //         //                 f2 ? (f2->GetParameter(1)/200-1)*100 : 0,
        //         //                 f3 ? (f3->GetParameter(1)/200-1)*100 : 0,
        //         //                 f1 ? (f1->GetParameter(2)/200*100) : 0,
        //         //                 f2 ? f2->GetParameter(2)/2 : 0,
        //         //                 f3 ? f3->GetParameter(2)/2 : 0);
        //         outfile << Form("%.2f & %.2f &",
        //                         f1 ? f1->GetParameter(2)/2 : 0,
        //                         f2 ? f2->GetParameter(2)/2 : 0);
        //     }
            
        //     outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
        //  }
                 outfile<<"binned error"<<endl;
         for(int i=0;i<nsize;++i){
            outfile <<size[i]<<" &";
            for(int j=0;j<nsize;++j){
                TF1* f1 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[0],size[j],size[i]));
                TF1* f2 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[1],size[j],size[i]));
                //TF1* f3 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[4],size[j],size[i]));
                // 打开文件用于追加写入
                
                // outfile << Form("Nb=%d, Ns=%d: URComb mean = %.2f,  URFit mean = %.2f,  KP10 mean = %.2f, URComb sigma=%.2f, URFit sigma=%.2f, KP10 sigma=%.2f\n",
                //                 size[i], size[j],
                //                 f1 ? (f1->GetParameter(1)/200-1)*100 : 0,
                //                 f2 ? (f2->GetParameter(1)/200-1)*100 : 0,
                //                 f3 ? (f3->GetParameter(1)/200-1)*100 : 0,
                //                 f1 ? (f1->GetParameter(2)/200*100) : 0,
                //                 f2 ? f2->GetParameter(2)/2 : 0,
                //                 f3 ? f3->GetParameter(2)/2 : 0);
                outfile << Form("%.2f & %.2f &",
                                sqrt(TMath::Power((f1->GetParameter(1)/200-1)*100,2)+TMath::Power((f2->GetParameter(2)/200)*100,2)),
                                sqrt(TMath::Power((f2->GetParameter(1)/200-1)*100,2)+TMath::Power((f2->GetParameter(2)/200)*100,2)));
            }

            outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
         }
          outfile<<"unbinned error"<<endl;
        for(int i=0;i<nsize;++i){
            outfile <<size[i]<<" &";
            for(int j=0;j<nsize;++j){
                TF1* f1 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[2],size[j],size[i]));
                TF1* f2 = (TF1*)f->Get(Form("f_%s_fit_%d_%d", names[3],size[j],size[i]));
                
                // 打开文件用于追加写入
                
                // outfile << Form("Nb=%d, Ns=%d: URComb mean = %.2f,  URFit mean = %.2f,  KP10 mean = %.2f, URComb sigma=%.2f, URFit sigma=%.2f, KP10 sigma=%.2f\n",
                //                 size[i], size[j],
                //                 f1 ? (f1->GetParameter(1)/200-1)*100 : 0,
                //                 f2 ? (f2->GetParameter(1)/200-1)*100 : 0,
                //                 f3 ? (f3->GetParameter(1)/200-1)*100 : 0,
                //                 f1 ? (f1->GetParameter(2)/200*100) : 0,
                //                 f2 ? f2->GetParameter(2)/2 : 0,
                //                 f3 ? f3->GetParameter(2)/2 : 0);
                outfile << Form("%.2f & %.2f &",
                                sqrt(TMath::Power((f1->GetParameter(1)/200-1)*100,2)+TMath::Power((f2->GetParameter(2)/200)*100,2)),
                                sqrt(TMath::Power((f2->GetParameter(1)/200-1)*100,2)+TMath::Power((f2->GetParameter(2)/200)*100,2)));
            }
            
            outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
         }
}