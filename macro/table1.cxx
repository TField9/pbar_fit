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

void table1(){
    // 配置
    std::vector<int> size = {10,15,20,25,30,40,50,75,100,150,200,300,400,500,700,1000,1500,2000,3000,4000,5000,8000};
    std::vector<int> seq = {0,2,6,8,13,15};
    
    // const char* names[] = {
    //   "RooComb","DA","URComb", "URFit" ,"KP1.0"};
    //    TFile *f =new TFile("/data03/tianye/pbar/root/bias_graph/bias_graph.root","read");
    //      int nsize=size.size();
    //      static std::ofstream outfile("/data03/tianye/pbar/syscode/pdf/means1.txt");

    //      outfile<<"binned bias(200:200 and 50:200 Nb)"<<endl;
    //   for(int i=0;i<seq.size();++i){
    //     outfile <<size[seq[i]]<<" &";
    //     for(int j=0;j<2;j++){
    //       TGraphErrors* g1 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[0], j));
    //       TGraphErrors* g2 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[1], j));
    //       outfile<< Form("%.2f & %.2f  &",
    //                             g1 ? (g1->GetY()[seq[i]]-1)*100 : 0,
    //                             g2 ? (g2->GetY()[seq[i]]-1)*100 : 0);
    //     }
    //     outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
    //      }

    //   outfile<<"binned bias(200:200 and 200:50 Ns)"<<endl;
    //   for(int i=0;i<seq.size();++i){
    //     outfile <<size[seq[i]]<<" &";
    //     for(int j=2;j<4;j++){
    //       TGraphErrors* g1 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[0], j));
    //       TGraphErrors* g2 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[1], j));
    //       outfile<< Form("%.2f & %.2f  &",
    //                             g1 ? (g1->GetY()[seq[i]]-1)*100 : 0,
    //                             g2 ? (g2->GetY()[seq[i]]-1)*100 : 0);
    //     }
    //     outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
    //      }


    //      outfile<<"unbinned bias(200:200 and 50:200 Nb)"<<endl;
    //   for(int i=0;i<seq.size();++i){
    //     outfile <<size[seq[i]]<<" &";
    //     for(int j=0;j<2;j++){
    //       TGraphErrors* g1 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[2], j));
    //       TGraphErrors* g2 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[3], j));
    //       TGraphErrors* g3 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[4], j));
    //       outfile<< Form("%.2f & %.2f & %.2f &",
    //                             g1 ? (g1->GetY()[seq[i]]-1)*100 : 0,
    //                             g2 ? (g2->GetY()[seq[i]]-1)*100 : 0,
    //                             g3 ? (g3->GetY()[seq[i]]-1)*100 : 0);
    //     }
    //     outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
    //   }

    //     outfile<<"unbinned bias(200:200 and 200:50 Ns)"<<endl;
    //     for(int i=0;i<seq.size();++i){
    //     outfile <<size[seq[i]]<<" &";
    //     for(int j=2;j<4;j++){
    //       TGraphErrors* g1 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[2], j));
    //       TGraphErrors* g2 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[3], j));
    //       TGraphErrors* g3 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[4], j));
    //       outfile<< Form("%.2f & %.2f & %.2f &",
    //                             g1 ? (g1->GetY()[seq[i]]-1)*100 : 0,
    //                             g2 ? (g2->GetY()[seq[i]]-1)*100 : 0,
    //                             g3 ? (g3->GetY()[seq[i]]-1)*100 : 0);
    //     }
    //     outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
    //      }


        const char* names[] = {
      "RooComb","DA","URComb", "URFit" ,"KP1.0"};
       TFile *f =new TFile("/data03/tianye/pbar/root/bias_graph/bias_graph.root","read");
         int nsize=size.size();
         static std::ofstream outfile("/data03/tianye/pbar/syscode/pdf/means1.txt");

         outfile<<"binned error(200:200 nb)"<<endl;
      for(int i=0;i<seq.size();++i){
        outfile <<size[seq[i]]<<" &";
        for(int j=0;j<1;j++){
          TGraphErrors* g0 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[0], j));
          TGraphErrors* g1 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[1], j));
          TGraphErrors* g2 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[2], j));
          TGraphErrors* g3 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[3], j));
          TGraphErrors* g4 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[4], j));

          TGraphErrors* g0_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[0], j));
          TGraphErrors* g1_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[1], j));
          TGraphErrors* g2_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[2], j));
          TGraphErrors* g3_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[3], j));
          TGraphErrors* g4_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[4], j));



          outfile << Form("%.2f & %.2f & %.2f & %.2f & %.2f &",
                sqrt(pow((g0->GetY()[seq[i]] - 1) * 100, 2) + pow(g0_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g1->GetY()[seq[i]] - 1) * 100, 2) + pow(g1_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g2->GetY()[seq[i]] - 1) * 100, 2) + pow(g2_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g3->GetY()[seq[i]] - 1) * 100, 2) + pow(g3_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g4->GetY()[seq[i]] - 1) * 100, 2) + pow(g4_sig->GetY()[seq[i]] * 100, 2)));
        }
        outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
         }

                  outfile<<"binned error(200:200 ns)"<<endl;
      for(int i=0;i<seq.size();++i){
        outfile <<size[seq[i]]<<" &";
        for(int j=2;j<3;j++){
          TGraphErrors* g0 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[0], j));
          TGraphErrors* g1 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[1], j));
          TGraphErrors* g2 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[2], j));
          TGraphErrors* g3 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[3], j));
          TGraphErrors* g4 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[4], j));

          TGraphErrors* g0_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[0], j));
          TGraphErrors* g1_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[1], j));
          TGraphErrors* g2_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[2], j));
          TGraphErrors* g3_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[3], j));
          TGraphErrors* g4_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[4], j));




outfile << Form("%.2f & %.2f & %.2f & %.2f & %.2f &",
                sqrt(pow((g0->GetY()[seq[i]] - 1) * 100, 2) + pow(g0_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g1->GetY()[seq[i]] - 1) * 100, 2) + pow(g1_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g2->GetY()[seq[i]] - 1) * 100, 2) + pow(g2_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g3->GetY()[seq[i]] - 1) * 100, 2) + pow(g3_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g4->GetY()[seq[i]] - 1) * 100, 2) + pow(g4_sig->GetY()[seq[i]] * 100, 2)));
        }
        outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
         }

        
         outfile<<"binned error(50:200 nb)"<<endl;

      for(int i=0;i<seq.size();++i){
        outfile <<size[seq[i]]<<" &";
        for(int j=1;j<2;j++){
          TGraphErrors* g0 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[0], j));
          TGraphErrors* g1 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[1], j));
          TGraphErrors* g2 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[2], j));
          TGraphErrors* g3 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[3], j));
          TGraphErrors* g4 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[4], j));

          TGraphErrors* g0_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[0], j));
          TGraphErrors* g1_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[1], j));
          TGraphErrors* g2_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[2], j));
          TGraphErrors* g3_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[3], j));
          TGraphErrors* g4_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[4], j));



outfile << Form("%.2f & %.2f & %.2f & %.2f & %.2f &",
                sqrt(pow((g0->GetY()[seq[i]] - 1) * 100, 2) + pow(g0_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g1->GetY()[seq[i]] - 1) * 100, 2) + pow(g1_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g2->GetY()[seq[i]] - 1) * 100, 2) + pow(g2_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g3->GetY()[seq[i]] - 1) * 100, 2) + pow(g3_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g4->GetY()[seq[i]] - 1) * 100, 2) + pow(g4_sig->GetY()[seq[i]] * 100, 2)));
        }
        outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
         }

        
         outfile<<"binned error(200:50 nb)"<<endl;

      for(int i=0;i<seq.size();++i){
        outfile <<size[seq[i]]<<" &";
        for(int j=3;j<4;j++){
          TGraphErrors* g0 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[0], j));
          TGraphErrors* g1 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[1], j));
          TGraphErrors* g2 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[2], j));
          TGraphErrors* g3 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[3], j));
          TGraphErrors* g4 = (TGraphErrors*)f->Get(Form("g_bias_%s_sys%d", names[4], j));

          TGraphErrors* g0_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[0], j));
          TGraphErrors* g1_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[1], j));
          TGraphErrors* g2_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[2], j));
          TGraphErrors* g3_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[3], j));
          TGraphErrors* g4_sig = (TGraphErrors*)f->Get(Form("g_sigma_%s_sys%d", names[4], j));




outfile << Form("%.2f & %.2f & %.2f & %.2f & %.2f &",
                sqrt(pow((g0->GetY()[seq[i]] - 1) * 100, 2) + pow(g0_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g1->GetY()[seq[i]] - 1) * 100, 2) + pow(g1_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g2->GetY()[seq[i]] - 1) * 100, 2) + pow(g2_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g3->GetY()[seq[i]] - 1) * 100, 2) + pow(g3_sig->GetY()[seq[i]] * 100, 2)),
                sqrt(pow((g4->GetY()[seq[i]] - 1) * 100, 2) + pow(g4_sig->GetY()[seq[i]] * 100, 2)));
        }
        outfile <<"\\\\"<<std::endl<<"\\hline"<< std::endl;
         }
         


}