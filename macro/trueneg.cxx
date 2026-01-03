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

void trueneg(){
   TString ifn = "/data03/selu/pbar/mdfil/templates.root";
   TFile  *f = new TFile(ifn,"read");
   TH2D *h2dneg= (TH2D*)f->Get("hist0024_neg");
   TH1D *hneg= h2dneg->ProjectionY("hneg",1,1);
   TCanvas *c1 = new TCanvas("c1","c1",800,600);
   hneg->SetTitle("True Negative Mass Distribution");
   hneg->GetXaxis()->SetTitle("Mass (GeV/c^2)");
   hneg->GetYaxis()->SetTitle("Counts");
   hneg->SetLineColor(kRed);
   hneg->SetLineWidth(2);
   hneg->Draw();
   c1->SaveAs("/data03/tianye/pbar/syscode/pdf/trueneg.pdf");
}