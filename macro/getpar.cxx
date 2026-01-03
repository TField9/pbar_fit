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

void getpar(){

TFile *f = new TFile("/data03/tianye/pbar/root/template_bin25/toymc_bin1_seq0001_npi0100.root","read");
TF1 * sig=(TF1*)f->Get("sig_template");
TF1 * bkg=(TF1*)f->Get("bkg_template");
if (sig) {
    std::cout << "Signal function parameters:" << std::endl;
    for (int i = 0; i < sig->GetNpar(); ++i) {
        std::cout << "Parameter " << i << ": " << sig->GetParameter(i) << std::endl;
    }
} else {
    std::cerr << "Signal function not found!" << std::endl;
}

if (bkg) {
    std::cout << "Background function parameters:" << std::endl;
    for (int i = 0; i < bkg->GetNpar(); ++i) {
        std::cout << "Parameter " << i << ": " << bkg->GetParameter(i) << std::endl;
    }
} else {
    std::cerr << "Background function not found!" << std::endl;
}


}