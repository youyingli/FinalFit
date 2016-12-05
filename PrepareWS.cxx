#include <iostream>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"

#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace RooFit;

int main(){

    TString inpath = "/wk_cms/youying/Analysis3/output/";
    const int nfiles = 6;
    TString filename[nfiles];
    //filename[0] = "output_VBFHToGG_M125_13TeV_amcatnlo_pythia8"                                    ;
    filename[0] = "output_VBFHToGG_M126_13TeV_amcatnlo_pythia8"                                    ;
    //filename[1] = "output_GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8"                             ;
    filename[1] = "output_GluGluHToGG_M126_13TeV_amcatnloFXFX_pythia8"                             ;
    filename[2] = "output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"  ;
    filename[3] = "output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8" ;
    filename[4] = "output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa"                                ;
    filename[5] = "output_Run2016BCD_ICHEP_PromptReco_v2"                                          ;  
    TString treename = "outTree";
    TFile* file[nfiles];
    TTree* tree[nfiles];
    for(int i=0;i<nfiles;i++){
        file[i]=TFile::Open( inpath + filename[i] + ".root");
        tree[i] = (TTree*)file[i] ->Get( treename );
    }
    //RooRealVar
    RooRealVar Mgg("Mgg","Mgg",100.,180.);
    RooRealVar BDT("BDT","BDT",0.4,1.0);
    RooRealVar NormWgt("NormWgt","NormWgt",-10000.,10000.);
    //RooDataSet(Signal)
    RooDataSet Hgg_vbf_126("Hgg_vbf_126","Hgg_vbf_126",RooArgSet(Mgg,BDT),Import(*tree[0]));
    RooDataSet Hgg_ggf_126("Hgg_ggf_126","Hgg_ggf_126",RooArgSet(Mgg,BDT),Import(*tree[1]));
    RooDataSet Hgg_vbf_126w("Hgg_vbf_126w","Hgg_vbf_126w",RooArgSet(Mgg,BDT,NormWgt),Import(*tree[0]),WeightVar(NormWgt));
    RooDataSet Hgg_ggf_126w("Hgg_ggf_126w","Hgg_ggf_126w",RooArgSet(Mgg,BDT,NormWgt),Import(*tree[1]),WeightVar(NormWgt));
    //RooDataSet(Background)
    RooDataSet gammajet1("gammajet1","gammajet1",RooArgSet(Mgg,BDT,NormWgt),Import(*tree[2]),WeightVar(NormWgt));
    RooDataSet gammajet2("gammajet2","gammajet2",RooArgSet(Mgg,BDT,NormWgt),Import(*tree[3]),WeightVar(NormWgt));
    RooDataSet dipho("dipho","dipho",RooArgSet(Mgg,BDT,NormWgt),Import(*tree[4]),WeightVar(NormWgt));
    RooDataSet bkg_13TeV(gammajet1,"bkg_13TeV");
    bkg_13TeV.append(gammajet2);
    bkg_13TeV.append(dipho);
    //RooDataSet(data)
    RooDataSet dataset("data","data",RooArgSet(Mgg,BDT),Import(*tree[5]));

    RooWorkspace ws("ws");
    ws.import(Mgg);
    ws.import(BDT);
    ws.import(Hgg_vbf_126);
    ws.import(Hgg_ggf_126);
    ws.import(Hgg_vbf_126w);
    ws.import(Hgg_ggf_126w);
    ws.import(bkg_13TeV);
    ws.import(dataset);

    ws.writeToFile("DataSet.root",true);
    ws.Print("v");

    return 0;
}
