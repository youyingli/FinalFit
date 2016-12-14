#include <iostream>
#include <vector>
#include <map>
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "RooExtendPdf.h"
#include "RooCBShape.h"
#include "RooMCStudy.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMath.h"

#include "../macro/rootlogon.C"
#include "PDFLIB.h"
#include "envelope.h"


using namespace std;
using namespace RooFit;

void makeSigPlot(RooRealVar* Mgg,vector<RooAbsPdf*> Gaussians ,string proc, RooDataSet* MC){

    rootlogon();

    int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
    RooPlot* sframe = Mgg->frame(Bins(40),Range(115.,135.));
    MC->plotOn(sframe,MarkerStyle(20));
    TObject *pdfLeg0 = sframe->getObject(int(sframe->numItems())-1);
    TLegend leg(0.7,0.7,1.,0.92);
    leg.AddEntry(pdfLeg0,Form("%s_126",proc.c_str()),"LP");

    int order = 1;
    for(vector<RooAbsPdf*>::iterator it = Gaussians.begin();it != Gaussians.end();++it){
	(*it)->plotOn(sframe,LineColor(color[order-1]));
	TObject *pdfLeg = sframe->getObject(int(sframe->numItems()-1));
	leg.AddEntry(pdfLeg,Form("%d order",order),"L");
	order++;
    }
    TCanvas canv("canv","canv",700,600);
    TPad pad("pad","",0.,0.,0.98,0.98);
    canv.cd();
    pad.SetLeftMargin(0.20);
    pad.SetRightMargin(0.03);
    pad.Draw();
    pad.cd();
    sframe->SetXTitle("M_{#gamma#gamma} (GeV)");
    sframe->GetYaxis()->SetTitleOffset(1.35);
    sframe->Draw();
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.Draw();
    canv.Print(Form("fitresult/%s.pdf",proc.c_str()));

}

double getGoodnessOfFit(RooRealVar* mass, RooAbsPdf* bkgPDF, RooDataSet* data){

    RooRealVar norm("norm","norm",data->sumEntries(),0,10e6);
    RooExtendPdf* pdf = new RooExtendPdf("ext","ext",*bkgPDF,norm);

    RooPlot *plot_chi2 = mass->frame();
    data->plotOn(plot_chi2,Binning(320),Name("data"));
    pdf->plotOn(plot_chi2,Name("pdf"));
    int np = pdf->getParameters(*data)->getSize();
    double chi2 = plot_chi2->chiSquare("pdf","data",np);
    double prob = TMath::Prob(chi2*(320-np),320-np);

return prob;
delete pdf;

}

map<string,int> plot(RooRealVar* Mgg,map<string,int>candPDF, RooDataSet* BkgDataSet, string* bestpdf, int* bestorder){

    rootlogon();

    map<string,int> tmp;
    map<TString,TObject*> multipdfleg;
    int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
    Mgg->setRange("Range1",100.,115.);
    Mgg->setRange("Range2",135.,180.);

    RooPlot* bframe = Mgg -> frame(Title("Mgg"),Bins(80));
    BkgDataSet->plotOn(bframe,MarkerStyle(20),CutRange("Range1,Range2"));
    TLegend leg(0.6,0.7,0.9,0.92);
    TObject *pdfLeg0 = bframe->getObject(int(bframe->numItems())-1);
    leg.AddEntry(pdfLeg0,"Data","LP");
    string bestfit = "";
    int best_order = 0;
    int j = 0;
    double mincorre = 10e8;
    for(map<string,int>::iterator it = candPDF.begin();it!=candPDF.end();++it){
        int order = 1;
        while(order<=it->second){
	    
            RooAbsPdf* bkgCandPDf = BkgPdfbuild(Mgg,it->first,order);

	    if(!bkgCandPDf){ order++; continue;}
            RooFitResult* res = bkgCandPDf->fitTo(*BkgDataSet,Minimizer("Minuit2","minimize"),Save(true),Verbose(false),SumW2Error(true));
	    double prob = getGoodnessOfFit(Mgg,bkgCandPDf,BkgDataSet);
	    if(order!=it->second && prob <0.01){ order++;continue; }
	    double correction = 2*res->minNll()+order;

            if(correction<mincorre){

               mincorre = correction;
	       bestfit = it->first;
	       best_order = order;

	    }

	    cout << it->first << " : " << order << " correction  =  " << correction <<endl;
            bkgCandPDf -> plotOn(bframe,LineColor(color[j]));
	    TObject *pdfLeg = bframe->getObject(int(bframe->numItems())-1);
            tmp.insert(make_pair(it->first,order));
	    multipdfleg.insert(make_pair(Form("%s_%d",(it->first).c_str(),order),pdfLeg));
            order++;
	    j++;
	}
    }

    *bestpdf = bestfit;
    *bestorder = best_order;
    cout << "Best fit : " << bestfit << best_order << endl;
    
    for(map<TString,TObject*>::iterator it = multipdfleg.begin();it!=multipdfleg.end();++it){
	if(it->first == Form("%s_%d",bestfit.c_str(),best_order) ) leg.AddEntry(it->second,it->first + "(best of fit)","L");
	else leg.AddEntry(it->second,it->first,"L");
    }
    TCanvas canv("canv","",700,600);
    canv.cd();
    bframe->Draw();
    bframe->SetXTitle("M_{#gamma#gamma} (GeV)");
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.Draw();
    canv.Print("fitresult/bkg.pdf");
    return tmp;

}

int main(){

//Option
//-------------------------------------------------------------------------------------
    double nBDTbin = 100;
    bool _doMCStudy = false;
    double LUMI = 30.;
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
//Input Information  
//-------------------------------------------------------------------------------------
    TString inpath = "/wk_cms/youying/Analysis3/fitter/";
    TString infile = "DataSet.root";
    system("mkdir -p fitresult");
    TFile*myfile = TFile::Open( inpath + infile );
    RooWorkspace* WS = (RooWorkspace*)myfile->Get("ws");

    //Observable(Mgg,BDT)
    //---------------------------------------------------------------------------------
    RooRealVar* Mgg = (RooRealVar*)WS->var("Mgg");
    RooRealVar* BDT = (RooRealVar*)WS->var("BDT");
    //DataSet
    //---------------------------------------------------------------------------------
    RooDataSet* Hgg_vbf_126 = (RooDataSet*)WS->data("Hgg_vbf_126");
    RooDataSet* Hgg_ggf_126 = (RooDataSet*)WS->data("Hgg_ggf_126");
    RooDataSet* Hgg_vbf_126w = (RooDataSet*)WS->data("Hgg_vbf_126w");
    RooDataSet* Hgg_ggf_126w = (RooDataSet*)WS->data("Hgg_ggf_126w");
    RooDataSet* bkg  = (RooDataSet*)WS->data("bkg_13TeV");
    RooDataSet* data = (RooDataSet*)WS->data("data");

//Signal Model Construction    
//-------------------------------------------------------------------------------------
    int MH =126;
    map<string,RooDataSet*> sigdataset;
    map<string,RooDataSet*> sigdatasetw;
    sigdataset["vbf"] = (RooDataSet*)(Hgg_vbf_126->reduce(RooArgSet(*Mgg)));
    sigdataset["ggf"] = (RooDataSet*)(Hgg_ggf_126->reduce(RooArgSet(*Mgg)));
    //To get expexted yields 
    sigdatasetw["vbf"] = (RooDataSet*)(Hgg_vbf_126w->reduce(RooArgSet(*Mgg)));
    sigdatasetw["ggf"] = (RooDataSet*)(Hgg_ggf_126w->reduce(RooArgSet(*Mgg)));

    //F-Test
    //--------------------------------------------------------------------------------------
    map<string, RooAbsPdf*> SignalPDF;
    for(map<string,RooDataSet*>::iterator it=sigdataset.begin();it!=sigdataset.end();++it){

        int pre_order = 0;
        int best_order = 0;
        double thisNll = 0;
        double preNll = 0;
        double chi2 = 0;
        double prob = 0;
	vector<RooAbsPdf*> Gaussians;
    
	for(int order=1;order<7;order++){
    
            RooAbsPdf* SigMggPdf = buildMultiGaussians(Mgg,order,MH,false,it->first,it->second);
            RooFitResult* fitres = SigMggPdf->fitTo(*(it->second),Minimizer("Minuit2","minimize"),Save(true),Verbose(false),SumW2Error(true),Range(115.,135.));
    
            thisNll = fitres->minNll();
            chi2 = 2.*(preNll - thisNll);
            if(order > 1 && chi2 < 0 ) chi2 = 0;
            int delta_dof = (2*order+(order-1))-(2*pre_order+(pre_order-1));
            prob = TMath::Prob(chi2,delta_dof);
            
            cout <<"[INFO]:"<<"This order = "<< order <<","<<"Pre Order = "<<pre_order<<","<< "p value = " << prob <<endl;
	    if(prob>0.05) break;
	    Gaussians.push_back(SigMggPdf);
            pre_order = order;
            preNll = thisNll;
	}
      
        best_order = pre_order;
        //best_order = 1;

        makeSigPlot(Mgg,Gaussians,it->first,it->second);
    
        RooAbsPdf* bestSigPDF = buildMultiGaussians(Mgg,best_order,MH,true/*Do fit*/,it->first,it->second);
        SignalPDF.insert(make_pair(it->first,bestSigPDF));
    }
    map<string, RooExtendPdf*> extSignalPDF;
    for(map<string,RooAbsPdf*>::iterator it = SignalPDF.begin();it!=SignalPDF.end();++it){
        RooRealVar* norm = new RooRealVar("norm","",sigdatasetw[it->first]->sumEntries(),-1000.,1000.);
        norm->setConstant(true);
        RooExtendPdf* extsigpdf = new RooExtendPdf(Form("%s_Extend",(it->first).c_str()),"",*(it->second),*norm);
        extSignalPDF.insert(make_pair(it->first,extsigpdf));
    }

//Background Shapes Construction
//-----------------------------------------------------------------------------------------------------------------------------------
    vector<string>PDFfamily;
    PDFfamily.push_back("Bernstein");
    PDFfamily.push_back("Laurent");
    PDFfamily.push_back("PowerLaw");
    //PDFfamily.push_back("PowerLaw2");
    PDFfamily.push_back("Exponential");
    //PDFfamily.push_back("Exponential2");

    //RooDataSet* BkgDataSet = (RooDataSet*)bkg->reduce(RooArgSet(*Mgg));

    //Use the dataSet of real data
    RooDataSet* BkgDataSet = (RooDataSet*)data->reduce(RooArgSet(*Mgg));

    //F-Test
    //---------------------------------------------------------------------------------------------------------------------------
    map<string,int>candPDF;
    for(vector<string>::iterator it = PDFfamily.begin();it!=PDFfamily.end();++it){

        int order = 1;
        int pre_order = 0;
        int cache_order = 0;
        int best_order = 0;
        double thisNll = 0.;
        double preNll = 10e8;
        double chi2 = 0.;
        double prob = 0.;
        int sts = 0;
        while(order<7 && prob < 0.05  ){

              RooAbsPdf* bkgCandPDf = BkgPdfbuild(Mgg,*it,order);
	    
              if(!bkgCandPDf){
		 order++;
              }else{
                 RooFitResult* fitres = bkgCandPDf->fitTo(*BkgDataSet,Minimizer("Minuit2","minimize"),Save(true),Verbose(false),SumW2Error(true));
                 thisNll = fitres->minNll();
                 chi2 = 2.*(preNll - thisNll);
                 if(order > 1 && chi2 < 0 ) chi2 = 0;
                 int delta_dof = order - pre_order;
                 prob = TMath::Prob(chi2,delta_dof);
	         sts = fitres ->status();

                 cout <<"[INFO]:"<<"This order = "<< order <<","<<"Pre Order = "<<pre_order<<", Type = "<<*it <<",p value = " << prob <<"  state = " << sts <<endl;
                 cache_order = pre_order;
                 pre_order = order;
                 preNll = thisNll;
                 order++;
	      }
	}
        best_order = cache_order;
        candPDF.insert(pair<string,int>(*it,best_order));

    }
    string bestPDF = "";
    int best_order = 0;
    map<string,int> multipdf = plot(Mgg,candPDF,BkgDataSet,&bestPDF,&best_order);
    RooAbsPdf* MggBkgShape = BkgPdfbuild(Mgg,bestPDF,best_order,true/*Do fit*/,BkgDataSet);
    envelope(Mgg,multipdf,BkgDataSet,bestPDF,best_order,extSignalPDF,true/*blind*/,true/*addSignal*/,false/*combine*/);

//BDT Shapes Test
//-------------------------------------------------------------------------------------------------------------------------------------------------------
    RooDataSet* vbfBDT = (RooDataSet*)Hgg_vbf_126->reduce(RooArgSet(*BDT));
    RooDataSet* ggfBDT = (RooDataSet*)Hgg_ggf_126->reduce(RooArgSet(*BDT));
    RooDataSet* bkgBDT = (RooDataSet*)bkg->reduce(RooArgSet(*BDT));
    RooDataSet* dataBDT = (RooDataSet*)data->reduce(RooArgSet(*BDT));
    RooAbsPdf* vbfBDTPDF;
    RooAbsPdf* ggfBDTPDF;
    RooAbsPdf* bkgBDTPDF;
    bool isBinnedBDT = false;
    if(isBinnedBDT){
       BDT->setBins(nBDTbin);
       RooDataHist*DHvbfBDT = vbfBDT->binnedClone();
       RooDataHist*DHggfBDT = ggfBDT->binnedClone();
       RooDataHist*DHbkgBDT = bkgBDT->binnedClone();
       RooHistPdf* HvbfBDTPDF = new RooHistPdf("HvbfBDTPDF","",RooArgSet(*BDT),*DHvbfBDT);
       RooHistPdf* HggfBDTPDF = new RooHistPdf("HvbfBDTPDF","",RooArgSet(*BDT),*DHggfBDT);
       RooHistPdf* HbkgBDTPDF = new RooHistPdf("HvbfBDTPDF","",RooArgSet(*BDT),*DHbkgBDT);
       vbfBDTPDF = (RooAbsPdf*)HvbfBDTPDF->Clone();
       ggfBDTPDF = (RooAbsPdf*)HggfBDTPDF->Clone();
       bkgBDTPDF = (RooAbsPdf*)HbkgBDTPDF->Clone();
    }
  
    //VBF Unbinned Shape
    //------------------------------------------------------------------------------------------------------
    RooRealVar* BDTmean = new RooRealVar("BDTmean","BDTmean",1.,0,10);	
    RooRealVar* sigma   = new RooRealVar("sigma","sigma"    ,0.2,0,10);	
    RooRealVar* alpha   = new RooRealVar("alpha","alpha"    ,1.,0,10);	
    RooRealVar* n       = new RooRealVar("n","n"            ,1.3,0,10);	
    vbfBDTPDF = new RooCBShape("vbfBDTPDF","vbfBDTPDF",*BDT,*BDTmean,*sigma,*alpha,*n);

    vbfBDTPDF->fitTo(*vbfBDT,Minimizer("Minuit2","minimize"),SumW2Error(true)); 
    BDTmean -> setConstant(true);
    sigma   -> setConstant(true);
    alpha   -> setConstant(true);
    n       -> setConstant(true);

    RooPlot* VBFBDTframe = BDT->frame(Bins(nBDTbin));
    vbfBDT->plotOn(VBFBDTframe,DataError(RooAbsData::SumW2));
    vbfBDTPDF->plotOn(VBFBDTframe);
    TCanvas canv1("c1","c1",650,600);
    canv1.cd();VBFBDTframe->Draw();canv1.Print("fitresult/VBFBDT.pdf");

    //GGF Unbinned Shape
    //------------------------------------------------------------------------------------------------------
    RooRealVar* gBDTmean = new RooRealVar("gBDTmean","gBDTmean",1,0,10);
    RooRealVar* gsigma   = new RooRealVar("gsigma"  ,"gsigma"  ,0.2,0,10);
    RooRealVar* galpha   = new RooRealVar("galpha"  ,"galpha"  ,1,0,10);
    RooRealVar* gn       = new RooRealVar("gn"      ,"gn"      ,1.3,0,10);
    ggfBDTPDF = new RooCBShape("ggfBDTPDF","ggfBDTPDF",*BDT,*gBDTmean,*gsigma,*galpha,*gn);

    ggfBDTPDF->fitTo(*ggfBDT,Minimizer("Minuit2","minimize2"),SumW2Error(true));
    gBDTmean->setConstant(true);
    gsigma  ->setConstant(true);
    galpha  ->setConstant(true);
    gn      ->setConstant(true);

    RooPlot* GGFBDTframe = BDT->frame(Bins(nBDTbin));
    ggfBDT->plotOn(GGFBDTframe,DataError(RooAbsData::SumW2));
    ggfBDTPDF->plotOn(GGFBDTframe);
    TCanvas canv2("c2","c2",650,600);
    canv2.cd();GGFBDTframe->Draw();canv2.Print("fitresult/GGFBDT.pdf");

    //Background Unbinned Shape
    //------------------------------------------------------------------------------------------------------
    bkgBDTPDF = BkgPdfbuild(BDT,"Bernstein",3,true,dataBDT);
    RooPlot* bkgBDTframe = BDT->frame(Title("BDT"),Bins(nBDTbin));
    dataBDT->plotOn(bkgBDTframe,DataError(RooAbsData::SumW2));
    bkgBDTPDF->plotOn(bkgBDTframe);
    TCanvas canv3("canv3","canv3",650,600);        
    canv3.cd();bkgBDTframe->Draw();canv3.Print("fitresult/ContBDT.pdf");

//Conbine Shapes(1D * 1D shapes)
//---------------------------------------------------------------------------------------------------------------------------------------
    //Calculate Scale Factor
    RooDataSet* sMCDataSet = (RooDataSet*)bkg->reduce(RooArgSet(*Mgg),"!(Mgg>120&&Mgg<130)");
    RooDataSet* sBkgDataSet = (RooDataSet*)data->reduce(RooArgSet(*Mgg),"!(Mgg>120&&Mgg<130)");
    double scale = sBkgDataSet->sumEntries()/sMCDataSet->sumEntries();

    //VBF 1D*1D
    double SM_VBF = sigdatasetw["vbf"]->sumEntries() /12.887*LUMI;
    RooRealVar* N_vbf = new RooRealVar("N_vbf","N_vbf", SM_VBF, -1000, 1000);
    RooProdPdf*h2dVBFPDF = new RooProdPdf("h2dVBFPDF","Mgg*BDT",*(SignalPDF["vbf"]),*vbfBDTPDF);
    //GGF 1D*1D
    double SM_GGF = sigdatasetw["ggf"]->sumEntries() /12.887*LUMI;
    RooRealVar* N_ggf = new RooRealVar("N_ggf","N_ggf", SM_GGF, -1000, 1000);
    RooProdPdf* h2dGGFPDF = new RooProdPdf("h2dGGFPDF","Mgg*BDT",*(SignalPDF["ggf"]),*ggfBDTPDF);
    //Bkg 1D*1D
    RooDataSet* MCDataSet = (RooDataSet*)bkg->reduce(RooArgSet(*Mgg));
    double SM_Bkg = MCDataSet->sumEntries() /12.887*LUMI*scale;
    RooRealVar* N_Bkg = new RooRealVar("N_gj","N_gj", SM_Bkg, 0, 20000);
    RooProdPdf* h2dBkgPDF = new RooProdPdf("h2dBkgPDF","Mgg*BDT",*MggBkgShape,*bkgBDTPDF);

    cout << SM_VBF <<"  "<< SM_GGF << "  " << BkgDataSet->sumEntries() <<endl;
    //Store the shapes into the same workspace file
    RooWorkspace outws("outws");
    outws.import(*h2dVBFPDF);
    outws.import(*h2dGGFPDF);
    outws.import(*h2dBkgPDF);
    outws.import(*data);
    outws.writeToFile("FinalResult.root",true);
    outws.Print();



//Make Plots about the Results of 1D*1D Fitting
//----------------------------------------------------------------------------------------------------------------------------------------

    rootlogon();

    //Ext PDF
    //------------------------------------------------------------------------------------------------------------
    RooExtendPdf* exth2dVBFPDF = new RooExtendPdf("exh2dVBFPDF","",*h2dVBFPDF,*N_vbf);
    RooExtendPdf* exth2dGGFPDF = new RooExtendPdf("exh2dGGFPDF","",*h2dGGFPDF,*N_ggf);
    RooExtendPdf* exth2dBkgPDF = new RooExtendPdf("exh2dBkgPDF","",*h2dBkgPDF,*N_Bkg);

    //Totally Conbime
    //------------------------------------------------------------------------------------------------------------
    RooAddPdf* model = new RooAddPdf("model","model",RooArgList(*exth2dVBFPDF,*exth2dGGFPDF,*exth2dBkgPDF));

    //Total 1D*1D PDF Plot
    //-------------------------------------------------------------------------------------------
    TCanvas canv4("canv4","",600,600);
    canv4.cd();	
    TH1* h2model = model->createHistogram("Mgg,BDT",80,100);
    h2model->SetLineColor(kBlue);
    h2model->Draw("surf");
    //h2model->Draw("histo");
    canv4.SetPhi(220);
    canv4.Print("fitresult/Fit2d.pdf");

    //Decouple to Mgg and BDT
    //-------------------------------------------------------------------------------------------
    int gen_num = SM_VBF + SM_GGF + SM_Bkg;
    RooDataSet* fakedata = model->generate(RooArgSet(*Mgg,*BDT),gen_num);

    TCanvas canv5("canv5","",1300,600);
    canv5.Divide(2,1);
    //Mgg
    //-------------------------------------
    canv5.cd(1);
    RooPlot* Mggframe = Mgg->frame(Bins(80),Range(100,180));
    fakedata->plotOn(Mggframe);
    TObject* dataleg = (TObject*)Mggframe->getObject(Mggframe->numItems()-1);
    model->plotOn(Mggframe) ;
    TObject* totalleg = (TObject*)Mggframe->getObject(Mggframe->numItems()-1);
    model->plotOn(Mggframe,Components(*exth2dVBFPDF),LineColor(kGreen),LineStyle(kDashed));
    TObject* vbfleg = (TObject*)Mggframe->getObject(Mggframe->numItems()-1);
    model->plotOn(Mggframe,Components(*exth2dGGFPDF),LineColor(kRed),LineStyle(kDashed));
    TObject* ggfleg = (TObject*)Mggframe->getObject(Mggframe->numItems()-1);
    model->plotOn(Mggframe,Components(*exth2dBkgPDF),LineColor(kBlack),LineStyle(kDashed));
    TObject* bkgleg = (TObject*)Mggframe->getObject(Mggframe->numItems()-1);
    Mggframe->SetXTitle("M_{#gamma#gamma} (GeV)");
    Mggframe->Draw();
    TLegend* leg1 = new TLegend(0.65, 0.75, 0.88, 0.93);
    leg1->AddEntry(dataleg,Form("pseudo data %.1f fb^{-1}",LUMI),"LP");
    leg1->AddEntry(totalleg,"MggPDF","L");
    leg1->AddEntry(vbfleg,"VBF","L");
    leg1->AddEntry(ggfleg,"GGF","L");
    leg1->AddEntry(bkgleg,"Continuum","L");
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->Draw(); 

    //BDT
    //-----------------------------------------------------------------------------------------------------------------
    canv5.cd(2);
    RooPlot* BDTframe = BDT->frame(Bins(100)) ;
    fakedata->plotOn(BDTframe);
    TObject* datalegb = (TObject*)BDTframe->getObject(BDTframe->numItems()-1);
    model->plotOn(BDTframe) ;
    TObject* totallegb = (TObject*)BDTframe->getObject(BDTframe->numItems()-1);
    model->plotOn(BDTframe,Components(*exth2dVBFPDF),LineColor(kGreen),LineStyle(kDashed),RooFit::Normalization(10));
    TObject* vbflegb = (TObject*)BDTframe->getObject(BDTframe->numItems()-1);
    model->plotOn(BDTframe,Components(*exth2dGGFPDF),LineColor(kRed),LineStyle(kDashed),RooFit::Normalization(10));
    TObject* ggflegb = (TObject*)BDTframe->getObject(BDTframe->numItems()-1);
    model->plotOn(BDTframe,Components(*exth2dBkgPDF),LineColor(kBlack),LineStyle(kDashed));
    TObject* bkglegb = (TObject*)BDTframe->getObject(BDTframe->numItems()-1);
    BDTframe->SetXTitle("dijet BDT");
    BDTframe->Draw();
    TLegend* leg2 = new TLegend(0.65, 0.75, 0.88, 0.93);
    leg2->AddEntry(datalegb,Form("pseudo data %.1f fb^{-1}",LUMI),"LP");
    leg2->AddEntry(totallegb,"BDTPDF","L");
    leg2->AddEntry(vbflegb,"VBF x10","L");
    leg2->AddEntry(ggflegb,"GGF x10","L");
    leg2->AddEntry(bkglegb,"Continuum","L");
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->Draw();
    canv5.Print("fitresult/BDT.pdf");


//MC Study
//--------------------------------------------------------------------------------------------------------------------------------------------

    if(_doMCStudy){
       RooMCStudy* mcstudy = new RooMCStudy(*model,RooArgSet(*Mgg,*BDT),Binned(false),Silence(),Extended(),FitOptions(Save(kTRUE),SumW2Error(true),Verbose(kFALSE),PrintEvalErrors(0)));
       //RooMCStudy* mcstudy = new RooMCStudy(*model,RooArgSet(*Mgg,*BDT),Binned(false),Silence(),Extended(),FitOptions(Save(kTRUE),Verbose(kFALSE),PrintEvalErrors(0)));
     
       mcstudy->generateAndFit(1000);
  
       //RooPlot* MCframe1 = mcstudy->plotParam(*N_vbf,Bins(100)) ;
       //RooPlot* MCframe2 = mcstudy->plotError(*N_vbf,Bins(100)) ;
       RooPlot* MCframe3 = mcstudy->plotPull(*N_vbf,Bins(50) ,FitGauss(kTRUE)) ;
       //RooPlot* MCframe4 = mcstudy->plotParam(*N_ggf,Bins(100)) ;
       //RooPlot* MCframe5 = mcstudy->plotError(*N_ggf,Bins(100)) ;
       RooPlot* MCframe6 = mcstudy->plotPull(*N_ggf,Bins(50) ,FitGauss(kTRUE)) ;
  
       TCanvas canv6("canv6","",1300,600);
       canv6.Divide(2,1);
       //canv6->cd();MCframe1->Draw();
       //canv6->cd();MCframe2->Draw();
       canv6.cd(1);MCframe3->Draw();
       //canv6->cd();MCframe4->Draw();
       //canv6->cd();MCframe5->Draw();
       canv6.cd(2);MCframe6->Draw();
  
       canv6.Print("fitresult/pull.pdf");
   }

return 0;
}
