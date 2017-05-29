#ifndef __ENVELOPE__
#define __ENVELOPE__


#include <iostream>
#include <map>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooMinimizer.h"
#include "RooHist.h"

#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TH1.h"
#include "TLatex.h"

#include "PDFLIB.h"

#define LUMIS 12.887

using namespace std;

typedef map<string,int> multiPDF_t;

double getNormNll(RooRealVar*mass, multiPDF_t multipdf, RooDataSet* data, double norm, double masslow, double masshigh){

    double bestFitNll=1.e9;

    for(multiPDF_t::iterator it = multipdf.begin();it!=multipdf.end();++it){
        mass->setRange("errRange",masslow,masshigh);
	RooRealVar normVar(Form("%s_norm",(it->first).c_str()),"",0.,1.e6);
	RooAbsPdf* bkgpdf = BkgPdfbuild(mass,it->first,it->second);
	RooExtendPdf extPdf(Form("%s_extPdf",(it->first).c_str()),"",*bkgpdf,normVar,"errRange");
	RooAbsReal* nll = extPdf.createNLL(*data,Extended());
	
        normVar.setConstant(false);
	normVar.setVal(norm);
	normVar.setConstant(true);

        RooMinimizer minim(*nll);
	minim.setStrategy(0);
	minim.migrad();
	double corrNll = 2*(nll->getVal())+it->second;
	if(corrNll<bestFitNll) bestFitNll = corrNll;
    }
return bestFitNll;
}

double bandvalue(RooRealVar* mass, multiPDF_t multipdf, RooDataSet* data, double bestPoint ,double nllbest, double boundary, double masslow, double masshigh, double diff){

    bool isDownSide;
    double UpEdge, DownEdge, value, valueNll;
    if(boundary<bestPoint){
       isDownSide = true;
       UpEdge = bestPoint;
       DownEdge = boundary;
    }else{
       isDownSide = false;
       UpEdge = boundary;
       DownEdge = bestPoint;
    }

    double distancetotrue = 1.e6;
    int nits=0;

    while(TMath::Abs(distancetotrue/diff)>0.05){

          value = (UpEdge-DownEdge)/2 + DownEdge;
	  valueNll = getNormNll(mass,multipdf,data,value,masslow,masshigh) - nllbest;
	  distancetotrue = diff - valueNll;

          if(isDownSide){
	     if(valueNll>diff)  DownEdge = value;
	     else UpEdge = value;
	  }else{
             if(valueNll>diff)  UpEdge = value;
	     else DownEdge = value;
	  }
	  nits++;

	  if(nits>20){
             return value;
	     DownEdge = TMath::Max(0.,DownEdge-20.);
	     UpEdge += 20.;
	     nits = 0;
	     if(TMath::Abs(valueNll)>2.e4) return 0;
	  }
    }
    return value;
}



void envelope(RooRealVar* mass, multiPDF_t multipdf, RooDataSet* data, string best_PDF, int best_order, map<string,RooExtendPdf*> sigPDF, bool isblind = false, bool addsignal = true, bool combine =true){
     
     
     RooAbsPdf* bestPDF = BkgPdfbuild(mass,best_PDF,best_order,true,data);
     RooPlot* plot = mass->frame();
     plot-> GetXaxis()->SetLabelSize(0.);
     data->plotOn(plot,Binning(80),Invisible());
     RooHist* plotdata = (RooHist*)plot->getObject(plot->numItems()-1);
     TObject* dataleg = (TObject*)plot->getObject(plot->numItems()-1);
     if(combine) bestPDF->plotOn(plot,LineColor(kRed),LineWidth(3),LineStyle(7));
     else bestPDF->plotOn(plot,LineColor(kRed),LineWidth(3));     
     RooCurve* bestBkgCurve = (RooCurve*)plot->getObject(plot->numItems()-1);

     TLegend leg(0.6,0.6,0.89,0.89);
     leg.AddEntry(dataleg,"data","LEP");
     leg.AddEntry(bestBkgCurve,"Bkg fit","L");
     leg.SetFillStyle(0);
     leg.SetBorderSize(0);

     TGraphAsymmErrors* oneSigmaband = new TGraphAsymmErrors();   //SetName
     TGraphAsymmErrors* twoSigmaband = new TGraphAsymmErrors();
     TGraphAsymmErrors* oneSigmaband_r = new TGraphAsymmErrors();
     TGraphAsymmErrors* twoSigmaband_r = new TGraphAsymmErrors();

     int p=0;
     //for(double MH = 100.;MH<180.5;MH+=0.5){
     for(double MH = 100.;MH<180.5;MH+=1.){

         double lowedge = MH-0.5;
	 double highedge = MH+0.5;
	 double center = MH;
         double normbkg = bestBkgCurve->interpolate(center);
         double bestnll = getNormNll(mass,multipdf,data,normbkg,lowedge,highedge);
         double lowRange = TMath::Max(0.,normbkg - 3*TMath::Sqrt(normbkg));
         double highRange = normbkg + 3*TMath::Sqrt(normbkg);

         double errdown1Value = bandvalue(mass,multipdf,data,normbkg,bestnll,lowRange,lowedge,highedge,1);
         double errdown2Value = bandvalue(mass,multipdf,data,normbkg,bestnll,lowRange,lowedge,highedge,4);
         double errUp1Value = bandvalue(mass,multipdf,data,normbkg,bestnll,highRange,lowedge,highedge,1);
         double errUp2Value = bandvalue(mass,multipdf,data,normbkg,bestnll,highRange,lowedge,highedge,4);

         double errdown1 = normbkg - errdown1Value;
         double errdown2 = normbkg - errdown2Value;
         double errUp1 = errUp1Value - normbkg;
         double errUp2 = errUp2Value - normbkg;


         oneSigmaband->SetPoint(p,center,normbkg);
         twoSigmaband->SetPoint(p,center,normbkg);
	 oneSigmaband->SetPointError(p,0.,0.,errdown1,errUp1);
	 twoSigmaband->SetPointError(p,0.,0.,errdown2,errUp2);

         oneSigmaband_r->SetPoint(p,center,0.);
         twoSigmaband_r->SetPoint(p,center,0.);
         oneSigmaband_r->SetPointError(p,0.,0.,errdown1,errUp1);
         twoSigmaband_r->SetPointError(p,0.,0.,errdown2,errUp2);


	 cout << " step = " << p << ", MH = " << MH <<endl;
         p++;
     }

     TCanvas canv("canv","",600,650);
     canv.cd();
     TPad pad1("pad1","",0.,0.24,1.,0.98);
     pad1.SetTopMargin(0.05);
     pad1.SetBottomMargin(0.019);
     pad1.Draw();
     pad1.cd();
     plot->Draw();

     twoSigmaband->SetLineColor(kYellow);
     twoSigmaband->SetFillColor(kYellow);
     twoSigmaband->SetMarkerColor(kYellow);
     twoSigmaband->Draw("L3 same");
     oneSigmaband->SetLineColor(kGreen);
     oneSigmaband->SetFillColor(kGreen);
     oneSigmaband->SetMarkerColor(kGreen);
     oneSigmaband->Draw("L3 same");
     leg.AddEntry(oneSigmaband,"#pm1#sigma","F");
     leg.AddEntry(twoSigmaband,"#pm2#sigma","F");

     twoSigmaband_r->SetLineColor(kYellow);
     twoSigmaband_r->SetFillColor(kYellow);
     twoSigmaband_r->SetMarkerColor(kYellow);
     oneSigmaband_r->SetLineColor(kGreen);
     oneSigmaband_r->SetFillColor(kGreen);
     oneSigmaband_r->SetMarkerColor(kGreen);

     RooPlot*splot = mass-> frame();

     if(addsignal){
        sigPDF["ggf"]->plotOn(plot,Normalization(1.0,RooAbsReal::RelativeExpected),LineColor(kBlue),LineWidth(1));
        RooCurve* ggfCurve = (RooCurve*)plot->getObject(plot->numItems()-1);
        sigPDF["ggf"]->plotOn(plot,Normalization(1.0,RooAbsReal::RelativeExpected),LineColor(kBlue),LineWidth(1),FillColor(38),FillStyle(3001),DrawOption("F"));
        TObject* ggfleg = (TObject*)plot->getObject(plot->numItems()-1);
        sigPDF["vbf"]->plotOn(plot,Normalization(1.0,RooAbsReal::RelativeExpected),LineColor(kRed),LineWidth(1));
        RooCurve* vbfCurve = (RooCurve*)plot->getObject(plot->numItems()-1);
        sigPDF["vbf"]->plotOn(plot,Normalization(1.0,RooAbsReal::RelativeExpected),LineColor(kRed),LineWidth(1),FillColor(46),FillStyle(3001),DrawOption("F"));
        TObject* vbfleg = (TObject*)plot->getObject(plot->numItems()-1);
	leg.AddEntry(vbfleg,"vbf(M_{H} = 126GeV)","F");
	leg.AddEntry(ggfleg,"ggf(M_{H} = 126GeV)","F");
	RooRealVar fvbf("fvbf","",vbfCurve->average(100.,180.)*80.,0.,1000.); fvbf.setConstant(true);
	RooRealVar fggf("fggf","",ggfCurve->average(100.,180.)*80.,0.,1000.); fggf.setConstant(true);
        RooAddPdf sigTotalPdf("total","",RooArgList(*sigPDF["ggf"],*sigPDF["vbf"]),RooArgList(fggf,fvbf));
        sigTotalPdf.plotOn(plot,Normalization(1.0,RooAbsReal::RelativeExpected),LineColor(kViolet+2),LineWidth(4),LineStyle(2));
        TObject* totalleg = (TObject*)plot->getObject(plot->numItems()-1);
	leg.AddEntry(totalleg,"vbf+ggf(M_{H} = 126GeV}","L");
        if(combine){
           sigTotalPdf.plotOn(splot,Normalization(1.0,RooAbsReal::RelativeExpected),LineColor(kRed),LineWidth(3),LineStyle(1));
	   RooRealVar fsig("fsig","",fvbf.getValV()+fggf.getValV(),0.,100000.); fsig.setConstant(true);
	   RooRealVar fbkg("fbkg","",bestBkgCurve->average(100.,180.)*80.,0.,100000.); fbkg.setConstant(true);
	   RooAddPdf TotalPDF("TotalPDF","",RooArgList(*bestPDF,sigTotalPdf),RooArgList(fbkg,fsig));
           TotalPDF.plotOn(plot,Normalization(1.0,RooAbsReal::RelativeExpected),LineColor(kRed),LineWidth(3),LineStyle(1));
	}
     }

     if(isblind){
        mass->setRange("unblind_low",100.,115.);
        mass->setRange("unblind_high",135.,180.);
        data->plotOn(plot,Binning(80),CutRange("unblind_low,unblind_high"),MarkerSize(0.9),XErrorSize(0));
     }else{
        data->plotOn(plot,Binning(80),MarkerSize(0.9));
     }

     plot->Draw("same");
     leg.Draw("same");
     canv.Update();

     TLatex* latex = new TLatex();
     latex -> SetNDC();
     latex -> SetTextFont(62);
     latex -> SetTextSize(0.055);
     latex -> DrawText(0.12,0.96,"CMS");
     latex -> SetTextFont(42);
     latex -> SetTextSize(0.045);
     latex -> DrawText(0.23,0.96,"Preliminary");
     latex -> SetTextSize(0.04);
     latex -> DrawLatex(0.6,0.96,Form("#sqrt{s}=13 TeV, L=%.3f fb^{-1}",LUMIS));
     latex -> SetTextSize(0.035);
     latex -> DrawLatex(0.3,0.88,"CMSSW 80X");
     latex -> DrawLatex(0.3,0.835,"Dijet BDT > 0.4");

     canv.cd();
     TPad pad2("pad2","",0.,0.01,1.,0.23);
     pad2.SetTopMargin(0.005);
     pad2.SetBottomMargin(0.35);
     pad2.Draw();
     pad2.cd();

     int npoint = plotdata->GetN();
     double xtmp,ytmp;

     TGraphAsymmErrors *hdatasub = new TGraphAsymmErrors(npoint);

     for(int i=0;i<npoint;++i){

         plotdata->GetPoint(i,xtmp,ytmp);
         double bkgvalue = bestBkgCurve->interpolate(xtmp);
	 if(isblind){
            if(xtmp>115.&&xtmp<135.) continue;
	 }
         double errorUp = plotdata -> GetErrorYhigh(i);
         double errorDown = plotdata -> GetErrorYlow(i);
	 hdatasub->SetPoint(i,xtmp,ytmp-bkgvalue);
	 hdatasub->SetPointError(i,0.,0.,errorDown,errorUp);

     }

     TH1 *hdata = new TH1D("hdata","",80,100,180);
     hdata->SetMaximum(hdatasub->GetHistogram()->GetMaximum()+1);
     hdata->SetMinimum(hdatasub->GetHistogram()->GetMinimum()-1);
     //hdata->GetYaxis()->SetTitle("data - best fit PDF");
     //hdata->GetYaxis()->SetTitleSize(0.12);
     hdata->GetYaxis()->SetLabelSize(0.14);
     hdata->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
     hdata->GetXaxis()->SetTitleSize(0.14);
     hdata->GetXaxis()->SetLabelSize(0.15);
     hdata->Draw("HIST");
     twoSigmaband_r->Draw("L3 SAME");
     oneSigmaband_r->Draw("L3 SAME");
     hdata->GetYaxis()->SetNdivisions(808);

     TLatex* latex2 = new TLatex();
     latex2 -> SetNDC();
     latex2 -> SetTextFont(42);
     latex2 -> SetTextSize(0.13);
     latex2 -> DrawText(0.7,0.86,"Bkg subtracted");
		  
     TLine *line3 = new TLine(100,0.,180,0.);
     line3->SetLineColor(kRed);
     line3->SetLineWidth(3);
     line3->Draw();
     if(combine){
	line3->SetLineStyle(7);
	splot->Draw("SAME");
     }
     hdatasub->Draw("PESAME");
     hdatasub->SetMarkerStyle(20);
     hdatasub->SetMarkerSize(0.9);


     canv.Update();

     canv.Print("fitresult/envelope.pdf");

}

#endif
