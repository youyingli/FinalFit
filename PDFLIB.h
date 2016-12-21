#ifndef __PDFLIB__
#define __PDFLIB__


#include <iostream>
#include <cmath>

#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooBernstein.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooAbsPdf.h"

using namespace std;
using namespace RooFit;


RooAbsPdf* buildMultiGaussians(RooRealVar* mass, int nGaussians, int mh, bool dofit, string signame ,RooDataSet* Mggdataset ){

    RooArgList* gaussians = new RooArgList();
    RooArgList* coeffs = new RooArgList();
    RooRealVar* MH = new RooRealVar("MH","MH",(double)mh,115.,135.);
    MH->setConstant(true);

    RooRealVar* dm[nGaussians];
    RooAbsReal* mean[nGaussians];
    RooRealVar* sigma[nGaussians];
    RooGaussian* gaus[nGaussians];
    RooRealVar* frac[nGaussians-1];
    for(int i=0;i<nGaussians;i++){
        dm[i]    = new RooRealVar(Form("%s_dm_mh%d_g%d",signame.c_str(),mh,i),Form("%s_dm_mh%d_g%d",signame.c_str(),mh,i),0.1,-8.,8.);
	mean[i]  = new RooFormulaVar(Form("%s_mean_mh%d_g%d",signame.c_str(),mh,i),Form("%s_mean_mh%d_g%d",signame.c_str(),mh,i),"@0+@1",RooArgList(*MH,*dm[i]));
	sigma[i] = new RooRealVar(Form("%s_sigma_mh%d_g%d",signame.c_str(),mh,i),Form("%s_sigma_mh%d_g%d",signame.c_str(),mh,i),2.,0.4,20.);
	gaus[i]  = new RooGaussian(Form("%s_gaus_mh%d_g%d",signame.c_str(),mh,i),Form("%s_gaus_mh%d_g%d",signame.c_str(),mh,i),*mass,*mean[i],*sigma[i]);
	gaussians->add(*gaus[i]);

	if(i<nGaussians-1){
           frac[i] = new RooRealVar(Form("%s_frac_mh%d_g%d",signame.c_str(),mh,i),Form("%s_frac_mh%d_g%d",signame.c_str(),mh,i),0.1,0.01,0.99);
	   coeffs->add(*frac[i]);
	}
    }
    //RooAddPdf* tmpMultiGaussians = new RooAddPdf(Form("PDF_%s%d",signame.c_str(),mh),Form("PDF_%s%d",signame.c_str(),mh),*gaussians,*coeffs,kTRUE);
    RooAddPdf* tmpMultiGaussians = new RooAddPdf(Form("PDF_%s%d",signame.c_str(),mh),Form("PDF_%s%d",signame.c_str(),mh),*gaussians,*coeffs,kTRUE);
    if(dofit){ 
       tmpMultiGaussians->fitTo(*Mggdataset,Minimizer("Minuit2","minimize"),Verbose(false),SumW2Error(true));
       for(int i=0;i<nGaussians;i++){
	   dm[i]->setConstant(true); 
	   sigma[i]->setConstant(true); 
	   if(i<nGaussians-1) frac[i]->setConstant(true);
       }
    }
    return tmpMultiGaussians;
    //delete tmpMultiGaussians;

}

RooAbsPdf* BkgPdfbuild(RooRealVar* mass, string pdfType, int order, bool dofit = false, RooDataSet* Mggdataset = NULL){

    if(pdfType == "Bernstein"){
       RooRealVar* param[order];
       //RooFormulaVar* form[order];
       RooArgList* coeffs = new RooArgList();
       for(int i=0;i<=order;i++){
           //param[i] = new RooRealVar(Form("param_%d",i),Form("param_%d",i),0.1*(i+1),-5.,5.);
           param[i] = new RooRealVar(Form("param_%d",i),Form("param_%d",i),0.1*(i+1),0.,25.);
           //form[i] = new RooFormulaVar(Form("form_%d",i),Form("form_%d",i),"@0*@0",RooArgList(*param[i]));
	   //coeffs->add(*form[i]);
	   coeffs->add(*param[i]);
       }
       RooBernstein* BernPdf = new RooBernstein(Form("%s_%d",pdfType.c_str(),order),Form("%s_%d",pdfType.c_str(),order),*mass,*coeffs);
       if(dofit){
          BernPdf->fitTo(*Mggdataset,Minimizer("Minuit2","minimize"),Verbose(false),SumW2Error(true));
	  for(int i=0;i<order;i++) param[i]->setConstant(true);
       }
       return BernPdf;
    }


    if(pdfType == "Laurent"){
       int nlower = int(ceil(order/2.));
       int nhigher = order - nlower;
       RooArgList* coeffs = new RooArgList();
       RooArgList* LaurPowList = new RooArgList();
       //first 0th order
       RooGenericPdf* LaurPow0Pdf = new RooGenericPdf("LaurPow_0","LaurPow_0","TMath::Power(@0,-4.)",RooArgList(*mass));
       LaurPowList->add(*LaurPow0Pdf);

       RooRealVar* efrac[nlower];
       RooRealVar* ofrac[nhigher];
       //even terms	
       for(int i=1;i<=nlower;i++){
	   RooGenericPdf* LaurPowePdf = new RooGenericPdf(Form("LaurPow_e%d",i),Form("LaurPow_e%d",i),Form("TMath::Power(@0,-4.-%d)",i),RooArgList(*mass));
	   efrac[i] = new RooRealVar(Form("efrac_%d",i),Form("efrac_%d",i),0.25/order,0.000001,0.99999);
           LaurPowList->add(*LaurPowePdf);
	   coeffs->add(*efrac[i]);
       }
       //odd terms
       for(int i=1;i<=nhigher;i++){
	   RooGenericPdf* LaurPowoPdf = new RooGenericPdf(Form("LaurPow_o%d",i),Form("LaurPow_o%d",i),Form("TMath::Power(@0,-4.+%d)",i),RooArgList(*mass));
	   ofrac[i] = new RooRealVar(Form("ofrac_%d",i),Form("ofrac_%d",i),0.25/order,0.000001,0.99999);
           LaurPowList->add(*LaurPowoPdf);
	   coeffs->add(*ofrac[i]);
       }
       RooAddPdf* LaurPdf = new RooAddPdf(Form("%s_%d",pdfType.c_str(),order),Form("%s_%d",pdfType.c_str(),order),*LaurPowList,*coeffs,true);
       if(dofit){
          LaurPdf->fitTo(*Mggdataset,Minimizer("Minuit2","minimize"),Verbose(false),SumW2Error(true));
          for(int i=1;i<=nlower;i++) efrac[i]->setConstant(true);
          for(int i=1;i<=nhigher;i++) ofrac[i]->setConstant(true);

       }

       return LaurPdf;
    }

    if(pdfType == "PowerLaw"){
       if(order%2 == 0){
       //if(order%2 == 1){
	  //cout << "[INFO] : order has to be odd number" << endl;
	  //cout << "[INFO] : order has to be even number" << endl;
	  return NULL;
       }else{

	  int nfrac = (order-1)/2;
	  //int nfrac = (order)/2;
	  int npow = order - nfrac;
          RooArgList* coeffs = new RooArgList();
          RooArgList* PowList = new RooArgList();
	  RooRealVar* expo[npow];
	  RooRealVar* frac[nfrac];
	  for(int i=1;i<=npow;i++){
	      expo[i] = new RooRealVar(Form("expo_%d",i),Form("expo_%d",i),TMath::Max(-9.,-1.*(i+1)),-9.,10.);
	      RooGenericPdf* PowPdf = new RooGenericPdf(Form("Pow_%d",i),Form("Pow_%d",i),"TMath::Power(@0,@1)",RooArgList(*mass,*expo[i]));
	      PowList->add(*PowPdf);
	  }
          for(int i=1;i<=nfrac;i++){
	      frac[i] = new RooRealVar(Form("frac_%d",i),Form("frac_%d",i),0.9-float(i-1)*1./nfrac,0.,1.);
	      coeffs->add(*frac[i]);
	  }

	  RooAddPdf* PowLowPdf = new RooAddPdf(Form("%s_%d",pdfType.c_str(),order),Form("%s_%d",pdfType.c_str(),order),*PowList,*coeffs,true);
	  //RooAddPdf* PowLowPdf = new RooAddPdf(Form("%s_%d",pdfType.c_str(),order),Form("%s_%d",pdfType.c_str(),order),*PowList,*coeffs);
	  if(dofit){
             PowLowPdf->fitTo(*Mggdataset,Minimizer("Minuit2","minimize"),Verbose(false),SumW2Error(true));
	     for(int i=1;i<=npow;i++)  expo[i]->setConstant(true);
             for(int i=1;i<=nfrac;i++) frac[i]->setConstant(true);
	  }

	  return PowLowPdf;
       }
    }
    if(pdfType == "Exponential"){
       if(order%2 == 0){
       //if(order%2 == 1){
	  //cout << "[INFO] : order has to be odd number" << endl;
	  //cout << "[INFO] : order has to be even number" << endl;
	  return NULL;
       }else{
	  int nfrac = (order-1)/2;
	  //int nfrac = (order)/2;
	  int npow = order - nfrac;
          RooArgList* coeffs = new RooArgList();
          RooArgList* ExpList = new RooArgList();
	  RooRealVar* expo[npow];
	  RooRealVar* frac[nfrac];
	  for(int i=1;i<=npow;i++){
	      expo[i] = new RooRealVar(Form("expo_%d",i),Form("expo_%d",i),TMath::Max(-1.,-0.04*(i+1)),-1.,0.);
	      RooExponential* exp = new RooExponential(Form("exp_%d",i),Form("exp_%d",i),*mass,*expo[i]);
	      ExpList->add(*exp);
	  }
          for(int i=1;i<=nfrac;i++){
	      frac[i] = new RooRealVar(Form("frac_%d",i),Form("frac_%d",i),0.9-float(i-1)*1./nfrac,0.,1.);
	      coeffs->add(*frac[i]);
	  }
	  RooAddPdf* ExpPdf = new RooAddPdf(Form("%s_%d",pdfType.c_str(),order),Form("%s_%d",pdfType.c_str(),order),*ExpList,*coeffs,true);
	  //RooAddPdf* ExpPdf = new RooAddPdf(Form("%s_%d",pdfType.c_str(),order),Form("%s_%d",pdfType.c_str(),order),*ExpList,*coeffs);
	  if(dofit){
             ExpPdf->fitTo(*Mggdataset,Minimizer("Minuit2","minimize"),Verbose(false),SumW2Error(true));
	     for(int i=1;i<=npow;i++)  expo[i]->setConstant(true);
             for(int i=1;i<=nfrac;i++) frac[i]->setConstant(true);
	  }
	      
	  return ExpPdf;
       }
    }
    //_______________________________________________________________________________________________________
    if(pdfType == "PowerLaw2"){

       if(order%2 != 0){
          return NULL;
       }else{
	       
          int nterm = order/2;
          string formula = "";
          int num = 1;
          RooArgList* params = new RooArgList();
	  params->add(*mass);
          RooRealVar* coeff[nterm];
          RooRealVar* expo[nterm];
          for(int i=1;i<=nterm;i++){
              coeff[i] = new RooRealVar(Form("coeff_%d",i),Form("coeff_%d",i),0.9-float(i-1)*1./nterm,0.,1.);
              expo[i] = new RooRealVar(Form("expo_%d",i),Form("expo_%d",i),TMath::Max(-9.,-1.*(i+1)),-9.,6.);
              if(i==1) formula += Form("@%d*pow(@0,@%d)",num,num+1);
              else formula += Form("+@%d*pow(@0,@%d)",num,num+1);
              params->add(*coeff[i]);
              params->add(*expo[i]);
              num+=2;
          }
          RooGenericPdf* PowLow2Pdf = new RooGenericPdf(Form("%s_%d",pdfType.c_str(),order),Form("%s_%d",pdfType.c_str(),order),formula.c_str(),*params);

          if(dofit){
             PowLow2Pdf->fitTo(*Mggdataset,Minimizer("Minuit2","minimize"),Verbose(false),SumW2Error(true));
             for(int i=1;i<=nterm;i++){  coeff[i]->setConstant(true); expo[i]->setConstant(true); }
          }
                                                                                                            
          return PowLow2Pdf;
       }
    }

    if(pdfType == "Exponential2"){
                                                                                                                                                            
       if(order%2 != 0){
          return NULL;
       }else{
               
          int nterm = order/2;
          string formula = "";
          int num = 1;
          RooArgList* params = new RooArgList();
	  params->add(*mass);
          RooRealVar* coeff[nterm];
          RooRealVar* expo[nterm];
          for(int i=1;i<=nterm;i++){
              coeff[i] = new RooRealVar(Form("coeff_%d",i),Form("coeff_%d",i),0.9-float(i-1)*1./nterm,0.,1.);
              expo[i] = new RooRealVar(Form("expo_%d",i),Form("expo_%d",i),TMath::Max(-1.,-0.04*(i+1)),-1.,0.);
              if(i==1) formula += Form("@%d*exp(@%d*@0)",num,num+1);
              else formula += Form("+@%d*exp(@%d*@0)",num,num+1);
              params->add(*coeff[i]);
              params->add(*expo[i]);
              num+=2;
          }
          RooGenericPdf* Exp2Pdf = new RooGenericPdf(Form("%s_%d",pdfType.c_str(),order),Form("%s_%d",pdfType.c_str(),order),formula.c_str(),*params);
                                                                                                                                                            
          if(dofit){
             Exp2Pdf->fitTo(*Mggdataset,Minimizer("Minuit2","minimize"),Verbose(false),SumW2Error(true));
             for(int i=1;i<=nterm;i++){  coeff[i]->setConstant(true); expo[i]->setConstant(true); }
          }
                                                                                                            
          return Exp2Pdf;
       }
    }

}

#endif
