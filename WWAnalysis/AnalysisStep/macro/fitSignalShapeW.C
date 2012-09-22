#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "./scales2.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "Math/MinimizerOptions.h"

using namespace RooFit;



#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
}

float weightTrue2011(float input);
float weightTrue2012(float input);


void fitSignalShapeW(int massBin=120, int id=1120, int ch=0, int year = 2011,
		     float lumi=10, bool doSfLepton=false,double rangeLow=105., double rangeHigh=130.,
		     double bwSigma=1.,
		     double fitValues[5]=0, double fitErrors[5]=0);

void all(int channels=0,int year=2011, bool doSfLepton=true){
  /*
    channels = 0 --> 4mu
    channels = 1 --> 4el
    channels = 2 --> 2mu2el

   */
  //init();

  double bwSigma[30];
  int mass[30]; int id[30]; double xLow[30]; double xHigh[30];  
  int maxMassBin;

  if(year==2011){
    init(true);
    mass[0] = 115; id[0]=200; xLow[0] = 100; xHigh[0] = 125; bwSigma[0] = 3.1/1000.;
    mass[1] = 120; id[1]=201; xLow[1] = 105; xHigh[1] = 130; bwSigma[1] = 3.5/1000.;
    mass[2] = 130; id[2]=202; xLow[2] = 115; xHigh[2] = 140; bwSigma[2] = 4.9/1000.;
    mass[3] = 140; id[3]=203; xLow[3] = 125; xHigh[3] = 150; bwSigma[3] = 8.1/1000.;
    mass[4] = 150; id[4]=204; xLow[4] = 135; xHigh[4] = 160; bwSigma[4] = 1.7/100.;
    mass[5] = 160; id[5]=205; xLow[5] = 145; xHigh[5] = 170; bwSigma[5] = 8.3/100.;
    mass[6] = 170; id[6]=206; xLow[6] = 155; xHigh[6] = 180; bwSigma[6] = 3.8/10.;
    mass[7] = 180; id[7]=207; xLow[7] = 165; xHigh[7] = 190; bwSigma[7] = 6.3/10.;
    mass[8] = 190; id[8]=208; xLow[8] = 170; xHigh[8] = 200; bwSigma[8] = 1.04;
    mass[9] = 200; id[9]=209; xLow[9] = 175; xHigh[9] = 210; bwSigma[9] = 1.43;
    mass[10] = 210; id[10]=210; xLow[10] = 185; xHigh[10] = 220; bwSigma[10] = 1.85;
    mass[11] = 220; id[11]=211; xLow[11] = 190; xHigh[11] = 230; bwSigma[11] = 2.31;
    maxMassBin = 12;
  }


  if(year==2012){
    init(false);
    mass[0] = 115; id[0]=1115; xLow[0] = 100; xHigh[0] = 125; bwSigma[0] = 3.1/1000.;
    mass[1] = 117; id[1]=1117; xLow[1] = 102; xHigh[1] = 127; bwSigma[1] = 3.2/1000.;
    mass[2] = 119; id[2]=1119; xLow[2] = 104; xHigh[2] = 129; bwSigma[2] = 3.4/1000.;
    mass[3] = 120; id[3]=1120; xLow[3] = 105; xHigh[3] = 130; bwSigma[3] = 3.5/1000.;
    mass[4] = 121; id[4]=1121; xLow[4] = 106; xHigh[4] = 131; bwSigma[4] = 3.6/1000.;
    mass[5] = 123; id[5]=1123; xLow[5] = 108; xHigh[5] = 133; bwSigma[5] = 3.8/1000.;
    mass[6] = 124; id[6]=1124; xLow[6] = 109; xHigh[6] = 134; bwSigma[6] = 3.9/1000.;
    mass[7] = 125; id[7]=1125; xLow[7] = 110; xHigh[7] = 135; bwSigma[7] = 4.0/1000.;
    mass[8] = 126; id[8]=1126; xLow[8] = 111; xHigh[8] = 136; bwSigma[8] = 4.2/1000.;
    mass[9] = 127; id[9]=1127; xLow[9] = 112; xHigh[9] = 137; bwSigma[9] = 4.3/1000.;
    mass[10] = 145; id[10]=1145; xLow[10] = 130; xHigh[10] = 155; bwSigma[10] = 1.1/100.;
    //
    mass[11] = 150; id[11]=1150; xLow[11] = 135; xHigh[11] = 160; bwSigma[11] = 1.7/100.;
    mass[12] = 180; id[12]=1180; xLow[12] = 150; xHigh[12] = 195; bwSigma[12] = 6.3/10.;
    mass[13] = 200; id[13]=1200; xLow[13] = 165; xHigh[13] = 220; bwSigma[13] = 1.4;
    mass[14] = 220; id[14]=1220; xLow[14] = 175; xHigh[14] = 250; bwSigma[14] = 2.3;
    maxMassBin = 15;
  }
  // -----------------------



  double massV[30],massE[30];
  for(int i=0; i<maxMassBin;++i){
    massV[i]=mass[i];
    massE[i]=0;
  }

  double aVal[30],aErr[30];
  double nVal[30],nErr[30];
  double meanCBVal[30],meanCBErr[30];
  double sigmaCBVal[30],sigmaCBErr[30];
  double meanBWVal[30],meanBWErr[30];

  double fitValues[5];
  double fitErrors[5];

  double extendL(1),extendH(1);

  if(channels==0) {extendL=1.0;extendH=1.0;}
  if(channels==1) {extendL=0.90;extendH=1.05;}
  if(channels==2) {extendL=0.95;extendH=1.04;}



  //fitSignalShapeW(mass[0],channels,10.,doSfLepton,max(95.,xLow[0]*0.90),xHigh[0]*1.1,bwSigma[12],
  // 		  fitValues,fitErrors);  Wait();
  
  for(int i=0; i<maxMassBin;++i){

    fitSignalShapeW(mass[i],id[i],channels,year,10.,doSfLepton,max(95.,xLow[i]*extendL),xHigh[i]*extendH,bwSigma[i],
		    fitValues,fitErrors);  
  
    cout << "a value,error: " << fitValues[0] << " , " << fitErrors[0] << endl; 
    aVal[i]=fitValues[0]; aErr[i]=fitErrors[0];

    cout << "n value,error: " << fitValues[3] << " , " << fitErrors[3] << endl; 
    nVal[i]=fitValues[3]; nErr[i]=fitErrors[3];

    cout << "meanCB value,error: " << fitValues[1] << " , " << fitErrors[1] << endl;
    meanCBVal[i]=fitValues[1]; meanCBErr[i]=fitErrors[1];
    
    cout << "sigmaCB value,error: " << fitValues[4] << " , " << fitErrors[4] << endl; 
    sigmaCBVal[i]=fitValues[4]; sigmaCBErr[i]=fitErrors[4];

    cout << "meanBW value,error: " << fitValues[2] << " , " << fitErrors[2] << endl; 
    meanBWVal[i]=fitValues[2]; meanBWErr[i]=fitErrors[2];

    //Wait();
  }
  

  TGraphErrors* gA = new TGraphErrors(maxMassBin,massV,aVal,massE,aErr);
  TGraphErrors* gN = new TGraphErrors(maxMassBin,massV,nVal,massE,nErr);
  TGraphErrors* gMeanCB = new TGraphErrors(maxMassBin,massV,meanCBVal,massE,meanCBErr);
  TGraphErrors* gSigmaCB = new TGraphErrors(maxMassBin,massV,sigmaCBVal,massE,sigmaCBErr);
  TGraphErrors* gMeanBW = new TGraphErrors(maxMassBin,massV,meanBWVal,massE,meanBWErr);

  gA->SetMarkerStyle(20);   gA->SetMarkerSize(1);
  gN->SetMarkerStyle(20);   gN->SetMarkerSize(1);
  gMeanCB->SetMarkerStyle(20);   gMeanCB->SetMarkerSize(1);
  gSigmaCB->SetMarkerStyle(20);   gSigmaCB->SetMarkerSize(1);
  gMeanBW->SetMarkerStyle(20);   gMeanBW->SetMarkerSize(1);
  

  gA->SetTitle("");
  gA->GetXaxis()->SetTitle("mass (GeV)");
  gA->GetYaxis()->SetTitle("CB a-parameter");

  gN->SetTitle("");
  gN->GetXaxis()->SetTitle("mass (GeV)");
  gN->GetYaxis()->SetTitle("CB n-parameter");

  gMeanCB->SetTitle("");
  gMeanCB->GetXaxis()->SetTitle("mass (GeV)");
  gMeanCB->GetYaxis()->SetTitle("CB mean");

  gSigmaCB->SetTitle("");
  gSigmaCB->GetXaxis()->SetTitle("mass (GeV)");
  gSigmaCB->GetYaxis()->SetTitle("CB sigma");

  gMeanBW->SetTitle("");
  gMeanBW->GetXaxis()->SetTitle("mass (GeV)");
  gMeanBW->GetYaxis()->SetTitle("BW mean");



  gA->Fit("pol0"); gA->Draw("Ap"); gPad->Update(); gPad->Print("aFit.pdf"); Wait();
  gN->Fit("pol1"); gN->Draw("Ap"); gPad->Update(); gPad->Print("nFit.pdf"); Wait();
  gMeanCB->Fit("pol1"); gMeanCB->Draw("Ap"); gPad->Update(); gPad->Print("meanFit.pdf"); Wait();
  gSigmaCB->Fit("pol1"); gSigmaCB->Draw("Ap"); gPad->Update(); gPad->Print("sigmaFit.pdf"); Wait();
  //gMeanBW->Fit("pol1"); gMeanBW->Draw("Ap"); gPad->Update(); Wait();


}



void fitSignalShapeW(int massBin,int id, int channels, int year, 
		     float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
		     double bwSigma,
		     double fitValues[5], double fitErrors[5]){
 // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  gStyle->SetOptStat("iourme");
  //gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(11);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 

  ROOT::Math::MinimizerOptions::SetDefaultTolerance( 1.E-7);

  TH1F* hPuWeights;


  stringstream ggFileName,vbfFileName;
  ggFileName << "/data/hzz4l/step2/ichep2012/" << year << "/MC/hzzTree_id" << id << ".root"; 
  //ggFileName << "/data/hzz4l/step2/ichep2012/2011/MC/hzzTree_id207.root";


  TFile* ggFile = new TFile(ggFileName.str().c_str()); 

  TTree* ggTree = (TTree*) ggFile->Get("zz4lTree/probe_tree");


  TH1F* ggMass = new TH1F("ggMass","ggMass",200,100,300.);


  float nvtx,channel;
  //float ntrueint;
  float nTrueInt,nObsInt,rho,mass,m4l,trueInt;

  float pt1,eta1,id1,pt2,eta2,id2,pt3,eta3,id3,pt4,eta4,id4;
  float sfLepton(1.);


int  nentries = ggTree->GetEntries();
 
 
  //--- ggTree part
  ggTree->SetBranchAddress("channel",&channel);
  ggTree->SetBranchAddress("numTrueInteractions",&nTrueInt);
  ggTree->SetBranchAddress("rho",&rho);
  ggTree->SetBranchAddress("mass",&mass);
  ggTree->SetBranchAddress("m4l",&m4l);
  ggTree->SetBranchAddress("l1pt",&pt1);ggTree->SetBranchAddress("l1eta",&eta1);ggTree->SetBranchAddress("l1pdgId",&id1);
  ggTree->SetBranchAddress("l2pt",&pt2);ggTree->SetBranchAddress("l2eta",&eta2);ggTree->SetBranchAddress("l2pdgId",&id2);
  ggTree->SetBranchAddress("l3pt",&pt3);ggTree->SetBranchAddress("l3eta",&eta3);ggTree->SetBranchAddress("l3pdgId",&id3);
  ggTree->SetBranchAddress("l4pt",&pt4);ggTree->SetBranchAddress("l4eta",&eta4);ggTree->SetBranchAddress("l4pdgId",&id4);

   /*
  for(int k=0; k<nentries; k++){
    ggTree->GetEvent(k);
    if(doSfLepton) sfLepton = getSF(pt1, eta1, id1,typeErrSfLep)
		     *getSF(pt2, eta2, id2,typeErrSfLep)
		     *getSF(pt3, eta3, id3,typeErrSfLep)
		     *getSF(pt4, eta4, id4,typeErrSfLep);
    else sfLepton = 1.0;
    ggMass->Fill(mass,weightTrue2012(nTrueInt)*xsecweights[1000+massBin]*lumi*sfLepton);
  }
  */





  /*TCanvas* myC = new TCanvas("myCanvas","myCanvas",0,0,500,500); 
   myC->cd(1);

   vbfMass->Draw(); gPad->Update(); Wait();
  */

  //cout << "yield gg: " << ggMass->Integral() << endl;
  //cout << "yield vbf: " << vbfMass->Integral() << endl;
  //cout << "yield total for mass " << massBin << " is: " << ggMass->Integral()<< endl;
  //ggMass->Draw(); gPad->Update(); Wait();




  //--- rooFit part
  double xMin,xMax,xInit;
  xInit = (double) massBin;
  xMin = rangeLow;
  xMax = rangeHigh ;
  //xMin = 114.;
  //xMax = 122 ;



  /*
  TCut cut = "channel == 0";
  RooRealVar x("mass","mass",xInit,xMin,xMax);
  RooRealVar ch("channel","channel",0,0,10);
  RooRealVar w("baseW","baseW",1.0,0.,1000.);
  RooArgSet ntupleVarSet(x,ch,w);
  RooDataSet dataset("mass4l","mass4l",ntupleVarSet,WeightVar(w),Import(*ggTree),Cut("channel==0"));
  */

  //---------  
  TCut cut = "channel == 0";
  RooRealVar x("mass","mass (GeV)",xInit,xMin,xMax);
  //RooRealVar x("m4l","m4l",xInit,xMin,xMax);
  RooRealVar w("myW","myW",1.0,0.,1000.);
  RooArgSet ntupleVarSet(x,w);
  RooDataSet dataset("mass4l","mass4l",ntupleVarSet,WeightVar("myW"));


  for(int k=0; k<nentries; k++){
    ggTree->GetEvent(k);
    if(doSfLepton) sfLepton = getSF(pt1, eta1, id1)
		     *getSF(pt2, eta2, id2)
		     *getSF(pt3, eta3, id3)
		     *getSF(pt4, eta4, id4);
    else sfLepton = 1.0;

    //if(channel==ch1 || channel==ch2){      
    if(channels==0){if(channel != 0) continue;}
    if(channels==1){if(channel != 1) continue;}
    if(channels==2){if(channel != 2 && channel != 3) continue;}

    ntupleVarSet.setRealValue("mass",mass);
    ntupleVarSet.setRealValue("m4l",m4l);
    //double localW = weightTrue2012(nTrueInt)*xsecweights[1000+massBin]*lumi*sfLepton;
    //double localW = weightTrue2012(nTrueInt)*xsecweights[1000+massBin]*lumi;
    //double localW = xsecweights[1000+massBin];
    //double localW = weightTrue2012(nTrueInt);
    //double localW = 1.;
    double localW(1);
    if(year==2012) localW = weightTrue2012(nTrueInt);
    else if(year==2011) localW = weightTrue2011(nTrueInt);    
    localW = localW*sfLepton;

    ntupleVarSet.setRealValue("myW",localW);
    if(x.getVal()>xMin && x.getVal()<xMax)
      dataset.add(ntupleVarSet, localW);
    
  }
  //---------

  //cout << "dataset n entries: " << dataset.sumEntries() << endl;
  //cout << "reduced data n entries: " << reducedData->sumEntries() << endl;


  //--- simple CrystalBall
  RooRealVar mean("mean","mean of gaussian",0,-5.0,5.0) ;
  RooRealVar sigma("sigma","width of gaussian",1.5,0.,5.); 
  RooRealVar a("a","a",1.46,0.,4.);
  RooRealVar n("n","n",1.92,0.,25.);   
  RooCBShape CBall("CBall","Crystal ball",x, mean,sigma, a,n);
  

  //--- simple Gaussian
  RooRealVar mean2("mean2","mean2 of gaussian",xInit,xInit*0.80,xInit*1.2) ;
  RooRealVar sigma2("sigma2","width2 of gaussian",10.,5.,500.); 
  RooGaussian tailCatcher("tailCatcher","tailCatcher",x,mean2,sigma2);
  RooRealVar fsig("fsig","signal fraction",0.95,0.7,1.);


  //--- Breit-Wigner
  RooRealVar mean3("mean3","mean3",xInit) ;
  RooRealVar sigma3("sigma3","width3",bwSigma); 
  RooRealVar scale3("scale3","scale3 ",1.); 
  RooRelBWUFParam bw("bw","bw",x,mean3,scale3);

  //RooAddPdf model("model","model",RooArgList(CBall,tailCatcher),fsig);
  //RooCBShape model("model","model",x, mean,sigma, a,n);
  x.setBins(10000,"fft");
  RooFFTConvPdf model("model","model",x,bw,CBall);

  

  model.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8));

  stringstream frameTitle;
  if(channels==0){frameTitle << "4mu, mH = ";}
  if(channels==1){frameTitle << "4el, mH = ";}
  if(channels==2){frameTitle << "2mu2el, mH = ";}
  frameTitle << massBin << " GeV";

  RooPlot* xframe = x.frame(Title(frameTitle.str().c_str() )) ;
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2) );
  model.plotOn(xframe);
  //model.paramOn(xframe);


  stringstream nameFile;
  nameFile << "fitM" << massBin << ".pdf";
  xframe->Draw(); gPad->Update(); gPad->Print(nameFile.str().c_str());


  if(fitValues!=0){
    fitValues[0] = a.getVal();
    fitValues[1] = mean.getVal();
    fitValues[2] = mean3.getVal();
    fitValues[3] = n.getVal();
    fitValues[4] = sigma.getVal();
  }

  if(fitErrors!=0){
    fitErrors[0] = a.getError();
    fitErrors[1] = mean.getError();
    fitErrors[2] = mean3.getError();
    fitErrors[3] = n.getError();
    fitErrors[4] = sigma.getError();
  }

}



float weightTrue2012(float input){
  float w[240];


w[0]= 1;
w[1]= 1;
w[2]= 0.222451;
w[3]= 0.0658851;
w[4]= 1;
w[5]= 1;
w[6]= 0.150902;
w[7]= 0.202205;
w[8]= 1;
w[9]= 1;
w[10]= 1.50116;
w[11]= 2.79375;
w[12]= 0.198341;
w[13]= 0.246893;
w[14]= 0.28116;
w[15]= 0.449377;
w[16]= 0.553276;
w[17]= 1.48919;
w[18]= 2.15249;
w[19]= 3.62415;
w[20]= 4.33041;
w[21]= 3.57192;
w[22]= 4.99603;
w[23]= 7.79303;
w[24]= 8.04276;
w[25]= 8.05557;
w[26]= 12.9364;
w[27]= 9.9036;
w[28]= 14.6975;
w[29]= 13.3387;
w[30]= 10.9734;
w[31]= 12.6077;
w[32]= 11.5617;
w[33]= 10.8107;
w[34]= 14.5043;
w[35]= 17.8497;
w[36]= 11.8817;
w[37]= 9.6805;
w[38]= 12.2255;
w[39]= 10.1117;
w[40]= 10.2482;
w[41]= 11.5398;
w[42]= 9.35737;
w[43]= 9.90259;
w[44]= 9.19216;
w[45]= 7.57377;
w[46]= 7.94847;
w[47]= 7.15578;
w[48]= 5.63016;
w[49]= 5.35972;
w[50]= 5.05791;
w[51]= 3.35313;
w[52]= 3.60582;
w[53]= 3.35256;
w[54]= 2.49496;
w[55]= 2.28219;
w[56]= 2.21227;
w[57]= 1.76362;
w[58]= 1.68533;
w[59]= 1.62149;
w[60]= 1.34263;
w[61]= 1.30646;
w[62]= 1.21918;
w[63]= 1.10347;
w[64]= 1.08544;
w[65]= 1.0251;
w[66]= 0.907123;
w[67]= 0.905997;
w[68]= 0.869217;
w[69]= 0.816708;
w[70]= 0.76043;
w[71]= 0.714367;
w[72]= 0.679723;
w[73]= 0.665294;
w[74]= 0.609956;
w[75]= 0.586386;
w[76]= 0.548999;
w[77]= 0.521088;
w[78]= 0.4929;
w[79]= 0.453545;
w[80]= 0.44546;
w[81]= 0.406266;
w[82]= 0.378486;
w[83]= 0.347898;
w[84]= 0.337097;
w[85]= 0.313674;
w[86]= 0.291392;
w[87]= 0.283346;
w[88]= 0.25272;
w[89]= 0.244178;
w[90]= 0.228673;
w[91]= 0.211327;
w[92]= 0.19084;
w[93]= 0.179408;
w[94]= 0.169234;
w[95]= 0.157131;
w[96]= 0.143818;
w[97]= 0.140968;
w[98]= 0.124021;
w[99]= 0.118273;
w[100]= 0.109751;
w[101]= 0.0977754;
w[102]= 0.0967206;
w[103]= 0.0870401;
w[104]= 0.0826372;
w[105]= 0.0746777;
w[106]= 0.0698592;
w[107]= 0.0656062;
w[108]= 0.0601853;
w[109]= 0.057892;
w[110]= 0.0517871;
w[111]= 0.0512109;
w[112]= 0.0465423;
w[113]= 0.0403982;
w[114]= 0.0443631;
w[115]= 0.0399185;
w[116]= 0.0338933;
w[117]= 0.0354274;
w[118]= 0.0310775;
w[119]= 0.0261122;
w[120]= 0.0280878;
w[121]= 0.0264014;
w[122]= 0.021299;
w[123]= 0.0245197;
w[124]= 0.0221076;
w[125]= 0.0189236;
w[126]= 0.0202148;
w[127]= 0.0177248;
w[128]= 0.0163634;
w[129]= 0.0188307;
w[130]= 0.0144512;
w[131]= 0.0134599;
w[132]= 0.0143315;
w[133]= 0.0130668;
w[134]= 0.0108666;
w[135]= 0.0162516;
w[136]= 0.0126035;
w[137]= 0.0102154;
w[138]= 0.0154442;
w[139]= 0.00959973;
w[140]= 0.0106827;
w[141]= 0.0146624;
w[142]= 0.0155156;
w[143]= 0.00761674;
w[144]= 0.0187999;
w[145]= 0.0135013;
w[146]= 0.0160794;
w[147]= 0.0180586;
w[148]= 0.00905508;
w[149]= 0.00869858;
w[150]= 0.0193968;
w[151]= 0.0209201;
w[152]= 0.0084405;
w[153]= 0.0407657;
w[154]= 0.0109116;
w[155]= 0.0262218;
w[156]= 0.0104767;
w[157]= 1;
w[158]= 1;
w[159]= 1;
w[160]= 1;
w[161]= 0.00658031;
w[162]= 0.0051814;
w[163]= 1;
w[164]= 1;
w[165]= 1;
w[166]= 1;
w[167]= 1;
w[168]= 1;
w[169]= 1;
w[170]= 1;
w[171]= 1;
w[172]= 1;
w[173]= 1;
w[174]= 1;
w[175]= 1;
w[176]= 1;
w[177]= 1;
w[178]= 1;
w[179]= 1;
w[180]= 1;
w[181]= 1;
w[182]= 1;
w[183]= 1;
w[184]= 1;
w[185]= 1;
w[186]= 1;
w[187]= 1;
w[188]= 1;
w[189]= 1;
w[190]= 1;
w[191]= 1;
w[192]= 1;
w[193]= 1;
w[194]= 1;
w[195]= 1;
w[196]= 1;
w[197]= 1;
w[198]= 1;
w[199]= 1;
w[200]= 1;
w[201]= 1;
w[202]= 1;
w[203]= 1;
w[204]= 1;
w[205]= 1;
w[206]= 1;
w[207]= 1;
w[208]= 1;
w[209]= 1;
w[210]= 1;
w[211]= 1;
w[212]= 1;
w[213]= 1;
w[214]= 1;
w[215]= 1;
w[216]= 1;
w[217]= 1;
w[218]= 1;
w[219]= 1;
w[220]= 1;
w[221]= 1;
w[222]= 1;
w[223]= 1;
w[224]= 1;
w[225]= 1;
w[226]= 1;
w[227]= 1;
w[228]= 1;
w[229]= 1;
w[230]= 1;
w[231]= 1;
w[232]= 1;
w[233]= 1;
w[234]= 1;
w[235]= 1;
w[236]= 1;
w[237]= 1;
w[238]= 1;
 w[239]= 1;



 TH1F h("boh","boh",240,0.,80.);
 
 for(int k=0;k<240;k++){
   h.SetBinContent(k+1,w[k]);
 }
 
 //h.Draw(); gPad->Update(); Wait();
 
 return h.GetBinContent(h.FindBin(input));

}






float weightTrue2011(float input){
  if(input>50) 
    return 1;

    
  float w[50];


w[0]= 0.212929;
w[1]= 0.0208114;
w[2]= 0.0584048;
w[3]= 0.538898;
w[4]= 1.357;
w[5]= 1.49913;
w[6]= 1.42247;
w[7]= 1.35904;
w[8]= 1.29946;
w[9]= 1.27925;
w[10]= 1.37845;
w[11]= 1.71246;
w[12]= 1.5291;
w[13]= 1.35234;
w[14]= 1.22215;
w[15]= 1.0155;
w[16]= 1.01137;
w[17]= 0.395465;
w[18]= 0.230984;
w[19]= 0.109883;
w[20]= 0.0433739;
w[21]= 0.0111497;
w[22]= 0.00408801;
w[23]= 0.00115678;
w[24]= 0.000365505;
w[25]= 0.000112391;
w[26]= 3.83894e-05;
w[27]= 1.60651e-05;
w[28]= 4.81412e-06;
w[29]= 1.39717e-06;
w[30]= 1.92368e-06;
w[31]= 4.10748e-06;
w[32]= 2.33157e-05;
w[33]= 4.0181e-05;
w[34]= 4.87786e-05;
w[35]= 0.00194128;
w[36]= 8.97414e-05;
w[37]= 1;
w[38]= 1;
w[39]= 0.000162709;
w[40]= 1;
w[41]= 1;
w[42]= 1;
w[43]= 1;
w[44]= 1;
w[45]= 1;
w[46]= 1;
w[47]= 1;
w[48]= 1;
w[49]= 1;


 TH1F h("boh","boh",50,0.,50.);
 
 for(int k=0;k<50;k++){
   h.SetBinContent(k+1,w[k]);
 }
 
 return h.GetBinContent(h.FindBin(input));

}
