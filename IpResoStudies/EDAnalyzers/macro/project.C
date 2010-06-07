#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TChain.h"
#include "TLeaf.h"

#include <sstream>
#include <iostream>

struct treeRaw{
  double pt;
  double eta;
  double phi;
  int nXLayers;
  int nMissedOut;
  int nMissedIn;
  int hasPXL;
  int    quality;
  double d0;
  double dz;
  double d0Err;
  double dzErr;
};



struct treeReso{
  double pt;
  double eta;
  double phi;
  int nXLayers;
  int nMissedOut;
  int nMissedIn;
  int hasPXL;
  int type;
  double dxyReso;
  double dzReso;
};

struct treeResp{
  double pt;
  double eta;
  double phi;
  int nXLayers;
  int nMissedOut;
  int nMissedIn;
  int hasPXL;
  int type;
  double dxyResp;
  double dzResp;
};


bool wantMore() {
                // ask if user wants more
                fprintf(stderr,"Type <CR> to continue or q to quit ==> ");
                // read first char
                int readch = getchar(), answer = readch;
                // poll out remaining chars from buffer
                while (readch != '\n' && readch != EOF) readch = getchar();
                // check first char
                return !(answer == 'q' || answer == 'Q');
}

void project(int dataset=1,
	     int variableCategory=1,
	     int variableType=1,
	     int projection=1) {
  using namespace std;

  // --- translate the input integers provided by the users into 
  //     strings to be used later in the code
  //
  TString prefix;
  switch (dataset){
  case 1:
    prefix="sim";
    break;
  case 2:
    prefix="data";
    break;
  default:
    cout << "ERROR: the dataset you specified (i.e. " << dataset 
	 << ") was not defined. Exit!" << endl;
    exit(1);
  }

  TString category;
  TString treeFolder;
  switch (variableCategory){
  case 1:
    category="raw";
    treeFolder="residuals/tree";
    break;
  case 2:
    category="reso";
    treeFolder="vertexResponsesAndTrueResolutions/tReso";
    break;
  case 3:
    category="resp";
    treeFolder="vertexResponsesAndTrueResolutions/tResp";
    break;
  default:
    cout << "ERROR: the variableCategory you specified (i.e. " << dataset 
	 << ") was not defined. Exit!" << endl;
    exit(1);
  }

  TString type;
  switch (variableType){
  case 1:
    type="d0";
    break;
  case 2:
    type="dz";
    break;
  case 3:
    type="dxyReso";
    break;
  case 4:
    type="dzReso";
    break;
  case 5:
    type="dxyResp";
    break;
  case 6:
    type="dzResp";
    break;
  case 7:
    type="d0Err";
    break;
  case 8:
    type="dzErr";
    break;
  default:
    cout << "ERROR: the variableType you specified (i.e. " << dataset 
	 << ") was not defined. Exit!" << endl;
    exit(1);
  }

  
  TString proj;
  switch (projection){
  case 1:
    proj="Pt";
    break;
  case 2:
    proj="Eta";
    break;
  case 3:
    proj="Phi";
    break;
  default:
    cout << "ERROR: the projection you specified (i.e. " << projection 
	 << ") was not defined. Exit!" << endl;
    exit(1);
  }
  // --- --- 





  // --- here the settings provided by the users are concateneted 
  //     outputFile name and histogram names are defined
  //
  TString histoName;
  TString outFileName;  
  outFileName = prefix+"."+category+"."+type+".vs"+proj+".SETBINSbins.root";
  histoName = prefix+"_"+category+"_"+type+"_vs"+proj+"_n";
  // --- ---


  TChain chain(treeFolder);
  chain.Add("SET_INPUTSIM");

  cout << "outFileName: " << outFileName << endl;
  cout << "chain has #files: " << chain.GetListOfFiles()->GetEntries() << endl;
  //wantMore();

  
  TFile output(outFileName,"RECREATE");

  const unsigned int nbins = SETBINS;
  TH1F* histos[nbins];
  TString hnames[nbins];

  double low,high,gap;
  if(projection == 1){//projection vs Pt
    gap = 0.025;
    low = 0.700;
  }
  if(projection == 2){//projection vs Eta
    gap = 0.1;
    low = -2.5;
  }

  if(projection == 3){//projection vs Phi
    gap = 2.*M_PI/nbins;
    low = -M_PI;
  }

  high = low+nbins*gap;

  cout << "xhigh,xlow,xgap: " 
       << high << " , "
       << low << " , "
       << gap << endl;

  for(unsigned int i=0; i<nbins; ++i){
    stringstream stream;  stream << i+1;
    TString counter = stream.str();
    hnames[i] = histoName+counter;
    //cout << "hnames[" << i << "]: " << hnames[i] << endl;
    histos[i] = new TH1F(hnames[i],hnames[i],SET_PBINS,SET_PLOW,SET_PHIGH);
  }

  treeRaw raw;
  treeReso reso;
  treeResp resp;
  
  TBranch* branch;
  if(category=="raw"){
    branch = chain.GetBranch("raw");
    branch->SetAddress(&raw);
  }else if(category=="reso"){
    branch = chain.GetBranch("reso");
    branch->SetAddress(&reso);
  }else if(category=="resp"){
    branch = chain.GetBranch("resp");
    branch->SetAddress(&resp);
  }else{
    cout << "ERROR: unable to set the tree branch address correctly" << endl;
    exit(1);
  }

  TH1F hBinSearch("hBinSearch","",nbins,low,high);

  //---- new faster/smarter implementation
  unsigned long int nentries = chain.GetEntries();
  
  cout << "going to loop over " << nentries << " tree entries" << endl;

  for(unsigned long int i=0; i<nentries; i++){
    //for(unsigned long int i=0; i<20000; i++){
    if(i % 10000 == 0) cout << "counter i: " << i << endl;
    chain.GetEntry(i);
    branch->GetEntry(i);

    double eta,phi,pt;
    eta = branch->GetLeaf("eta")->GetValue();
    phi = branch->GetLeaf("phi")->GetValue();
    pt  = branch->GetLeaf("pt")->GetValue();
    bool hasPXL;
    hasPXL = branch->GetLeaf("hasPXL")->GetValue();

    double var2proj = branch->GetLeaf(type)->GetValue();
    /*
    cout << "type: " << type << endl;
    cout << "projection: " << projection << endl;
    cout << "eta,phi,pt: " 
	 << eta << " , " 
	 << phi << " , " 
	 << pt  << endl;
    */

    
    // --- selection for projection vs pt
    if(projection == 1){
      bool selection = 	eta >SETETAMIN && eta < SETETAMAX  && hasPXL;
      if(!selection) continue;
      int bin = hBinSearch.FindBin(pt);
      
      if(bin>=1 && bin <=nbins) {
	histos[bin-1]->Fill(var2proj);
      }
    }
    

    // --- selection for projection vs eta
    if(projection == 2){
      bool selection = pt>SETPTMIN && pt < SETPTMAX && hasPXL;
      if(!selection) continue;
      int bin = hBinSearch.FindBin(eta);
      
      if(bin>=1 && bin <=nbins) {
	histos[bin-1]->Fill(var2proj);
      }
    }

    // --- selection for projection vs phi
    if(projection == 3){
      bool selection = pt >SETPTMIN && pt < SETPTMAX 
	&& eta >SETETAMIN && eta < SETETAMAX
	&& hasPXL;
      if(!selection) continue;
      int bin = hBinSearch.FindBin(phi);
      
      if(bin>=1 && bin <=nbins) {
	histos[bin-1]->Fill(var2proj);
      }
    }

  }
  
  for(int i=0; i<nbins; i++){
    histos[i]->Write();
  }
}


