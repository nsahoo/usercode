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
  double p;
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
  double p;
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
  double p;
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
  case 4:
    proj="P";
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
  if(projection==1 || projection==4)
    outFileName = prefix+"."+category+"."+type+".vs"+proj+".root";
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


  if(projection == 1 || projection ==4){
    //const int nbinsPT=100; //solution1
    //const int nbinsPT=113;   //up to 10
    //const int nbinsPT=118;   //solution3 (up to 20)
    //const int nbinsPT=116;   //solution4 (up to 20)
    //const int nbinsPT=116;   //solution5 (up to 20)
    //const int nbinsPT=115;   //solution6 (up to 20)
    const int nbinsPT=114;   //solution7 (up to 20)
    double xbins[nbinsPT+1];
    double lowestEdge = 0.7;
    double firstSectionGap = 0.025;
    for(int i=1; i<=92; i++){
      xbins[i-1]=lowestEdge+firstSectionGap*(i-1);
    }//first 92 bins set

    /*
    //solution1
    xbins[93-1] = 3.0;
    xbins[94-1] = 4.0;
    xbins[95-1] = 5.0;
    xbins[96-1] = 6.0;
    xbins[97-1] = 7.0;
    xbins[98-1] = 8.0;
    xbins[99-1] = 9.0;
    xbins[100-1]=10.0;
    xbins[101-1]=99.0;    
    */

    
    xbins[93-1] = 3.0;
    xbins[94-1] = 3.1; 
    xbins[95-1] = 3.2;    
    xbins[96-1] = 3.3;
    xbins[97-1] = 3.4;
    xbins[98-1] = 3.5; 
    xbins[99-1] = 3.6;    
    xbins[100-1] = 3.7;
    xbins[101-1] = 3.8;
    xbins[102-1] = 3.9; //next 10 bins set

    xbins[103-1] = 4.0; 
    xbins[104-1] = 4.2;    
    xbins[105-1] = 4.4;
    xbins[106-1] = 4.6;
    xbins[107-1] = 4.8; //next 5 bins set 

    xbins[108-1] = 5.0;
    xbins[109-1] = 5.5; //next 2 bins set

    /*
    //solution3
    xbins[110-1] = 6.0;
    xbins[111-1] = 7.0;
    xbins[112-1] = 8.0;
    xbins[113-1] = 9.0;
    xbins[114-1] = 10.0; //last line for plot up to 10 GeV
    
    
    xbins[115-1] = 11.0;
    xbins[116-1] = 12.0;

    xbins[117-1] = 13.0;
    xbins[118-1] = 15.0;

    xbins[119-1] = 20.0;//last bin upper edge
    */

 
    /*
    //solution4
    xbins[110-1] = 6.0;
    xbins[111-1] = 7.0;
    xbins[112-1] = 8.0;
    xbins[113-1] = 9.0;
    xbins[114-1] = 11.0; 
    xbins[115-1] = 13.0;
    xbins[116-1] = 15.0;
    xbins[117-1] = 20.0;//last bin upper edge
    */

    /*
    //solution5
    xbins[110-1] = 6.0;
    xbins[111-1] = 7.0;
    xbins[112-1] = 8.0;
    xbins[113-1] = 9.0;
    xbins[114-1] = 10.0;
    xbins[115-1] = 12.0; 
    xbins[116-1] = 14.0;
    xbins[117-1] = 20.0;
    */

    /*
    //solution6
    xbins[110-1] = 6.0;
    xbins[111-1] = 7.0;
    xbins[112-1] = 8.0;
    xbins[113-1] = 9.0;
    xbins[114-1] = 11.0;
    xbins[115-1] = 14.0;
    xbins[116-1] = 20.0;
    */

    //solution7
    xbins[110-1] = 6.0;
    xbins[111-1] = 7.0;
    xbins[112-1] = 8.0;
    xbins[113-1] = 9.0;
    xbins[114-1] = 14.0;
    xbins[115-1] = 20.0;

    for(int i=0; i<=nbinsPT; i++){
      cout << "xbin, content: " << i+1 << " , " <<xbins[i] << endl;
    }

    hBinSearch.SetBins(nbinsPT,xbins);
  }

  //---- new faster/smarter implementation
  unsigned long int nentries = chain.GetEntries();
  
  cout << "going to loop over " << nentries << " tree entries" << endl;

  for(unsigned long int i=0; i<nentries; i++){
    //for(unsigned long int i=0; i<20000; i++){
    if(i % 10000 == 0) cout << "counter i: " << i << endl;
    chain.GetEntry(i);
    branch->GetEntry(i);

    double eta,phi,pt,p;
    eta = branch->GetLeaf("eta")->GetValue();
    phi = branch->GetLeaf("phi")->GetValue();
    pt  = branch->GetLeaf("pt")->GetValue();
    p   = branch->GetLeaf("p")->GetValue();
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

    if(projection == 4){
      bool selection = 	eta >SETETAMIN && eta < SETETAMAX  && hasPXL;
      if(!selection) continue;
      int bin = hBinSearch.FindBin(p);
      
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


