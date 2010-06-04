#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TH2F.h"
#include "TTree.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TChain.h"

#include <sstream>
#include <iostream>

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
  high = low+gap;

  for(unsigned int i=0; i<nbins; ++i){
    stringstream stream;  stream << i+1;
    TString counter = stream.str();
    hnames[i] = histoName+counter;
    //cout << "hnames[" << i << "]: " << hnames[i] << endl;
    histos[i] = new TH1F(hnames[i],hnames[i],SET_PBINS,SET_PLOW,SET_PHIGH);

    stringstream lowerCut;
    stringstream higherCut;
    
    lowerCut << low;
    higherCut << high;

    TString selection;

    // --- selection for projection vs pt
    if(projection == 1){
      selection = "abs(eta)<0.4 && pt>";
      selection += lowerCut.str();
      selection += " && pt<";
      selection += higherCut.str();
      selection += " && hasPXL"; //additional selection to improve purity of prompt tracks
    }

    // --- selection for projection vs eta
    if(projection == 2){
      selection = "pt>SETPTMIN && pt < SETPTMAX && eta>";
      selection += lowerCut.str();
      selection += " && eta<";
      selection += higherCut.str();
      selection += " && hasPXL"; //additional selection to improve purity of prompt tracks
    }

    // --- selection for projection vs phi
    if(projection == 3){
      //STILL TO BE DEFINED
    }   

    TString variableToProject =category+"."+type;
    cout << "selection: " << selection  
	 << " . VariableToProject: " << variableToProject << endl;

    chain.Project(hnames[i],variableToProject,selection);


    //histos[i]->Draw(); gPad->Update(); 
    histos[i]->Write();

    low += gap;
    high += gap;
    //wantMore();
  }

}
