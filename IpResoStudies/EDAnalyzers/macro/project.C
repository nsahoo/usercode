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

void project(int type=0) {
  using namespace std;
  TChain chain("residuals/tree");
  
  TString histoName;
  TString outFileName;
  
  /*
  if(type==0)  outFileName = "data.rawD0.vsPt.60bins.root";
  if(type==2)  outFileName = "data.rawDz.vsPt.60bins.root";
  if(type==0)  histoName = "data_rawD0_vsPt_n";
  if(type==2)  histoName = "data_rawDz_vsPt_n"; 
  */

  if(type==1)  outFileName = "sim.rawD0.vsPt.SETBINSbins.root";
  if(type==2)  outFileName = "sim.rawDz.vsPt.SETBINSbins.root";

  if(type==3)  outFileName = "sim.rawD0.vsEta.SETBINSbins.root";
  if(type==4)  outFileName = "sim.rawDz.vsEta.SETBINSbins.root";

  if(type==5)  outFileName = "sim.rawD0.vsPhi.SETBINSbins.root";
  if(type==6)  outFileName = "sim.rawDz.vsPhi.SETBINSbins.root";


  // --- 
  if(type==1)  histoName = "sim_rawD0_vsPt_n";
  if(type==2)  histoName = "sim_rawDz_vsPt_n";

  if(type==3)  histoName = "sim_rawD0_vsEta_n";
  if(type==4)  histoName = "sim_rawDz_vsEta_n";

  if(type==5)  histoName = "sim_rawD0_vsPhi_n";
  if(type==6)  histoName = "sim_rawDz_vsPhi_n";


  chain.Add("SET_INPUTSIM");


  cout << "type: " << type << endl;
  cout << "chain has #files: " << chain.GetListOfFiles()->GetEntries() << endl;
  //wantMore();


  //gDirectory->ls();
  //TTree* myTree = (TTree*) file->Get("residuals/tree");
  
  TFile output(outFileName,"RECREATE");

  const unsigned int nbins = SETBINS;

  TH1F* histos[nbins];
  TString hnames[nbins];

  double low,high,gap;

  if(type == 1 || type == 2){
    gap = 0.025;
    low = 0.500;
  }

  if(type == 3 || type == 4){
    gap = 0.1;
    low = -2.5;
  }

  high = low+gap;

  for(unsigned int i=0; i<nbins; ++i){
    stringstream stream;  stream << i+1;
    TString counter = stream.str();
    hnames[i] = histoName+counter;
    //cout << "hnames[" << i << "]: " << hnames[i] << endl;
    histos[i] = new TH1F(hnames[i],hnames[i],200,-2000,2000);

    stringstream lowerCut;
    stringstream higherCut;
    
    lowerCut << low;
    higherCut << high;

    TString selection;

    // selection for projection vs pt
    if(type == 1 || type == 2){
      selection = "abs(eta)<0.4 && pt>";
      selection += lowerCut.str();
      selection += " && pt<";
      selection += higherCut.str();
    }

    // selection for projection vs eta
    if(type == 3 || type == 4){
      selection = "pt>SETPTMIN && pt < SETPTMAX && eta>";
      selection += lowerCut.str();
      selection += " && eta<";
      selection += higherCut.str();
    }

    // selection for projection vs phi
    if(type == 5 || type == 6){
      //STILL TO BE DEFINED
    }

    
    cout << "selection: " << selection << endl;

    if(type==1 || type==3 || type==5)
      chain.Project(hnames[i],"d0",selection);

    if(type==2 || type==4 || type==6)
      chain.Project(hnames[i],"dz",selection);

    //histos[i]->Draw(); gPad->Update(); 
    //int pippo; cin >> pippo;
    histos[i]->Write();

    low += gap;
    high += gap;
    //wantMore();
  }

}
