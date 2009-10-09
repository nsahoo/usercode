#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#if !defined(__CINT__) && !defined(__MAKECINT__)
//Headers for the data items
#include <DataFormats/MuonReco/interface/Muon.h>
//#include <BrezzaQuarkonia/DataFormats/interface/QuarkoniaHypothesis.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

//
// --------- temporary workaround -----------
//#include "OctoberXTracking/DataFormats/interface/SilentMLscribe.h"

#endif

#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "Math/VectorUtil.h"
//#include "SetAlias.C"




void setInput(std::vector<std::string>& fileNames,int sample);

//void massMacro(int option=0) {
void spectrum(string outHistoName, int sample) {
  
  /*
  // --------- temporary workaround -----------
  SilentMLscribe scribe;
  edm::MessageLoggerQ::setMLscribe_ptr(&scribe);
  edm::MessageDrop::instance()->messageLoggerScribeIsRunning = edm::MLSCRIBE_RUNNING_INDICATOR;
  // --------- temporary workaround -----------
  */  


  /*
    option = 0 --> to be defined
   */


  // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  //gStyle->SetOptStat("iourme");
  gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(1111);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 

  //using namespace reco;
  using namespace std;
  using namespace ROOT::Math::VectorUtil;


  TCanvas* myC = new TCanvas("myCanvas","myCanvas",600,600); myC->cd(1);

  const int etaBins(25);
  const int ptBins(8);
  const int phiBins(20);

  //--- for JPis ---
  const double simMass(3.096916);
  const double massRangeMin(0.0);  const double massRangeMax(20.0);
  const double massRangeMinb(2.5);  const double massRangeMaxb(3.6);
  const double massRangeMin2(0.0);  const double massRangeMax2(20.0);
  const double biasRangeMin(3.08);  const double biasRangeMax(3.12);
  const double sigmaRangeMin(0.0);  const double sigmaRangeMax(0.075);


  TH1F* recoMass0 = new TH1F("recoJmass0","reco jpsi mass 0",400, massRangeMin,massRangeMax);
  TH1F* recoMass1 = new TH1F("recoJmass1","reco jpsi mass 1",400, massRangeMin,massRangeMax);
  TH1F* recoMass2 = new TH1F("recoJmass2","reco jpsi mass 2",400, massRangeMin,massRangeMax);
  TH1F* recoMass3 = new TH1F("recoJmass3","reco jpsi mass 3",400, massRangeMin,massRangeMax);
  TH1F* recoMass4 = new TH1F("recoJmass4","reco jpsi mass 4",400, massRangeMin,massRangeMax);
  TH1F* recoMass5 = new TH1F("recoJmass5","reco jpsi mass 5",400, massRangeMin,massRangeMax);
  TH1F* recoMass6 = new TH1F("recoJmass6","reco jpsi mass 6",400, massRangeMin,massRangeMax);

  TH1F* recoMass0b = new TH1F("recoJmass0b","reco jpsi mass 0b",100, massRangeMinb,massRangeMaxb);
  TH1F* recoMass1b = new TH1F("recoJmass1b","reco jpsi mass 1b",100, massRangeMinb,massRangeMaxb);
  TH1F* recoMass2b = new TH1F("recoJmass2b","reco jpsi mass 2b",100, massRangeMinb,massRangeMaxb);
  TH1F* recoMass3b = new TH1F("recoJmass3b","reco jpsi mass 3b",100, massRangeMinb,massRangeMaxb);
  TH1F* recoMass4b = new TH1F("recoJmass4b","reco jpsi mass 4b",100, massRangeMinb,massRangeMaxb);
  TH1F* recoMass5b = new TH1F("recoJmass5b","reco jpsi mass 5b",100, massRangeMinb,massRangeMaxb);
  TH1F* recoMass6b = new TH1F("recoJmass6b","reco jpsi mass 6b",100, massRangeMinb,massRangeMaxb);


  TH1F* vtxNChi2 = new TH1F("vtxNChi2","vtx Norm Chi2",50, 0.,10.);
  TF1 *myGaus = new TF1("myGaus","gaus(0)"); 

  /*
  TFile* file = new 
    TFile("/home/users/mangano/Quarkonia/samples/myOutput/test4/ppMuXLoose/myOutputFile_1.root");
  fwlite::Event ev(file);
  */
  
  vector<string> fileNames;
  setInput(fileNames,sample);
  fwlite::ChainEvent ev(fileNames);

  int counterEvents(0),counterPassedEvents(0);
  int counterCandidates(0),counterPassedCandidates(0);

  int counter(0);
  for( ev.toBegin(); ! ev.atEnd(); ++ev) {
    counter++; 
    counterEvents++;

    if(counter == 20000) break;

    fwlite::Handle<pat::CompositeCandidateCollection > collH;
    collH.getByLabel(ev,"myJPsiAnalysisPAT");

    
    if(counter % 10000 == 1) {   
      cout << "file name: " << ev.getTFile()->GetName() << endl;
    }
    

    bool passedEvent(false);
    for(vector<pat::CompositeCandidate>::const_iterator it=collH.ptr()->begin();
	it!=collH.ptr()->end();++it){
      //if(it->recoJPsi().mass() >0.) counterCandidates++;

      bool step1(true),step2(true),step3(true),step4(true),step5(true),step6(true);
      
      /*
      if(it->recoLeg(1).charge() == it->recoLeg(2).charge()  ) step1=false;
      if(!(it->recoLeg(1).isGlobalMuon() && it->recoLeg(2).isGlobalMuon())) step2=false;
      if(!(it->recoLegTk(1).found() >=11 && it->recoLegTk(2).found() >=11 )) step3=false;
      if( it->recoLeg(1).pt()< 2.5 || it->recoLeg(2).pt()<2.5 ) step4=false;
      if( it->recoVertex().normalizedChi2() > 4.) step5=false;
      if( it->recoJPsi().mass()<2.6 || it->recoJPsi().mass()>3.5) step6=false;
      //if(step >=7 && it->simJPsi().mass()>0.) continue;
      */


      const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(it->daughter("muon1"));
      const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(it->daughter("muon2"));
      
      if(muon1==0 || muon2==0) {cout << "ERROR: dynamic cast failed" << endl; continue;}
      
      if(muon1->charge() == muon2->charge()) step1=false;
      if(!(muon1->isGlobalMuon() && muon2->isGlobalMuon())) step2=false;
      if(!(muon1->innerTrack()->found() >11  && muon2->innerTrack()->found() >=11 )) step3=false;
      if( muon1->innerTrack()->pt()< 2.5 || muon2->innerTrack()->pt()<2.5 ) step4=false;
      if( it->userFloat("vNChi2") > 4.) step5=false;
      if( it->mass()<2.6 || it->mass()>3.5) step6=false;


      recoMass0->Fill(it->mass()); 
      recoMass0b->Fill(it->mass());
      if(step1) {
	recoMass1->Fill(it->mass());
	recoMass1b->Fill(it->mass());
      }
      if(step1 && step2) {
	recoMass2->Fill(it->mass());
	recoMass2b->Fill(it->mass());
      }
      if(step1 && step2 && step3) {
	recoMass3->Fill(it->mass());
	recoMass3b->Fill(it->mass());
      }
      if(step1 && step2 && step3 && step4) {
	recoMass4->Fill(it->mass());
	recoMass4b->Fill(it->mass());
      }
      if(step1 && step2 && step3 && step4 && step5){
	recoMass5->Fill(it->mass());
	recoMass5b->Fill(it->mass());
      }
      if(step1 && step2 && step3 && step4 && step5 && step6) {
	recoMass6->Fill(it->mass());
	recoMass6b->Fill(it->mass());
      }

      if(it->mass() >0.) counterPassedCandidates++;
      if(!passedEvent){
	counterPassedEvents++;
	passedEvent = true;
      }
      

      vtxNChi2->Fill(it->userFloat("vNChi2"));

    }	     
  }// end loop over events 
    

    /*
      cout << "counterEvents: " << counterEvents << endl;
      cout << "counterCandidates: " << counterCandidates << endl;
      cout << endl;
      cout << "counterPassedEvents: " << counterPassedEvents << endl;
      cout << "counterPassedCandidates: " << counterPassedCandidates << endl;
    */

  //int puppa;
  //cout << "here recoMass is!!" << endl;
  /*
  recoMass0->Draw(); gPad->Update(); cin >> puppa;
  recoMass1->Draw(); gPad->Update(); cin >> puppa;
  recoMass2->Draw(); gPad->Update(); cin >> puppa;
  recoMass3->Draw(); gPad->Update(); cin >> puppa;
  recoMass4->Draw(); gPad->Update(); cin >> puppa;
  recoMass5->Draw(); gPad->Update(); cin >> puppa;
  recoMass6->Draw(); gPad->Update(); cin >> puppa;
  */

  /*
  recoMass0b->Draw(); gPad->Update(); cin >> puppa;
  recoMass1b->Draw(); gPad->Update(); cin >> puppa;
  recoMass2b->Draw(); gPad->Update(); cin >> puppa;
  recoMass3b->Draw(); gPad->Update(); cin >> puppa;
  recoMass4b->Draw(); gPad->Update(); cin >> puppa;
  recoMass5b->Draw(); gPad->Update(); cin >> puppa;
  recoMass6b->Draw(); gPad->Update();
  */

  //cin >> puppa ; vtxNChi2->Draw(); gPad->Update(); 
  //gPad->Print("recoMass.pdf"); cin >>puppa;  


  // write histograms on disk
  char *fileOutName0,*fileOutName1,*fileOutName2,*fileOutName3,*fileOutName4,*fileOutName5,*fileOutName6;
  fileOutName0 = "spectrum.0.root";
  fileOutName1 = "spectrum.1.root";
  fileOutName2 = "spectrum.2.root";
  fileOutName3 = "spectrum.3.root";
  fileOutName4 = "spectrum.4.root";
  fileOutName5 = "spectrum.5.root";
  fileOutName6 = "spectrum.6.root";

  TFile output0(fileOutName0,"update");  recoMass0->Write(outHistoName.c_str()); recoMass0b->Write((outHistoName+"b").c_str());
  TFile output1(fileOutName1,"update");  recoMass1->Write(outHistoName.c_str()); recoMass1b->Write((outHistoName+"b").c_str());
  TFile output2(fileOutName2,"update");  recoMass2->Write(outHistoName.c_str()); recoMass2b->Write((outHistoName+"b").c_str());
  TFile output3(fileOutName3,"update");  recoMass3->Write(outHistoName.c_str()); recoMass3b->Write((outHistoName+"b").c_str());
  TFile output4(fileOutName4,"update");  recoMass4->Write(outHistoName.c_str()); recoMass4b->Write((outHistoName+"b").c_str());
  TFile output5(fileOutName5,"update");  recoMass5->Write(outHistoName.c_str()); recoMass5b->Write((outHistoName+"b").c_str());
  TFile output6(fileOutName6,"update");  recoMass6->Write(outHistoName.c_str()); recoMass6b->Write((outHistoName+"b").c_str());

}


 

void setInput(vector<string>& fileNames,int sample){  
  if(sample==0){  
    fileNames.push_back("/home/users/mangano/OctoberExercise/CMSSW_3_1_4/src/OctoberXTracking/JPsiAnalysis/test/myOutputFile.root");
  }else if(sample==1){  
  //ls -lhtr|grep root|awk '{print "     fileNames.push_back(\"/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/"$9"\");"}'
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_16.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_4.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_38.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_3.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_49.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_26.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_59.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_34.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_5.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_23.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_50.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_54.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_39.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_24.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_35.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_55.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_52.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_11.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_43.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_36.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_30.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_13.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_51.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_8.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_33.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_62.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_56.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_9.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_14.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_28.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_47.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_7.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_17.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_48.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_2.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_10.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_22.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_12.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_25.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_44.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_53.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_27.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_1.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_19.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_42.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_31.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_45.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_37.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_18.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_20.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_60.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_29.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_46.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_15.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_58.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_21.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_32.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_57.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/JPsiMuMu/TRK-JPsi_JPsiMuMu-SKIM_test4/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_61.root");
  }  else if(sample==2){  
  //ls -lh|grep root|awk '{print "     fileNames.push_back(\"/hadoop/cms/store/group/tracker-pog/ppMuX/TRK-JPsi_ppMuX-SKIM_test2/b6694ead8530d72ae99d3d30b3d3cb8b/"$9"\");"}'
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/ppMuX/TRK-JPsi_ppMuX-SKIM_test2/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_1.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/ppMuX/TRK-JPsi_ppMuX-SKIM_test2/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_2.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/ppMuX/TRK-JPsi_ppMuX-SKIM_test2/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_3.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/ppMuX/TRK-JPsi_ppMuX-SKIM_test2/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_4.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/ppMuX/TRK-JPsi_ppMuX-SKIM_test2/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_5.root");
  }  else if(sample==3){  
  //ls -lh|grep root|awk '{print "     fileNames.push_back(\"/hadoop/cms/store/group/tracker-pog/ppMuMuX/TRK-JPsi_ppMuMuX-SKIM_test2/b6694ead8530d72ae99d3d30b3d3cb8b/"$9"\");"}'ls -lh|grep root|awk '{print "     fileNames.push_back(\"/hadoop/cms/store/group/tracker-pog/ppMuMuX/TRK-JPsi_ppMuMuX-SKIM_test2/b6694ead8530d72ae99d3d30b3d3cb8b/"$9"\");"}'
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/ppMuMuX/TRK-JPsi_ppMuMuX-SKIM_test2/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_1.root");
     fileNames.push_back("/hadoop/cms/store/group/tracker-pog/ppMuMuX/TRK-JPsi_ppMuMuX-SKIM_test2/b6694ead8530d72ae99d3d30b3d3cb8b/myOutputFile_2.root");
  } else{
    cout << "ERROR: select a different input file" << endl;
  }   
}




