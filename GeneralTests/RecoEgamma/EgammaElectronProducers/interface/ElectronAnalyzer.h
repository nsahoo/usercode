#ifndef ElectronAnalyzer_h
#define ElectronAnalyzer_h

// -*- C++ -*-
//
// Package:    ElectronAnalyzer
// Class:      ElectronAnalyzer
// 
/**\class ElectronAnalyzer ElectronAnalyzer.cc Demo/ElectronAnalyzer/src/ElectronAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Alessandro Palma
//         Created:  Thu Sep 21 11:41:35 CEST 2006
// $Id: ElectronAnalyzer.h,v 1.9 2007/03/14 09:13:54 rahatlou Exp $
//
//


// system include files
//#include <memory>
#include<string>
//#include "math.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "FWCore/Utilities/interface/Exception.h"


class TH1F;
class TFile;

//
// class declaration
//

class ElectronAnalyzer : public edm::EDAnalyzer {
 public:
  explicit ElectronAnalyzer(const edm::ParameterSet&);
  ~ElectronAnalyzer();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  double minElePt_;
  double REleCut_;
  TFile*  rootFile_;  
  std::string outputFile_;
  edm::InputTag electronProducer_;
  std::string mcProducer_;
  edm::InputTag scProducer_;

  TH1F* h1_nEleReco_;
  TH1F* h1_recoEleEnergy_;
  TH1F* h1_recoElePt_;
  TH1F* h1_recoEleEta_;
  TH1F* h1_recoElePhi_;
  TH1F* h1_RMin_;
  TH1F* h1_eleERecoOverEtrue_;
  TH1F* h1_eleRecoTrackChi2_;
  TH1F* h1_recoElePtRes_;
  TH1F* h1_recoEleDeltaEta_;
  TH1F* h1_recoEleDeltaPhi_;

  std::string islandBarrelBasicClusterCollection_;
  std::string islandBarrelBasicClusterProducer_;
  std::string islandBarrelBasicClusterShapes_;

  TH1F* h1_islandEBBC_e2x2_Over_e3x3_;
  TH1F* h1_islandEBBC_e3x3_Over_e5x5_;


};

#endif
