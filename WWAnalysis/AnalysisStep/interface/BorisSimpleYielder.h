#include <map>
#include <string>

#include "TH1.h"
#include "PhysicsTools/UtilAlgos/interface/BasicAnalyzer.h"
#include "WWAnalysis/AnalysisStep/interface/TriggerBitChecker.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent.h"

const int numberCuts(20);

class BorisSimpleYielder : public edm::BasicAnalyzer {
  
 public:
  /// default constructor
  BorisSimpleYielder(const edm::ParameterSet& cfg, TFileDirectory& fs);
  /// default destructor
  virtual ~BorisSimpleYielder();
  /// everything that needs to be done before the event loop
  void beginJob(){};
  /// everything that needs to be done after the event loop
  void endJob(){};
  /// everything that needs to be done during the event loop
  void analyze(const edm::EventBase& event);

 private:
  void printout(const edm::EventBase& event,std::vector<reco::SkimEvent>::const_iterator it);

  // ----- datamembers
 private:  
  edm::InputTag eventHypo0_;
  edm::InputTag eventHypoID_;
  edm::InputTag eventHypoISO_;
  edm::InputTag eventHypoNOCONV_;
  edm::InputTag eventHypoIP_;
  double minMet;
  double minMll;
  double minDiffMz;
  int  nCentralJet;
  double jetVetoEt;
  double jetVetoEta;
  double maxMll;
  double ptMinHigh;
  double ptMinLow;
  double minProjMet;
  double maxDPhiLLDegrees;
  bool printEventNumber;
  /// input tag for mouns
  //edm::InputTag skimEventsTag_;
  //int pdType_;
  int yieldCounter[numberCuts];
  long int lastEventBookkeeper[numberCuts];

  /*
  bool cut(std::vector<reco::SkimEvent>::const_iterator it, 
	   reco::SkimEvent::primaryDatasetType pdType,
	   edm::TriggerResults triggerResults,
	   const edm::EventBase& event);
  */

  /// histograms
  std::map<std::string, TH1*> hists_;
};
