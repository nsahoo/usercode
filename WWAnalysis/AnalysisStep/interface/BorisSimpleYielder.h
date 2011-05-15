#include <map>
#include <string>

#include "TH1.h"
#include "PhysicsTools/UtilAlgos/interface/BasicAnalyzer.h"
#include "WWAnalysis/AnalysisStep/interface/TriggerBitChecker.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent.h"

class BorisSimpleYielder : public edm::BasicAnalyzer {

 public:
  /// default constructor
  BorisSimpleYielder(const edm::ParameterSet& cfg, TFileDirectory& fs);
  /// default destructor
  virtual ~BorisSimpleYielder(){};
  /// everything that needs to be done before the event loop
  void beginJob(){};
  /// everything that needs to be done after the event loop
  void endJob(){};
  /// everything that needs to be done during the event loop
  void analyze(const edm::EventBase& event);

 private:
  /// input tag for mouns
  edm::InputTag skimEventsTag_;
  int pdType_;
  void printout(const edm::EventBase& event,std::vector<reco::SkimEvent>::const_iterator it);
  bool cut(std::vector<reco::SkimEvent>::const_iterator it, 
	   reco::SkimEvent::primaryDatasetType pdType,
	   edm::TriggerResults triggerResults,
	   const edm::EventBase& event);


  /// histograms
  std::map<std::string, TH1*> hists_;
};
