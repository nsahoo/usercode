#include <map>
#include <string>

#include "TH1.h"
#include "PhysicsTools/UtilAlgos/interface/BasicAnalyzer.h"


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
  
  /// histograms
  std::map<std::string, TH1*> hists_;
};
