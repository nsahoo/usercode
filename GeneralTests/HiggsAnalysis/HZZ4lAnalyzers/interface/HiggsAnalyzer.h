// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"


// Track Specifc include Files
//#include "DataFormats/TrackReco/interface/Track.h"

// root specific include files
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

// analyzer specific class
#include "HiggsAnalysis/DataFormats/interface/HiggsEvent.h"

class HiggsAnalyzer : public edm::EDAnalyzer
{
 public:
  explicit HiggsAnalyzer(const edm::ParameterSet&);
  ~HiggsAnalyzer();


 private:
  virtual void beginJob(const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  HepMC::GenParticle* GenParticleMother(HepMC::GenParticle* MCParticle);
  HepMC::GenParticle* GenParticleGrandMother(HepMC::GenParticle* MCParticle);
  
  HiggsEvent*	higgsEvent;
  
  TFile*	outputRootFile;
  TTree*	higgsEventTree;
  //TH1I	*hNumEvents;
  int		numEventsProcessed;
		
  std::string outputRootFileName;
  std::string  hepMCLabel;
  
};

//define this as a plug-in
DEFINE_FWK_MODULE(HiggsAnalyzer);
