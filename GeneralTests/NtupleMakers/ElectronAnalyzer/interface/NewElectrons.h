#ifndef NewElectrons_H
#define NewElectrons_H

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GlobalCtfElectron.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>
#include <Math/Point3D.h>
#include <vector>
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CLHEP/HepMC/GenParticle.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class TFile;
class TTree;

class NewElectrons: public edm::EDAnalyzer {
public:

  NewElectrons(const edm::ParameterSet& pset);

  virtual ~NewElectrons();

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  void beginJob(const edm::EventSetup& eventSetup);
  void endJob();
  bool inCrack(float eta);
  int mother(HepMC::GenParticle *p); 
  void R9_25_gsf(const edm::Event & event, const reco::PixelMatchGsfElectron*,
                 float&, float&, float&, float&, float&); 
  void R9_25_ctf(const edm::Event & event, const reco::GlobalCtfElectron*,
                 float&, float&, float&, float&, float&);

protected:

private:
  std::string baselineEleCollName;
  std::string customEleCollName;
  std::string fileName;

  TFile *file;
  TTree *tree;
  int nMC;
  float mc_pt[10], mc_eta[10], mc_phi[10], mc_e[10];
  float sc_e[10], sc_eta[10], sc_phi[10], sc_dr[10], sc_et[10];
  float tk_pt[10], tk_eta[10], tk_phi[10], tk_dr[10];
  float el_pt[10], el_eta[10], el_phi[10], el_dr[10], el_e[10];
  float el_eopin[10], el_eopout[10], el_hoe[10], el_detain[10], el_dphiin[10];
  float el_fbrem[10], el_eseed[10], el_e3x3[10], el_detaout[10], el_dphiout[10];
  float el_e5x5[10], el_spp[10], el_see[10], el_pout[10];	
  float el1_pt[10], el1_eta[10], el1_phi[10], el1_dr[10], el1_e[10];
  float el1_eopin[10], el1_eopout[10], el1_hoe[10], el1_detain[10], el1_dphiin[10];
  float el1_fbrem[10], el1_eseed[10], el1_e3x3[10], el1_detaout[10], el1_dphiout[10];
  float el1_e5x5[10], el1_spp[10], el1_see[10], el1_pout[10];	
  int mc_id[10], sc_type[10], mc_mother[10], mc_crack[10];
  int tk_layer[10], tk_subdet[10], tk_nhit[10];
  int el_class[10], el1_class[10];

  HepMC::GenEvent* myGenEvent;
};
#endif

