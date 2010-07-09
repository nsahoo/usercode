#ifndef NewElectrons_H
#define NewElectrons_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "NtupleMakers/ElectronAnalyzer/interface/UserData.h"



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

protected:

private:
  std::string baselineEleCollName;
  std::string customEleCollName;


  TFile *file;
  TTree *tree;
  int nMC;
  float mc_pt[10], mc_eta[10], mc_phi[10], mc_e[10];
  float sc_e[10], sc_eta[10], sc_phi[10], sc_dr[10], sc_et[10];
  float tk_pt[10], tk_eta[10], tk_phi[10], tk_dr[10];
  float el_pt[10], el_eta[10], el_phi[10], el_dr[10], el_e[10];
  float el1_pt[10], el1_eta[10], el1_phi[10], el1_dr[10], el1_e[10];
  int mc_id[10], sc_type[10], mc_mother[10], mc_crack[10];
  int tk_layer[10], tk_subdet[10], tk_nhit[10];

  HepMC::GenEvent* myGenEvent;

  //UserDataInt1D idMC;
  std::vector<math::XYZTLorentzVector> mc_p4_vector;
};
#endif

