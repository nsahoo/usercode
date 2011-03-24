//
// $Id: PATJetSlimmer.cc,v 1.2 2010/01/26 12:32:00 gpetrucc Exp $
//

/**
  \class    pat::PATJetSlimmer PATJetSlimmer.h "PhysicsTools/PatAlgos/interface/PATJetSlimmer.h"
  \brief    Matcher of reconstructed objects to L1 Muons 
            
  \author   Giovanni Petrucciani
  \version  $Id: PATJetSlimmer.cc,v 1.2 2010/01/26 12:32:00 gpetrucc Exp $
*/


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#define protected public
#include "DataFormats/PatCandidates/interface/Jet.h"
#undef protected

namespace pat {

  class PATJetSlimmer : public edm::EDProducer {
    public:
      explicit PATJetSlimmer(const edm::ParameterSet & iConfig);
      virtual ~PATJetSlimmer() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      edm::InputTag src_;
      
      /// clear mJetArea, mPassNumber, mPileupEnergy
      bool clearJetVars_;
      /// reset daughters to an empty vector
      bool clearDaughters_;
//       /// reduce GenJet to a bare 4-vector
//       bool slimGenJet_;
      /// drop the Calo or PF specific
      bool dropSpecific_;
//       /// drop the JetCorrFactors (but keep the jet corrected!)
//       bool dropJetCorrFactors_;
  };

} // namespace

pat::PATJetSlimmer::PATJetSlimmer(const edm::ParameterSet & iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    clearJetVars_(iConfig.getParameter<bool>("clearJetVars")),
    clearDaughters_(iConfig.getParameter<bool>("clearDaughters")),
//     slimGenJet_(iConfig.getParameter<bool>("slimGenJet")),
    dropSpecific_(iConfig.getParameter<bool>("dropSpecific"))
//     dropJetCorrFactors_(iConfig.getParameter<bool>("dropJetCorrFactors"))
{
    produces<std::vector<pat::Jet> >();
}

void 
pat::PATJetSlimmer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<pat::Jet> >      src;
    iEvent.getByLabel(src_, src);

    auto_ptr<vector<pat::Jet> >  out(new vector<pat::Jet>());
    out->reserve(src->size());

    for (View<pat::Jet>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
        out->push_back(*it);
        pat::Jet & jet = out->back();

        if (clearJetVars_) {
//             jet.setJetArea(0); 
            jet.setNPasses(0);
//             jet.setPileup(0);
        }
        if (clearDaughters_) {
            jet.clearDaughters();
            jet.pfCandidatesFwdPtr_.clear();
            jet.caloTowersFwdPtr_.clear();
        }
//         if (slimGenJet_) {
//             const reco::GenJet * genjet = it->genJet();
//             if (genjet) {
//                 std::vector<reco::GenJet> tempGenJet(1, reco::GenJet(genjet->p4(), reco::Particle::Point(), reco::GenJet::Specific()));
//                 jet.setGenJet(reco::GenJetRef(&tempGenJet,0), true);
//             }
//         }
        if (dropSpecific_) {
            // FIXME add method in pat::Jet
            jet.specificCalo_.clear();    
            jet.specificPF_.clear();    
        }
//         if (dropJetCorrFactors_) {
//             // FIXME add method in pat::Jet
//             jet.jetEnergyCorrections_.clear();
//         }
    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace pat;
DEFINE_FWK_MODULE(PATJetSlimmer);
