#include "WWAnalysis/DataFormats/interface/SkimEvent.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace
{
  struct dictionary {
    reco::SkimEvent dummy1;
    std::vector<reco::SkimEvent> dummy2;
    edm::Wrapper<reco::SkimEvent> dummy3;
    edm::Wrapper<std::vector<reco::SkimEvent> > dummy4;
  };
}

