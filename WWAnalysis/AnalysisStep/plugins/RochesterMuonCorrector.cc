//
// $Id: MuonCleanerBySegments.cc,v 1.1 2012/09/25 22:07:27 gpetrucc Exp $
//

/**
  \class    modules::RochesterMuonCorrectorT RochesterMuonCorrectorT.h "MuonAnalysis/MuonAssociators/interface/RochesterMuonCorrectorT.h"
  \brief    Applies Rochester corrections to muons            
  \author   Giovanni Petrucciani
  \version  $Id: MuonCleanerBySegments.cc,v 1.1 2012/09/25 22:07:27 gpetrucc Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "WWAnalysis/AnalysisStep/src/rochcor.h"

namespace modules {

  template<typename T>
  class RochesterMuonCorrectorT : public edm::EDProducer {
    public:
      explicit RochesterMuonCorrectorT(const edm::ParameterSet & iConfig);
      virtual ~RochesterMuonCorrectorT() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      /// Labels for input collections
      edm::InputTag src_;

      /// Rochester corrector
      rochcor corrector_;
  };

} // namespace

template<typename T>
modules::RochesterMuonCorrectorT<T>::RochesterMuonCorrectorT(const edm::ParameterSet & iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src"))
{
    produces<std::vector<T> >(); 
}

template<typename T>
void 
modules::RochesterMuonCorrectorT<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<T> > src;
    iEvent.getByLabel(src_, src);

    unsigned int nsrc = src->size();
    auto_ptr<vector<T> > out(new vector<T>());
    out->reserve(nsrc);

    unsigned int run = iEvent.id().run(); 

    TLorentzVector p4;
    for (unsigned int i = 0; i < nsrc; ++i) {
        const T &mu = (*src)[i];
        p4.SetPtEtaPhiM(mu.pt(), mu.eta(), mu.phi(), mu.mass());
        if (run == 1) {
            corrector_.momcor_mc(p4, mu.charge(), 0.0, 0);
        } else {
            corrector_.momcor_data(p4, mu.charge(), 0.0, run <= 173692 ? 0 : 1);
        }
        out->push_back(mu);
        out->back().setP4(reco::Particle::PolarLorentzVector(p4.Pt(), p4.Eta(), p4.Phi(), p4.M()));
    }

    iEvent.put(out);
}


namespace modules {
    //typedef modules::RochesterMuonCorrectorT<reco::Muon>  RochesterMuonCorrector;
    typedef modules::RochesterMuonCorrectorT<pat::Muon>   RochesterPATMuonCorrector;
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace modules;
//DEFINE_FWK_MODULE(RochesterMuonCorrector);
DEFINE_FWK_MODULE(RochesterPATMuonCorrector);
