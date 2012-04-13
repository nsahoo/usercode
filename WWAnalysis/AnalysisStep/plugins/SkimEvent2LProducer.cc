#include <vector>
#include <algorithm>
#include <cmath>

#include "WWAnalysis/DataFormats/interface/SkimEvent2L.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"


class SkimEvent2LProducer : public edm::EDProducer {
    public:
        SkimEvent2LProducer(const edm::ParameterSet &iConfig) ;
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;
    private:
        edm::InputTag src_;
        edm::InputTag pfMet_;
        edm::InputTag vertices_;
        bool isMC;
};

SkimEvent2LProducer::SkimEvent2LProducer(const edm::ParameterSet &iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    pfMet_(iConfig.getParameter<edm::InputTag>("pfMet")),
    vertices_(iConfig.getParameter<edm::InputTag>("vertices")),
    isMC(iConfig.existsAs<bool>("isMC")?iConfig.getParameter<bool>("isMC"):false)
{
    produces<std::vector<reco::SkimEvent2L> >();
}

void
SkimEvent2LProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<reco::CandidateView> src; 
    iEvent.getByLabel(src_, src);

    edm::Handle<reco::VertexCollection> vertices; 
    iEvent.getByLabel(vertices_, vertices);
    edm::Handle<reco::PFMETCollection> pfMet; 
    iEvent.getByLabel(pfMet_, pfMet);
    edm::Handle<std::vector<PileupSummaryInfo> > puH;
    if (isMC) iEvent.getByLabel("addPileupInfo", puH);

    std::auto_ptr<std::vector<reco::SkimEvent2L> > out(new std::vector<reco::SkimEvent2L>());
    for (reco::CandidateView::const_iterator it = src->begin(), ed= src->end(); it != ed; ++it) {
        const reco::CompositeCandidate *srcz = dynamic_cast<const reco::CompositeCandidate *>(&*it);
        if (srcz == 0) throw cms::Exception("CorruptData") << "Found something that is not a composite candidate" << std::endl;
        const pat::CompositeCandidate *patz = dynamic_cast<const pat::CompositeCandidate *>(srcz);
        if (patz == 0) {
            const reco::Candidate *d0 = srcz->daughter(0);
            if (d0 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)\n";
            out->push_back( reco::SkimEvent2L(*srcz) );
        } else {
            const reco::Candidate *d0 = patz->daughter(0);
            if (d0 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)\n";
            out->push_back( reco::SkimEvent2L(*patz) );
        }
        reco::SkimEvent2L &z = out->back();
        z.setVertices(vertices);
        z.setPFMet(pfMet);
        if (isMC) z.setPileupInfo(puH);
    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SkimEvent2LProducer);
