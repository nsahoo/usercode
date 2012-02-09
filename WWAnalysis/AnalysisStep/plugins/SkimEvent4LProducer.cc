#include <vector>
#include <algorithm>
#include <cmath>

#include "WWAnalysis/DataFormats/interface/SkimEvent4L.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"


class SkimEvent4LProducer : public edm::EDProducer {
    public:
        SkimEvent4LProducer(const edm::ParameterSet &iConfig) ;
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;
    private:
        edm::InputTag src_;

        edm::InputTag jets_;
        edm::InputTag pfMet_;
        edm::InputTag vertices_;

        bool          isMC_;
        edm::InputTag mcMatch_;
};

SkimEvent4LProducer::SkimEvent4LProducer(const edm::ParameterSet &iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    jets_(iConfig.getParameter<edm::InputTag>("jets")),
    pfMet_(iConfig.getParameter<edm::InputTag>("pfMet")),
    vertices_(iConfig.getParameter<edm::InputTag>("vertices")),
    isMC_(iConfig.getParameter<bool>("isMC")),
    mcMatch_(isMC_ ? iConfig.getParameter<edm::InputTag>("mcMatch") : edm::InputTag("FAKE"))
{
    produces<std::vector<reco::SkimEvent4L> >();
}

void
SkimEvent4LProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<reco::CandidateView> src; 
    iEvent.getByLabel(src_, src);

    edm::Handle<reco::VertexCollection> vertices; 
    iEvent.getByLabel(vertices_, vertices);
    edm::Handle<pat::JetCollection> jets; 
    iEvent.getByLabel(jets_, jets);
    edm::Handle<reco::PFMETCollection> pfMet; 
    iEvent.getByLabel(pfMet_, pfMet);

    edm::Handle<edm::Association<reco::GenParticleCollection> > mcMatch; 
    if (isMC_) iEvent.getByLabel(mcMatch_, mcMatch);

    std::auto_ptr<std::vector<reco::SkimEvent4L> > out(new std::vector<reco::SkimEvent4L>());
    for (reco::CandidateView::const_iterator it = src->begin(), ed= src->end(); it != ed; ++it) {
        const reco::CompositeCandidate *srczz = dynamic_cast<const reco::CompositeCandidate *>(&*it);
        if (srczz == 0) throw cms::Exception("CorruptData") << "Found something that is not a composite candidate" << std::endl;

        out->push_back( reco::SkimEvent4L(*srczz) );
        reco::SkimEvent4L &zz = out->back();

        zz.setVertex(vertices);
        zz.setPFMet(pfMet);
        zz.setJets(jets);

        if (isMC_) zz.setGenMatches(*mcMatch);
    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SkimEvent4LProducer);
