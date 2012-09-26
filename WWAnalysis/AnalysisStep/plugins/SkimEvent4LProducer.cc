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
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "WWAnalysis/AnalysisStep/interface/HZZ4lMelaDiscriminator.h"
#include "WWAnalysis/AnalysisStep/interface/CompositeCandMassResolution.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

class SkimEvent4LProducer : public edm::EDProducer {
    public:
        SkimEvent4LProducer(const edm::ParameterSet &iConfig) ;
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;
    private:
        edm::InputTag src_;

        edm::InputTag reducedPFCands_;
        edm::InputTag jets_;
        edm::InputTag pfMet_;
        edm::InputTag vertices_;
        edm::InputTag rho_;

        bool          isMC_;
        bool          isSignal_;
        bool doswap;
        edm::InputTag mcMatch_;
        bool doMELA_;
        bool doAnglesWithFSR_;
        bool doMassRes_;

        std::auto_ptr<HZZ4LMelaDiscriminator> melaSMH_, melaPSH_, melaQQZZ_;
        CompositeCandMassResolution massRes_;
};

SkimEvent4LProducer::SkimEvent4LProducer(const edm::ParameterSet &iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    reducedPFCands_(iConfig.getParameter<edm::InputTag>("reducedPFCands")),
    jets_(iConfig.getParameter<edm::InputTag>("jets")),
    pfMet_(iConfig.getParameter<edm::InputTag>("pfMet")),
    vertices_(iConfig.getParameter<edm::InputTag>("vertices")),
    isMC_(iConfig.getParameter<bool>("isMC")),
    isSignal_(iConfig.existsAs<bool>("isSignal")?iConfig.getParameter<bool>("isSignal"):false),
    doswap(iConfig.existsAs<bool>("doswap")?iConfig.getParameter<bool>("doswap"):true),
    mcMatch_(isSignal_ ? iConfig.getParameter<edm::InputTag>("mcMatch") : edm::InputTag("FAKE")),
    doMELA_(iConfig.existsAs<bool>("doMELA")?iConfig.getParameter<bool>("doMELA"):false),
    doAnglesWithFSR_(iConfig.existsAs<bool>("doAnglesWithFSR")?iConfig.getParameter<bool>("doAnglesWithFSR"):true),
    doMassRes_(iConfig.existsAs<bool>("doMassRes")?iConfig.getParameter<bool>("doMassRes"):false)
{
    if (doMELA_) {
        std::string spath = edm::FileInPath(iConfig.getParameter<std::string>("melaQQZZHistos")).fullPath();
        const char *cpath = spath.c_str();
        melaSMH_.reset(new HZZ4LMelaDiscriminator(HZZ4LMelaDiscriminator::SMHiggs, NULL));
        melaPSH_.reset(new HZZ4LMelaDiscriminator(HZZ4LMelaDiscriminator::PSHiggs, NULL));
        melaQQZZ_.reset(new HZZ4LMelaDiscriminator(HZZ4LMelaDiscriminator::QQZZ,   cpath));
    }
    produces<std::vector<reco::SkimEvent4L> >();
}

void
SkimEvent4LProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    if (doMassRes_) massRes_.init(iSetup);

    edm::Handle<reco::CandidateView> src; 
    iEvent.getByLabel(src_, src);

    edm::Handle<reco::VertexCollection> vertices; 
    iEvent.getByLabel(vertices_, vertices);
    edm::Handle<pat::JetCollection> jets; 
    iEvent.getByLabel(jets_, jets);
    edm::Handle<reco::PFMETCollection> pfMet; 
    iEvent.getByLabel(pfMet_, pfMet);
    edm::Handle<std::vector<reco::LeafCandidate> > pfleaves; 
    iEvent.getByLabel(reducedPFCands_, pfleaves);

    edm::Handle<edm::Association<reco::GenParticleCollection> > mcMatch; 
    if (isSignal_) iEvent.getByLabel(mcMatch_, mcMatch);

    edm::Handle<std::vector<PileupSummaryInfo> > puH; 
    if (isMC_) iEvent.getByLabel("addPileupInfo", puH);

    std::auto_ptr<std::vector<reco::SkimEvent4L> > out(new std::vector<reco::SkimEvent4L>());
    for (reco::CandidateView::const_iterator it = src->begin(), ed= src->end(); it != ed; ++it) {
        const reco::CompositeCandidate *srczz = dynamic_cast<const reco::CompositeCandidate *>(&*it);
        if (srczz == 0) throw cms::Exception("CorruptData") << "Found something that is not a composite candidate" << std::endl;
        const pat::CompositeCandidate *patzz = dynamic_cast<const pat::CompositeCandidate *>(srczz);
        if (patzz == 0) {
            const reco::Candidate *d0 = srczz->daughter(0);
            if (d0 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)\n";
            const reco::Candidate *d00 = srczz->daughter(0)->daughter(0);
            if (d00 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)->daughter(0)\n";
            out->push_back( reco::SkimEvent4L(*srczz, doswap) );
        } else {
            const reco::Candidate *d0 = patzz->daughter(0);
            if (d0 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)\n";
            const reco::Candidate *d00 = patzz->daughter(0)->daughter(0);
            if (d00 == 0) throw cms::Exception("CorruptData") << "Null daughter(0)->daughter(0)\n";
            out->push_back( reco::SkimEvent4L(*patzz, doswap) );
        }
        reco::SkimEvent4L &zz = out->back();

        zz.setVertex(vertices);
        zz.setPFMet(pfMet);
        zz.setJets(jets);
        zz.setPFLeaves(pfleaves);
        zz.setNumRecoVertices(vertices);
        zz.setAngles(doAnglesWithFSR_);
        if (doMELA_) {
            zz.addUserFloat("melaSMH",  melaSMH_->get( zz.mass(), zz.mz(0), zz.mz(1), zz.getCosThetaStar(), zz.getCosTheta1(), zz.getCosTheta2(), zz.getPhi(), zz.getPhi1()));
            zz.addUserFloat("melaPSH",  melaPSH_->get( zz.mass(), zz.mz(0), zz.mz(1), zz.getCosThetaStar(), zz.getCosTheta1(), zz.getCosTheta2(), zz.getPhi(), zz.getPhi1()));
            zz.addUserFloat("melaQQZZ", melaQQZZ_->get(zz.mass(), zz.mz(0), zz.mz(1), zz.getCosThetaStar(), zz.getCosTheta1(), zz.getCosTheta2(), zz.getPhi(), zz.getPhi1()));
        }
        if (doMassRes_) zz.addUserFloat("massErr", massRes_.getMassResolution(zz));
        if (isMC_) zz.setPileupInfo(*puH);
        if (isSignal_) zz.setGenMatches(*mcMatch);
    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SkimEvent4LProducer);
