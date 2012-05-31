#include <vector>
#include <algorithm>
#include <cmath>

#include "WWAnalysis/DataFormats/interface/SkimEvent2L.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// bunch of anonymous tools
namespace {
    typedef std::pair<edm::ProductID,size_t> edmID;
    template<typename T, typename C, typename F>  edmID id(const edm::Ref<T,C,F> &r) { return edmID(r.id(), r.key()); }
    template<typename T> edmID id(const edm::Ptr<T> &r) { return edmID(r.id(), r.key()); }
    template<typename T> edmID id(const edm::RefToBase<T> &r) { return edmID(r.id(), r.key()); }
    edmID id(const reco::Candidate *c) { return c->hasMasterClonePtr() ? id(c->masterClonePtr()) : id(c->masterClone()); }
    
    struct SortIndexByPt {
        public:
            SortIndexByPt(const edm::View<pat::PFParticle> &pho, bool byPt) : pho_(pho), byPt_(byPt) {}
            bool operator()(int i, int j) const { 
                if (byPt_ || (pho_[i].pt() > 4 || pho_[j].pt() > 4)) {
                    return pho_[i].pt() > pho_[j].pt(); 
                } else {
                    return deltaR(pho_[i], *pho_[i].overlaps("goodLepNoIso")[0]) < deltaR(pho_[j], *pho_[j].overlaps("goodLepNoIso")[0]);
                }
            }
        private:
            const edm::View<pat::PFParticle> &pho_;
            bool  byPt_;
    };
}

class SkimEvent2LFsrCollector : public edm::EDProducer {
    public:
        SkimEvent2LFsrCollector(const edm::ParameterSet &iConfig) ;
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;
    private:
        edm::InputTag zll_, photons_;
        StringCutObjectSelector<pat::PFParticle> photonSelection_;
        std::string photonMatch_;
        std::string isolationLabel_;
        double      isolationCut_;
        bool doFsrRecovery_;
        //bool isZ1_, 
        bool sortPhotonsByPt_;

        bool passIsoFSR(const reco::Candidate &lep, const pat::PFParticle &pho) const ;
        bool passIsoNoFSR(const reco::Candidate &lep) const ;
        bool passIsoFSR(const pat::Muon &lep, const pat::PFParticle &pho) const ;
        bool passIsoNoFSR(const pat::Muon &lep) const ;
        bool passIsoFSR(const pat::Electron &lep, const pat::PFParticle &pho) const ;
        bool passIsoNoFSR(const pat::Electron &lep) const ;
};

SkimEvent2LFsrCollector::SkimEvent2LFsrCollector(const edm::ParameterSet &iConfig) :
    zll_(iConfig.getParameter<edm::InputTag>("zll")),
    photons_(iConfig.getParameter<edm::InputTag>("photons")),
    photonSelection_(iConfig.getParameter<std::string>("photonSelection")),
    photonMatch_(iConfig.getParameter<std::string>("photonMatch")),
    isolationLabel_(iConfig.getParameter<std::string>("isolationLabel")),
    isolationCut_(iConfig.getParameter<double>("isolationCut")),
    doFsrRecovery_(iConfig.getParameter<bool>("doFsrRecovery")),
    //isZ1_(iConfig.getParameter<bool>("isZ1")),
    sortPhotonsByPt_(iConfig.getParameter<bool>("sortPhotonsByPt"))
{
    produces<std::vector<reco::SkimEvent2L> >();
}

void
SkimEvent2LFsrCollector::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<std::vector<reco::SkimEvent2L> > zll; 
    iEvent.getByLabel(zll_, zll);
    edm::Handle<edm::View<pat::PFParticle> > photons; 
    if (doFsrRecovery_) iEvent.getByLabel(photons_, photons);

    std::auto_ptr<std::vector<reco::SkimEvent2L> > out(new std::vector<reco::SkimEvent2L>());

    // Associate to each lepton a list of photons which have that lepton as closest one
    std::map<edmID,std::vector<int> > leptonToPhotonMap; int idx = 0;
    if (doFsrRecovery_) {
        for (edm::View<pat::PFParticle>::const_iterator it = photons->begin(), ed = photons->end(); it != ed; ++it, ++idx) {
            if (!photonSelection_(*it)) continue;
            const reco::CandidatePtrVector &leps = it->overlaps(photonMatch_);
            if (leps.empty()) continue;
            leptonToPhotonMap[id(leps[0])].push_back(idx);
        }        
    }
    // Now loop on Z candidates
    foreach(const reco::SkimEvent2L &z0, *zll) {
        // Retrieve FSR photons, make a list
        std::vector<int> myphotons;
        if (doFsrRecovery_)  {
            foreach(const int &i, leptonToPhotonMap[id(z0.daughter(0))]) myphotons.push_back(i);
            foreach(const int &i, leptonToPhotonMap[id(z0.daughter(1))]) myphotons.push_back(i);
            // now sort them by pt
            std::sort(myphotons.begin(), myphotons.end(), SortIndexByPt(*photons, sortPhotonsByPt_));
        }
        // now loop and pick the first one that works, if any
        const pat::PFParticle *match = 0; int imatch = -1;
        foreach(const int &i, myphotons) {
            const pat::PFParticle *pho = &photons->at(i);
            // Kinematic cuts
            reco::Particle::LorentzVector p4ll = z0.p4(), p4llg = p4ll + pho->p4();
            double mll = p4ll.M(), mllg = p4llg.M();
            if (mllg > 100 || mllg <  4) continue;
            if (fabs(mllg - 91.188) >= fabs(mll - 91.188)) continue;
            match = pho; imatch = i; break;
        }
        // Now test isolation
        bool passIso = true;
        if (match) {
            int inear = deltaR2(*z0.daughter(0), *match) <=  deltaR2(*z0.daughter(1), *match) ? 0 : 1;
            passIso = passIsoFSR(z0.l(inear), *match) && passIsoNoFSR(z0.l(1-inear));
        } else {
            passIso = passIsoNoFSR(z0.l(0)) && passIsoNoFSR(z0.l(1));
        }
        if (passIso) {
            out->push_back(z0);
            if (match) {
                reco::ShallowClonePtrCandidate clone(photons->ptrAt(imatch), 0, match->p4(), match->vertex());
                clone.setPdgId(22); clone.setStatus(match->status());
                out->back().addDaughter(clone);
                out->back().setP4(out->back().p4() + clone.p4());
            }
        } 
    }

    iEvent.put(out);
}

bool SkimEvent2LFsrCollector::passIsoFSR(const reco::Candidate &lep, const pat::PFParticle &pho) const {
    return abs(lep.pdgId()) == 13 ? 
                passIsoFSR(dynamic_cast<const pat::Muon     &>(lep), pho) :
                passIsoFSR(dynamic_cast<const pat::Electron &>(lep), pho);
}
bool SkimEvent2LFsrCollector::passIsoNoFSR(const reco::Candidate &lep) const {
    return abs(lep.pdgId()) == 13 ? 
                passIsoNoFSR(dynamic_cast<const pat::Muon     &>(lep)) :
                passIsoNoFSR(dynamic_cast<const pat::Electron &>(lep));
}
bool SkimEvent2LFsrCollector::passIsoFSR(const pat::Muon &lep, const pat::PFParticle &pho) const {
    double dr = deltaR(lep,pho);
    if (dr < 0.4 && dr >= 0.01 && pho.pt() > 0.5) {
        float chiso = lep.userFloat(isolationLabel_+"ChHad");
        float nhiso = lep.userFloat(isolationLabel_+"NHad");
        float phiso = lep.userFloat(isolationLabel_+"Pho");
        float eatot = lep.userFloat(isolationLabel_+"EAtot");
        float rho   = lep.userFloat(isolationLabel_+"Rho");
        if (phiso < pho.pt()) throw cms::Exception("LogicError") << "Cannot subtract photon of pt " << pho.pt() << " from photon iso " << phiso << ", for muon of pt " << lep.pt() << ", eta " << lep.eta() << ", phi " << lep.phi() << ", deltaR = " << dr << std::endl;
        phiso -= pho.pt();
        double iso = chiso + std::max<double>(0., nhiso + phiso - rho*eatot);
        return iso/lep.pt() < isolationCut_;
    } else {
        return passIsoNoFSR(lep); 
    }

}
bool SkimEvent2LFsrCollector::passIsoNoFSR(const pat::Muon &lep) const {
    return lep.userFloat(isolationLabel_)/lep.pt() < isolationCut_;
}

bool SkimEvent2LFsrCollector::passIsoFSR(const pat::Electron &lep, const pat::PFParticle &pho) const {
    double dr = deltaR(lep,pho);
    if (dr < 0.4 && (fabs(lep.superCluster()->eta()) < 1.479 || dr > 0.08)) {
        float chiso = lep.userFloat(isolationLabel_+"ChHad");
        float nhiso = lep.userFloat(isolationLabel_+"NHad");
        float phiso = lep.userFloat(isolationLabel_+"Pho");
        float eatot = lep.userFloat(isolationLabel_+"EAtot");
        float rho   = lep.userFloat(isolationLabel_+"Rho");
        if (phiso < pho.pt()) throw cms::Exception("LogicError") << "Cannot subtract photon of pt " << pho.pt() << " from photon iso " << phiso << ", for electron of pt " << lep.pt() << ", eta " << lep.eta() << ", phi " << lep.phi() << ", deltaR = " << dr << std::endl;
        phiso -= pho.pt();
        double iso = chiso + std::max<double>(0., nhiso + phiso - rho*eatot);
        return iso/lep.pt() < isolationCut_;
    } else {
        return passIsoNoFSR(lep); 
    }
}
bool SkimEvent2LFsrCollector::passIsoNoFSR(const pat::Electron &lep) const {
    return lep.userFloat(isolationLabel_)/lep.pt() < isolationCut_;
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SkimEvent2LFsrCollector);
