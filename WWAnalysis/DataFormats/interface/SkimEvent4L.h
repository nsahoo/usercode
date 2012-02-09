#ifndef AnalysisDataFormats_SkimEvent4L_h
#define AnalysisDataFormats_SkimEvent4L_h

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <vector>


namespace utils {
    struct FunctionCache;
    struct SelectorCache;
}

namespace reco {

   class SkimEvent4L : public pat::CompositeCandidate {
        public:
            enum hypoType {undefined = 0, ZZ4EL=3, ZZ4MU=6, ZZ2EL2MU=4, ZZ2MU2EL=5 };

            SkimEvent4L() {}
            virtual ~SkimEvent4L() ;
            explicit SkimEvent4L(const pat::CompositeCandidate &src) ;
            explicit SkimEvent4L(const reco::CompositeCandidate &src) ;

            const reco::Candidate & z(unsigned int i) const { return *daughter(i); }

            const reco::Candidate & l(unsigned int iz, unsigned int il) const { 
                const reco::Candidate &proxy = * daughter(iz)->daughter(il);
                return  proxy.hasMasterClonePtr() ? * proxy.masterClonePtr() :
                       (proxy.hasMasterClone()    ? * proxy.masterClone() :
                        proxy);  
            }

            const reco::Vertex & vtx() const { return *vtx_; }
            const reco::MET    & pfMet() const { return *pfMet_; }
            const int            njetany() const { return jets_.size(); }
            const pat::Jet     & anyjet(unsigned int i) const { return *jets_[i]; }

            const reco::Candidate & zByMass(unsigned int rank) const ;
            const reco::Candidate & lByPt(unsigned int rank) const ;

            const int hypo() const { return hypo_; }
            const int channel() const ;

            float mz(unsigned int iz) const { return z(iz).mass(); } 
            reco::Particle::LorentzVector p4ll(unsigned int iz1, unsigned int il1, unsigned int iz2, unsigned int il2) const {
                return (daughter(iz1)->daughter(il1)->p4() + daughter(iz2)->daughter(il2)->p4());
            }
            float mll(unsigned int iz1, unsigned int il1, unsigned int iz2, unsigned int il2) const { 
                return p4ll(iz1, il1, iz2, il2).mass();
            } 

            float lpt(unsigned int iz, unsigned int il)  const { return l(iz,il).pt();  }
            float leta(unsigned int iz, unsigned int il) const { return l(iz,il).eta(); }
            float lphi(unsigned int iz, unsigned int il) const { return l(iz,il).phi(); }
            float lq(unsigned int iz, unsigned int il)   const { return l(iz,il).charge(); }
            float lpdgId(unsigned int iz, unsigned int il)  const { return l(iz,il).pdgId(); }
            float ldz(unsigned int iz, unsigned int il)  const { return luserFloat(iz,il,"dzPV"); }
            float lip2d(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"tip2"); }
            float lip3d(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"ip2"); }
            float lsip2d(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"tip2")/luserFloat(iz,il,"tipErr2"); }
            float lsip3d(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"ip2")/luserFloat(iz,il,"ipErr2"); }
            float lisoPfCh(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"rePFIsoCh"); }
            float lisoPfNeu(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"rePFIsoNeu"); }
            float lisoPfPho(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"rePFIsoPho"); }
            float lisoPfAll(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"rePFIsoCh") + luserFloat(iz,il,"rePFIsoNeu") + luserFloat(iz,il,"rePFIsoPho"); }
            float luserFloat(unsigned int iz, unsigned int il, const char *label)  const ;
            float luserFloat(unsigned int iz, unsigned int il, const std::string &label)  const ;
            bool  lgood(unsigned int iz, unsigned int il, const char *muId, const char *eleId) const ;
            bool  lgood(unsigned int iz, unsigned int il, const std::string &muId, const std::string &eleId) const ;
            float lval(unsigned int iz, unsigned int il, const char *expr) const ;
            float lval(unsigned int iz, unsigned int il, const std::string &expr) const ;
            float lval(unsigned int iz, unsigned int il, const char *muExpr, const char *eleExpr) const ;
            float lval(unsigned int iz, unsigned int il, const std::string &muExpr, const std::string &eleExpr) const ;

            reco::GenParticleRef genz(int iz) const ;
            reco::GenParticleRef genl(int iz, int il) const ;

            const float elePtMin() const ;
            const float muPtMin() const ;
            const float lPtMin() const ;
            const float eleEtaMax() const ;
            const float muEtaMax() const ;
            const float lEtaMax() const ;
            const int nGoodLeptons(const char *muId, const char *eleId) const ;
            const int nGoodLeptons(const std::string &muId, const std::string &eleId) const ;
            const int nMatchedLeptons() const { return nMatchedLeptonsInZ(0) + nMatchedLeptonsInZ(1); }
            const int nMatchedLeptonsInZ(int iz) const { return genl(iz,0).isNonnull() + genl(iz,1).isNonnull(); }

            bool hasCorrectAssignment() const {
                reco::GenParticleRef z1 = genz(0), z2 = genz(1);
                return (z1.isNonnull() && z2.isNonnull() && z1 != z2);
            }

            using reco::LeafCandidate::setVertex;
            void setVertex(const edm::Handle<reco::VertexCollection> &);
            void setPFMet(const edm::Handle<reco::PFMETCollection> &);
            void setJets(const edm::Handle<pat::JetCollection> &, double ptMin=-1);

            void setGenMatches(const edm::Association<reco::GenParticleCollection> &genMatch) ;
        private:
            void init() ;

            hypoType hypo_;

            reco::VertexRef  vtx_;
            reco::PFMETRef   pfMet_;
            pat::JetRefVector jets_;

            std::vector<reco::GenParticleRef> z1GenMatches_, z2GenMatches_;

            bool passesSingleMuData_;
            bool passesSingleElData_;
            bool passesDoubleMuData_;
            bool passesDoubleElData_;
            bool passesMuEGData_    ;
            bool passesSingleMuMC_  ;
            bool passesSingleElMC_  ;
            bool passesDoubleMuMC_  ;
            bool passesDoubleElMC_  ;
            bool passesMuEGMC_      ;
    };

}

#endif

