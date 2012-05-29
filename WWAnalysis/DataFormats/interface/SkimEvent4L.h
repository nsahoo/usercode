#ifndef AnalysisDataFormats_SkimEvent4L_h
#define AnalysisDataFormats_SkimEvent4L_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "WWAnalysis/DataFormats/interface/hzziso.h"

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
            explicit SkimEvent4L(const pat::CompositeCandidate &src, bool doswap = true) ;
            explicit SkimEvent4L(const reco::CompositeCandidate &src, bool doswap = true) ;

            const reco::Candidate & z(unsigned int i) const { return *daughter(i); }

            const reco::Candidate & l(unsigned int iz, unsigned int il) const { 
                const reco::Candidate &proxy = * lproxy(iz,il);
                return  proxy.hasMasterClonePtr() ? * proxy.masterClonePtr() :
                       (proxy.hasMasterClone()    ? * proxy.masterClone() :
                        proxy);  
            }

            class Pair : public reco::LeafCandidate {
                public:
                    Pair() {}
                    Pair(const Candidate *c1, const Candidate *c2) {
                        addDaughter(c1);
                        addDaughter(c2);
                    }
                    virtual size_t numberOfDaughters() const { return daughters_.size(); }
                    virtual const Candidate * daughter( size_type i ) const { return (i < daughters_.size() ? daughters_[i] : 0); }
                    float deltaR(size_type i1=0, size_type i2=1) const {
                        return reco::deltaR(*daughters_[i1], *daughters_[i2]);
                    }
                    float deltaPhi(size_type i1=0, size_type i2=1) const {
                        return reco::deltaPhi(daughters_[i1]->phi(), daughters_[i2]->phi());
                    }
                    void addDaughter(const reco::Candidate *c) {
                        daughters_.push_back(c);
                        setP4(p4()+c->p4());
                        setCharge(charge()+c->charge());
                    }
                    void clear() {
                        daughters_.clear();
                        setP4(reco::Particle::LorentzVector(0,0,0,0));
                        setCharge(0);
                    }
                private:
                    std::vector<const reco::Candidate *> daughters_;
            };
            const Pair pair(unsigned int iz1, unsigned int il1, unsigned int iz2, unsigned int il2) const {
                return Pair( & l(iz1,il1), & l(iz2,il2) );
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
                return (lproxy(iz1,il1)->p4() + lproxy(iz2,il2)->p4());
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
            float lip2d(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"tip"); }
            float lip3d(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"ip"); }
            float lsip2d(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"tip")/luserFloat(iz,il,"tipErr"); }
            float lsip3d(unsigned int iz, unsigned int il)  const  { return luserFloat(iz,il,"ip")/luserFloat(iz,il,"ipErr"); }

            float lisoTrkCustom(unsigned int iz, unsigned int il) const ;
            float lisoNeuCustom(unsigned int iz, unsigned int il) const ;
            float lisoPhoCustom(unsigned int iz, unsigned int il) const ;
            float lisoTrkBaseline(unsigned int iz, unsigned int il) const { return luserFloat(iz,il,"tkZZ4L"); }
            float lisoCombinedCorrCustom(unsigned int iz, unsigned int il, float deltaR, int type);

            bool lfiresTrigger(unsigned int iz, unsigned int il) const {
                std::string mustr1 = "triggerObjectMatchesByPath(\"HLT_Mu13_Mu8_v*\").size()";
                std::string elstr1 = "triggerObjectMatchesByPath(\"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*\").size()";
                std::string mustr2 = "triggerObjectMatchesByPath(\"HLT_Mu17_Mu8_v*\").size()";
                std::string elstr2 = "triggerObjectMatchesByPath(\"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*\").size()";
                return lval(iz, il, mustr1, elstr1) > 0 || lval(iz, il, mustr2, elstr2) > 0;
            }

            float lisoEcalBaseline(unsigned int iz, unsigned int il) const {
                return HZZIso::getEcalIso(luserFloat(iz,il,"ecalZZ4L"), lval(iz, il, "userFloat('rhoMu')", "userFloat('rhoEl')"), lpt(iz, il), abs(leta(iz, il)), abs(lpdgId(iz, il)) == 13);
            }

            float lisoHcalBaseline(unsigned int iz, unsigned int il) const {
                return HZZIso::getHcalIso(luserFloat(iz,il,"hcalZZ4L"), lval(iz, il, "userFloat('rhoMu')", "userFloat('rhoEl')"), lpt(iz, il), abs(leta(iz, il)), abs(lpdgId(iz, il)) == 13);
            }
            float lisoCombRelBaseline(unsigned int iz, unsigned int il) const {
                return (lisoTrkBaseline(iz,il) + lisoEcalBaseline(iz,il) + lisoHcalBaseline(iz,il))/lpt(iz,il);
            }

            float lisoDirectional(unsigned int iz, unsigned int il, float coneDR=0.4, bool falloff=false) const;

            // return the wose sum of comb rel iso, for baseline cut
            float worsePairCombRelIsoBaseline() const ;

            float lisoPf(unsigned int iz, unsigned int il, const char *name)  const ; 
            float lisoPf(unsigned int iz, unsigned int il, const std::string &name)  const {
                return lisoPf(il,iz,name.c_str()); 
            }

            float luserFloat(unsigned int iz, unsigned int il, const char *label)  const ;
            float luserFloat(unsigned int iz, unsigned int il, const std::string &label)  const ;
            int   luserInt(unsigned int iz, unsigned int il, const char *label)  const ;
            int   luserInt(unsigned int iz, unsigned int il, const std::string &label)  const ;
            bool  lgood(unsigned int iz, unsigned int il, const char *muId, const char *eleId) const ;
            bool  lgood(unsigned int iz, unsigned int il, const std::string &muId, const std::string &eleId) const ;
            float lval(unsigned int iz, unsigned int il, const char *expr) const ;
            float lval(unsigned int iz, unsigned int il, const std::string &expr) const ;
            float lval(unsigned int iz, unsigned int il, const char *muExpr, const char *eleExpr) const ;
            float lval(unsigned int iz, unsigned int il, const std::string &muExpr, const std::string &eleExpr) const ;

            reco::GenParticleRef genh() const {
                reco::GenParticleRef z1 = genz(0), z2 = genz(1);
                if (z1.isNonnull() && z2.isNonnull() && z1 != z2) {
                    return z1->motherRef();
                } else {
                    return reco::GenParticleRef();
                }
            }
            reco::Particle::LorentzVector gen4lp4() const ;
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

            /// cut on the number of pairs that satisfy a cut.
            /// each pair is a reco::SkimEvent4L::Pair with daughters being the two candidates (no shallow clones)
            /// if anySign = true,  there are 6 pairs: (1,2), (1,3), (1,4), (2,3), (3,4), (3,4)
            /// if anySign = false, there are 4 OS pairs: (1,2), (3,4), (1,4), (2,3)
            const int nGoodPairs(const char *pairCut, bool anySign) const ;
            const int nGoodPairs(const std::string &pairCut, int anySign) const ;

            const int nMatchedLeptons() const { return nMatchedLeptonsInZ(0) + nMatchedLeptonsInZ(1); }
            const int nMatchedLeptonsInZ(int iz) const { return genl(iz,0).isNonnull() + genl(iz,1).isNonnull(); }

            bool hasCorrectAssignment() const {
                reco::GenParticleRef z1 = genz(0), z2 = genz(1);
                return (z1.isNonnull() && z2.isNonnull() && z1 != z2);
            }

            unsigned int intimeSimVertices() const {return numsimvertices_;}
            float numTrueInteractions() const {return numTrueInteractions_;}
            unsigned int numRecoVertices() const {return numrecovertices_;}

            float getCosThetaStar() const {return costhetastar_;}
            float getCosTheta1() const {return costheta1_;}
            float getCosTheta2() const {return costheta2_;}
            float getPhi() const {return phi_;}
            float getPhi1() const {return phi1_;}
            float getPhi2() const {return phi2_;}
            float getPhiStar1() const {return phistar1_;}
            float getPhiStar2() const {return phistar2_;}

            float getRho() const {
                if (abs(lpdgId(0, 0)) == 13) return luserFloat(0, 0,"rhoMu");
                else return luserFloat(0, 0,"rhoEl");
            }    
            float getNumRecoVertices() const {return numrecovertices_;}

            using reco::LeafCandidate::setVertex;
            void setVertex(const edm::Handle<reco::VertexCollection> &);
            void setPFMet(const edm::Handle<reco::PFMETCollection> &);
            void setJets(const edm::Handle<pat::JetCollection> &, double ptMin=-1);
            void setPFLeaves(const edm::Handle<std::vector<reco::LeafCandidate> >&);
            void setNumRecoVertices(const edm::Handle<reco::VertexCollection> & vtxH) {numrecovertices_ = vtxH->size();}

            void setAngles();

            void setGenMatches(const edm::Association<reco::GenParticleCollection> &genMatch) ;

            void setPileupInfo(std::vector<PileupSummaryInfo>);
        protected:
            /// return the proxy of a lepton (ShallowCloneCandidate or ShallowClonePtrCandidate)
            const reco::Candidate * lproxy(unsigned int iz, unsigned int il) const {
                return daughter(iz)->masterClone()->daughter(il);
            }

        private:
            void init() ;

            hypoType hypo_;

            unsigned int numsimvertices_, numrecovertices_;
	    float numTrueInteractions_;

            reco::VertexRef  vtx_;
            reco::PFMETRef   pfMet_;
            pat::JetRefVector jets_;
            std::vector<reco::LeafCandidate> pfLeaves_;

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

            float costheta1_; 
            float costheta2_; 
            float phi_; 
            float costhetastar_; 
            float phistar1_; 
            float phistar2_; 
            float phi1_; 
            float phi2_;
    };

}

#endif

