#ifndef AnalysisDataFormats_SkimEvent_h
#define AnalysisDataFormats_SkimEvent_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <vector>
#include <utility>

using std::size_t;

const double Mw = 80.450;
const double Mz = 91.188;

namespace reco {

    struct mwlSortByPtClass {
        bool operator() ( edm::RefToBase<Candidate> a, edm::RefToBase<Candidate> b) { return a->pt() > b->pt(); }
    };

    class SkimEvent : public LeafCandidate {
        //class SkimEvent {

        public:
            enum hypoType {undefined = 0, WELNU = 1, WMUNU=2, WWELEL=3, WWELMU=4, WWMUMU=5,hypoTypeSize=6};

            static const std::string hypoTypeNames[];

            /// static functions used to convert to&from hypoType and strings
            static std::string hypoTypeName(SkimEvent::hypoType);
            static hypoType hypoTypeByName(const std::string &name);


            SkimEvent();
            SkimEvent(const hypoType &);

            //Lepton variables
            const int hypo() const { return hypo_; }
            const int nLep(float a = -1) const;
            const int nExtraLep(float a = -1) const;
            const int nSoftMu(float a = -1) const;
            //const int nMu(float a = -1) const;
            //const int nExtraMu(float a = -1) const;
            //const int nExtraEl(float a = -1) const;
            //const int nEl(float a = -1) const;
            const int pdgId(size_t a = 0) const;
            /*
               const float sigmaIEtaIEta(size_t a = 0) const;
               const float deltaEtaIn(size_t a = 0) const;
               const float deltaPhiIn(size_t a = 0) const;
               const float hOverE(size_t a = 0) const;
               const int classification(size_t a = 0) const;
               const float isoTracker(size_t a = 0) const;
               const float isoEcal(size_t a = 0) const;
               const float isoHcal(size_t a = 0) const;
               const float isoSum(size_t a = 0) const;
             */
            const float pt(size_t a = 0) const;
            const float ptMax() const {return std::max(pt(0),pt(1));}
            const float ptMin() const {return std::min(pt(0),pt(1));}
            const float eta(size_t a = 0) const;
            const float phi(size_t a = 0) const;
            //const int missingHits(size_t a = 0) const;
            const float dcot(size_t a = 0) const;
            const float dist(size_t a = 0) const;
            /*
               const bool isEE(size_t a = 0) const;
               const bool isEB(size_t a = 0) const;
               const float eSeedOverPout(size_t a = 0) const;
               const bool isTracker(size_t a = 0) const;
               const bool isGlobal(size_t a = 0) const;
               const int nTrackerHits(size_t a = 0) const;
               const int nPixelHits(size_t a = 0) const;
               const int nChambers(size_t a = 0) const;
               const int nMatches(size_t a = 0) const;
               const int ndof(size_t a = 0) const;
               const float chi2ndof(size_t a = 0) const;
               const int muonHits(size_t a = 0) const;
               const float dXY(size_t a = 0) const;
               const float dZ(size_t a = 0) const;
               const float dZ(size_t i, const math::XYZPoint &p) const;
               const bool isTMOneStation(size_t a = 0) const;
               const bool isTM2DLoose(size_t a = 0) const;
             */
            const int q(size_t a = 0) const;

            //Jet variables
            const int nJets(float a = 30, bool = false) const;
            const int nCentralJets(float pt = 20,float eta=3.0,bool applyCorrection=false) const;
            const float jetPt(size_t a = 0,bool = false) const;
            static void setupJEC(const std::string&, const std::string&, const std::string&);

            //Event variables
            const bool leptEtaCut(float maxAbsEtaMu=2.4,float maxAbsEtaEl=2.5) const;
	    const bool triggerMatchingCut() const;
            const bool eleExpHitCut(bool isNew=1) const;
            const float pfMet() const;
            const float tcMet() const;
            const float mll() const;
            const float pTll() const;
            const float dPhill() const;
            const float dRll() const;
            const float dEtall() const;
            const float etall() const;
            const float yll() const;
            const float dPhillPfMet() const;
            const float dPhillTcMet() const;
            //const float dPhilPfMet(size_t a = 0) const;
            //const float dPhilTcMet(size_t a = 0) const;
            const float projPfMet() const;
            const float projTcMet() const;
            //const float pfMT(size_t a = 0) const;
            //const float tcMT(size_t a = 0) const;
            //const float pfMTll(size_t a = 0, size_t b = 1) const;
            //const float tcMTll(size_t a = 0, size_t b = 1) const;
            //const float llType(size_t a = 0, size_t b = 1) const;
            const float nTracks() const;
            //const float dPhiJetTcMet(size_t a = 0, size_t b = 1) const;
            //const float dPhiJetPfMet(size_t a = 0, size_t b = 1) const;
            const float cosThetaStar(size_t a = 0) const;
            const float afb(size_t a = 0) const;
            //const size_t dPhiPfMetMin(size_t a=0, size_t b=0) const;
            //const size_t dPhiTcMetMin(size_t a=0, size_t b=0) const;

            /*
               void setMuons(const edm::Handle<pat::MuonCollection> &);
            ///To be called only after setMuons()			       
            void setExtraMuons(const edm::Handle<pat::MuonCollection> & extra); 
            void setAllMuons(const edm::Handle<pat::MuonCollection> & base,
            const edm::Handle<pat::MuonCollection> & extra);

            void setElectrons(const edm::Handle<pat::ElectronCollection> &);
            ///To be called only after setElectronss()			       
            void setExtraElectrons(const edm::Handle<pat::ElectronCollection> & extra);				   
            void setAllElectrons(const edm::Handle<pat::ElectronCollection> & base,
            const edm::Handle<pat::ElectronCollection> & extra);
             */

            //Selection Functions
            const bool hasGoodVertex() const;
            const double d0SnT(size_t a=0) const;
            const double dZSnT(size_t a=0) const;
            const double d0Reco(size_t a=0) const;
            const double d0RecoSPT2(size_t a=0) const;
            const double dZReco(size_t a=0) const;
            const bool passesIDV1(size_t a=0) const;
            const bool passesConversion(size_t a=0) const;
            const bool isSTA(size_t a=0) const;
            const bool isMuTriggered(size_t a=0) const;

            void setLepton(const pat::Electron& ele);
            void setLepton(const pat::Muon& mu);

            void setSoftMuon(const pat::Muon& mu);

            void setExtraLepton(const pat::Electron& ele);
            void setExtraLepton(const pat::Muon& mu);


            void setJets(const edm::Handle<pat::JetCollection> &);
            void setTagJets(const edm::Handle<pat::JetCollection> &);
            void setTCMet(const edm::Handle<reco::METCollection> &);
            void setPFMet(const edm::Handle<reco::PFMETCollection> &);
            void setVertex(const edm::Handle<reco::VertexCollection> &);
            void setVtxSumPts(const edm::Handle<edm::ValueMap<float> > &s);
            void setVtxSumPt2s(const edm::Handle<edm::ValueMap<float> > &s);

            //void sortLepsByPt()     { std::sort(leps_.begin(),    leps_.end(),    mwlSortByPt); }
            //void sortMusByPt()      { std::sort(mus_.begin(),     mus_.end(),     mwlSortByPt); }
            //void sortElsByPt()      { std::sort(els_.begin(),     els_.end(),     mwlSortByPt); }
            //void sortExtraMusByPt() { std::sort(extraMus_.begin(),extraMus_.end(),mwlSortByPt); }
            //void sortExtraElsByPt() { std::sort(extraEls_.begin(),extraEls_.end(),mwlSortByPt); }
            void sortJetsByPt()     { std::sort(jets_.begin(),    jets_.end(),    mwlSortByPt); }
            void sortTagJetsByPt()     { std::sort(tagJets_.begin(),    tagJets_.end(),    mwlSortByPt); }


            const bool passesVtxSel(size_t a=0) const;
            const reco::Vertex highestPtVtx() const;
            const int bTaggedJetsUnder(const float&, const float&) const;


            //Iso Functions
            const bool isEB(size_t a = 0) const;
            const float tkPt(size_t a = 0) const;

            const size_t indexByPt(size_t a = 0) const;
            const size_t indexByIso(size_t a = 0) const;

            const float tkIso(size_t a = 0) const;
            const float ecalIso(size_t a = 0) const;
            const float hcalIso(size_t a = 0) const;
            const float allIso(size_t a = 0) const;

            const float tkVeto(size_t a = 0) const;
            const float ecalVeto(size_t a = 0) const;
            const float hcalVeto(size_t a = 0) const;
            const float allVeto(size_t a = 0) const;

            const float tkIsoByPt   (size_t i = 0) const { return tkIso  (indexByPt (i)); }
            const float tkIsoByIso  (size_t i = 0) const { return tkIso  (indexByIso(i)); }
            const float ecalIsoByPt (size_t i = 0) const { return ecalIso(indexByPt (i)); }
            const float ecalIsoByIso(size_t i = 0) const { return ecalIso(indexByIso(i)); }
            const float hcalIsoByPt (size_t i = 0) const { return hcalIso(indexByPt (i)); }
            const float hcalIsoByIso(size_t i = 0) const { return hcalIso(indexByIso(i)); }
            const float allIsoByPt  (size_t i = 0) const { return allIso (indexByPt (i)); }
            const float allIsoByIso (size_t i = 0) const { return allIso (indexByIso(i)); }
            const float tkVetoByPt   (size_t i = 0) const { return tkVeto  (indexByPt (i)); }
            const float tkVetoByIso  (size_t i = 0) const { return tkVeto  (indexByIso(i)); }
            const float ecalVetoByPt (size_t i = 0) const { return ecalVeto(indexByPt (i)); }
            const float ecalVetoByIso(size_t i = 0) const { return ecalVeto(indexByIso(i)); }
            const float hcalVetoByPt (size_t i = 0) const { return hcalVeto(indexByPt (i)); }
            const float hcalVetoByIso(size_t i = 0) const { return hcalVeto(indexByIso(i)); }
            const float allVetoByPt  (size_t i = 0) const { return allVeto (indexByPt (i)); }
            const float allVetoByIso (size_t i = 0) const { return allVeto (indexByIso(i)); }
            const int   pdgIdByPt   (size_t i = 0) const { return pdgId  (indexByPt (i)); }
            const int   pdgIdByIso  (size_t i = 0) const { return pdgId  (indexByIso(i)); }
            const float ptByPt      (size_t i = 0) const { return pt     (indexByPt (i)); }
            const float ptByIso     (size_t i = 0) const { return pt     (indexByIso(i)); }
            const float etaByPt     (size_t i = 0) const { return eta    (indexByPt (i)); }
            const float etaByIso    (size_t i = 0) const { return eta    (indexByIso(i)); }
            const float phiByPt     (size_t i = 0) const { return phi    (indexByPt (i)); }
            const float phiByIso    (size_t i = 0) const { return phi    (indexByIso(i)); }
            const int   qByPt       (size_t i = 0) const { return q      (indexByPt (i)); }
            const int   qByIso      (size_t i = 0) const { return q      (indexByIso(i)); }
            const bool  isEBByPt    (size_t i = 0) const { return isEB   (indexByPt (i)); }
            const bool  isEBByIso   (size_t i = 0) const { return isEB   (indexByIso(i)); }
            const float tkPtByPt    (size_t i = 0) const { return tkPt   (indexByPt (i)); }
            const float tkPtByIso   (size_t i = 0) const { return tkPt   (indexByIso(i)); }


        private:
            static mwlSortByPtClass mwlSortByPt;
            static std::vector<std::string> jecFiles_;

            int hypo_;
//             math::XYZPoint vtxPoint_;
//             reco::VertexRef vtxRef_;
//             reco::VertexCollection vtxs_;
//             std::vector<reco::VertexRef> vtxs_;
            reco::VertexRefVector vtxs_;
            std::vector<double> sumPts_;
            std::vector<double> sumPt2s_;
            reco::METRef tcMet_;
            reco::PFMETRef pfMet_;
            //std::vector<edm::RefToBase<Candidate> > leps_;
            //std::vector<edm::RefToBase<Candidate> > mus_;
            //std::vector<edm::RefToBase<Candidate> > els_;
            //std::vector<edm::RefToBase<Candidate> > softMus_;
            //std::vector<edm::RefToBase<Candidate> > extraLeps_;
            edm::OwnVector<reco::RecoCandidate> leps_;
            edm::OwnVector<reco::RecoCandidate> extraLeps_;
            std::vector<pat::Muon> softMuons_;
            std::vector<edm::RefToBase<Candidate> > jets_;
            std::vector<edm::RefToBase<Candidate> > tagJets_;

            //JEC


    };

    }

#endif

