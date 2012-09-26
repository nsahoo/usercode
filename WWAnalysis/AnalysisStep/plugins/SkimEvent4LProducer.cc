#include <vector>
#include <algorithm>
#include <cmath>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"

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

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Config.h"
#include "TMVA/MethodBDT.h"

class SkimEvent4LProducer : public edm::EDProducer {
    public:
        SkimEvent4LProducer(const edm::ParameterSet &iConfig) ;
        ~SkimEvent4LProducer() ;
        virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;
    private:
        edm::InputTag src_;

        edm::InputTag reducedPFCands_;
        edm::InputTag jets_;
        edm::InputTag pfMet_;
        edm::InputTag vertices_;
        edm::InputTag rho_;

        bool          isMC_;
        edm::InputTag gensTag_;
        std::string higgsmassweightfile_;
        HiggsMassWeightProvider* hmwp;
        bool          isSignal_;
        bool doswap;
        edm::InputTag mcMatch_;
        bool doMELA_;
        bool doAnglesWithFSR_;
        bool doMassRes_;
        bool doBDT_;

        std::auto_ptr<HZZ4LMelaDiscriminator> melaSMH_, melaPSH_, melaQQZZ_;
        CompositeCandMassResolution massRes_;

        //BDT input variables
        TMVA::Reader *ScalarVsBkgBDTReader;
        std::string weightfileScalarVsBkg_;
        float MVAInputVar_costheta1;
        float MVAInputVar_costheta2;
        float MVAInputVar_costhetastar;
        float MVAInputVar_Phi;
        float MVAInputVar_Phi1;
        float MVAInputVar_mZ1;
        float MVAInputVar_mZ2;
        float MVAInputVar_ZZpt;
        float MVAInputVar_zzdotz1;
        float MVAInputVar_zzdotz2;
        float MVAInputVar_ZZptZ1ptCosDphi;
        float MVAInputVar_ZZptZ2ptCosDphi;
        float MVAInputVar_reduced_Z1pt;
        float MVAInputVar_reduced_Z2pt;
        float MVAInputVar_ZZy;
        float MVAInputVar_m4l;
};

SkimEvent4LProducer::SkimEvent4LProducer(const edm::ParameterSet &iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    reducedPFCands_(iConfig.getParameter<edm::InputTag>("reducedPFCands")),
    jets_(iConfig.getParameter<edm::InputTag>("jets")),
    pfMet_(iConfig.getParameter<edm::InputTag>("pfMet")),
    vertices_(iConfig.getParameter<edm::InputTag>("vertices")),
    isMC_(iConfig.getParameter<bool>("isMC")),
    gensTag_(isMC_ ? (iConfig.existsAs<bool>("gensTag") ? iConfig.getParameter<edm::InputTag>("gensTag") : edm::InputTag("prunedGen")) : edm::InputTag("NOGENHERE")),
    higgsmassweightfile_(iConfig.existsAs<std::string>("higgsMassWeightFile")?iConfig.getParameter<std::string>("higgsMassWeightFile"):""),
    isSignal_(iConfig.existsAs<bool>("isSignal")?iConfig.getParameter<bool>("isSignal"):false),
    doswap(iConfig.existsAs<bool>("doswap")?iConfig.getParameter<bool>("doswap"):true),
    mcMatch_(isSignal_ ? iConfig.getParameter<edm::InputTag>("mcMatch") : edm::InputTag("FAKE")),
    doMELA_(iConfig.existsAs<bool>("doMELA")?iConfig.getParameter<bool>("doMELA"):false),
    doAnglesWithFSR_(iConfig.existsAs<bool>("doAnglesWithFSR")?iConfig.getParameter<bool>("doAnglesWithFSR"):true),
    doMassRes_(iConfig.existsAs<bool>("doMassRes")?iConfig.getParameter<bool>("doMassRes"):false),
    doBDT_(iConfig.existsAs<bool>("doBDT")?iConfig.getParameter<bool>("doBDT"):false),
    weightfileScalarVsBkg_(iConfig.existsAs<std::string>("weightfile_ScalarVsBkgBDT")?iConfig.getParameter<std::string>("weightfile_ScalarVsBkgBDT"):"")
{
    if (doMELA_) {
        std::string spath = edm::FileInPath(iConfig.getParameter<std::string>("melaQQZZHistos")).fullPath();
        const char *cpath = spath.c_str();
        melaSMH_.reset(new HZZ4LMelaDiscriminator(HZZ4LMelaDiscriminator::SMHiggs, NULL));
        melaPSH_.reset(new HZZ4LMelaDiscriminator(HZZ4LMelaDiscriminator::PSHiggs, NULL));
        melaQQZZ_.reset(new HZZ4LMelaDiscriminator(HZZ4LMelaDiscriminator::QQZZ,   cpath));
    }
    if (doBDT_) {
      ScalarVsBkgBDTReader = new TMVA::Reader( "V" );

      ScalarVsBkgBDTReader->AddVariable( "costheta1",		&MVAInputVar_costheta1);
      ScalarVsBkgBDTReader->AddVariable( "costheta2",		&MVAInputVar_costheta2);
      ScalarVsBkgBDTReader->AddVariable( "costhetastar",      &MVAInputVar_costhetastar);
      ScalarVsBkgBDTReader->AddVariable( "Phi",		&MVAInputVar_Phi);
      ScalarVsBkgBDTReader->AddVariable( "Phi1",		&MVAInputVar_Phi1);
      ScalarVsBkgBDTReader->AddVariable( "mZ1",		&MVAInputVar_mZ1);
      ScalarVsBkgBDTReader->AddVariable( "mZ2",		&MVAInputVar_mZ2);
      ScalarVsBkgBDTReader->AddVariable( "ZZpt/m4l",		&MVAInputVar_ZZpt);
      ScalarVsBkgBDTReader->AddVariable( "ZZdotZ1/(m4l*mZ1)",	&MVAInputVar_zzdotz1);
      ScalarVsBkgBDTReader->AddVariable( "ZZdotZ2/(m4l*mZ2)",	&MVAInputVar_zzdotz2);
      ScalarVsBkgBDTReader->AddVariable( "ZZptCosDphiZ1pt",	&MVAInputVar_ZZptZ1ptCosDphi);
      ScalarVsBkgBDTReader->AddVariable( "ZZptCosDphiZ2pt",	&MVAInputVar_ZZptZ2ptCosDphi);
      ScalarVsBkgBDTReader->AddVariable( "Z1pt/m4l",		&MVAInputVar_reduced_Z1pt);
      ScalarVsBkgBDTReader->AddVariable( "Z2pt/m4l",		&MVAInputVar_reduced_Z2pt);
      ScalarVsBkgBDTReader->AddVariable( "ZZy",		&MVAInputVar_ZZy);
      // add spectators
      ScalarVsBkgBDTReader->AddSpectator("m4l",               &MVAInputVar_m4l);
      // initialize
      ScalarVsBkgBDTReader->BookMVA("BDTG", edm::FileInPath(weightfileScalarVsBkg_.c_str()).fullPath());
      
    }

    const char *higgsweightbase=getenv("CMSSW_BASE");
    std::string higgsweightpath(higgsweightbase);
    higgsweightpath += "/src/";
    if (higgsmassweightfile_ != "") higgsweightpath += higgsmassweightfile_;
    else higgsweightpath = "";

    hmwp = new HiggsMassWeightProvider(higgsweightpath);

    produces<std::vector<reco::SkimEvent4L> >();
}

SkimEvent4LProducer::~SkimEvent4LProducer() {
    delete hmwp;
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

    edm::Handle<reco::GenParticleCollection> gensH;
    if (isMC_) iEvent.getByLabel(gensTag_, gensH);

    float genhiggsmass = 0.0;
    if (isMC_) {
        reco::GenParticleCollection gens = *gensH;
        for (std::size_t i = 0; i < gens.size(); i++) {
            if (gens[i].pdgId() == 25 && gens[i].status() == 3) {
                genhiggsmass = gens[i].mass();
            }
        }
    }

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
        zz.setGenHiggsMass(genhiggsmass);
        zz.setHiggsMassWeight(hmwp);

        if (doMELA_) {
            zz.addUserFloat("melaSMH",  melaSMH_->get( zz.mass(), zz.mz(0), zz.mz(1), zz.getCosThetaStar(), zz.getCosTheta1(), zz.getCosTheta2(), zz.getPhi(), zz.getPhi1()));
            zz.addUserFloat("melaPSH",  melaPSH_->get( zz.mass(), zz.mz(0), zz.mz(1), zz.getCosThetaStar(), zz.getCosTheta1(), zz.getCosTheta2(), zz.getPhi(), zz.getPhi1()));
            zz.addUserFloat("melaQQZZ", melaQQZZ_->get(zz.mass(), zz.mz(0), zz.mz(1), zz.getCosThetaStar(), zz.getCosTheta1(), zz.getCosTheta2(), zz.getPhi(), zz.getPhi1()));
        }

        if (doBDT_) {
          
          MVAInputVar_costheta1 = zz.getCosTheta1();
          MVAInputVar_costheta2 = zz.getCosTheta2();
          MVAInputVar_costhetastar =  zz.getCosThetaStar();
          MVAInputVar_Phi = zz.getPhi();
          MVAInputVar_Phi1 = zz.getPhi1();
          MVAInputVar_mZ1 = zz.mz(0);
          MVAInputVar_mZ2 = zz.mz(1);
          MVAInputVar_ZZpt = zz.pt();
          MVAInputVar_zzdotz1 = zz.p4().Dot(zz.z(0).p4());
          MVAInputVar_zzdotz2 = zz.p4().Dot(zz.z(1).p4());


          TLorentzVector Z1vec; Z1vec.SetPtEtaPhiM( zz.z(0).p4().Pt(), zz.z(0).p4().Eta(), zz.z(0).p4().Phi(), zz.z(0).p4().M() );
          TLorentzVector Z2vec; Z2vec.SetPtEtaPhiM( zz.z(1).p4().Pt(), zz.z(1).p4().Eta(), zz.z(1).p4().Phi(), zz.z(1).p4().M() );
          TLorentzVector ZZvec = Z1vec + Z2vec;
          TVector3 zzvec3 = Z1vec.Vect();
          TVector3 z1vec3 = Z2vec.Vect();
          TVector3 z2vec3 = ZZvec.Vect();
          TVector2 zzvecxy = zzvec3.XYvector();
          TVector2 z1vecxy = z1vec3.XYvector();
          TVector2 z2vecxy = z2vec3.XYvector();
          MVAInputVar_ZZptZ1ptCosDphi = TMath::Cos(zzvecxy.DeltaPhi( z1vecxy));
          MVAInputVar_ZZptZ2ptCosDphi = TMath::Cos(zzvecxy.DeltaPhi( z2vecxy));
          MVAInputVar_reduced_Z1pt = zz.z(0).p4().Pt() / zz.mass();
          MVAInputVar_reduced_Z2pt = zz.z(1).p4().Pt() / zz.mass();
          MVAInputVar_ZZy = zz.p4().Rapidity();
          MVAInputVar_m4l = zz.mass();
          zz.addUserFloat("BDT_ScalarVsBkg_125", ScalarVsBkgBDTReader->EvaluateMVA("BDTG"));
        }

        if (doMassRes_) zz.addUserFloat("massErr", massRes_.getMassResolution(zz));
        if (isMC_) zz.setPileupInfo(*puH);
        if (isSignal_) zz.setGenMatches(*mcMatch);
    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SkimEvent4LProducer);
