// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"

#include <vector>
#include <string>

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaPhi.h"


//
// class declaration
//

class PatElectronEffAreaIso : public edm::EDProducer {
    public:
        explicit PatElectronEffAreaIso(const edm::ParameterSet&);
        ~PatElectronEffAreaIso() {}
    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);
        // ----------member data ---------------------------
        edm::InputTag electronTag_;
        ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_;
        ElectronEffectiveArea::ElectronEffectiveAreaType   effAreaGamma_, effAreaNeutralHad_;
        std::string rho_; 
        std::string deltaR_; 
        std::string label_; 
};

PatElectronEffAreaIso::PatElectronEffAreaIso(const edm::ParameterSet& iConfig) :
    electronTag_(iConfig.getParameter<edm::InputTag>("src")),
    rho_(iConfig.getParameter<std::string>("rho")),
    deltaR_(iConfig.getParameter<std::string>("deltaR")),
    label_(iConfig.getParameter<std::string>("label"))
{
  produces<pat::ElectronCollection>();  
  std::string eaTarget = iConfig.getParameter<std::string>("effectiveAreaTarget");
  if      (eaTarget == "NoCorr") effAreaTarget_ = ElectronEffectiveArea::kEleEANoCorr;
  else if (eaTarget == "Data2011") effAreaTarget_ = ElectronEffectiveArea::kEleEAData2011;
  else if (eaTarget == "Summer11MC") effAreaTarget_ = ElectronEffectiveArea::kEleEASummer11MC;
  else if (eaTarget == "Fall11MC") effAreaTarget_ = ElectronEffectiveArea::kEleEAFall11MC;
  else throw cms::Exception("Configuration") << "Unknown effective area " << eaTarget << "\n";
  if (deltaR_ == "03") {
      effAreaGamma_   = ElectronEffectiveArea::kEleGammaIso03;
      effAreaNeutralHad_ = ElectronEffectiveArea::kEleNeutralHadronIso03;
  } else if (deltaR_ == "04") {
      effAreaGamma_   = ElectronEffectiveArea::kEleGammaIso04;
      effAreaNeutralHad_ = ElectronEffectiveArea::kEleNeutralHadronIso04;
  } else throw cms::Exception("Configuration") << "Unsupported deltaR " << deltaR_ << "\n";
}


void PatElectronEffAreaIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<pat::Electron> > electrons;
    iEvent.getByLabel(electronTag_,electrons);

    std::auto_ptr<pat::ElectronCollection> pOut(new pat::ElectronCollection);

    // ----- here is the real loop over the electrons ----
    for(edm::View<pat::Electron>::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){    
      pat::Electron clone = *ele;
      float rho = clone.userFloat(rho_), abseta = abs(clone.eta());
      float eff_area_ga = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
      float eff_area_nh = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
      clone.addUserFloat(label_+"EApho",    eff_area_ga);
      clone.addUserFloat(label_+"EAneuHad", eff_area_nh);
      float iso = 0;
      iso += clone.userFloat("electronPFIsoChHad"+deltaR_);
      iso += max<float>(0., clone.userFloat("electronPFIsoNHad"+deltaR_)   - eff_area_nh*rho);
      iso += max<float>(0., clone.userFloat("electronPFIsoPhoton"+deltaR_) - eff_area_ga*rho);
      clone.addUserFloat(label_, iso);
      pOut->push_back(clone);
    }
    iEvent.put(pOut);
}

class PatMuonEffAreaIso : public edm::EDProducer {
    public:
        explicit PatMuonEffAreaIso(const edm::ParameterSet&);
        ~PatMuonEffAreaIso() {}
    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);
        // ----------member data ---------------------------
        edm::InputTag muonTag_;
        MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_;
        MuonEffectiveArea::MuonEffectiveAreaType   effAreaGamma_, effAreaNeutralHad_, effAreaTotal_;
        std::string rho_; 
        std::string deltaR_; 
        std::string label_; 
};

PatMuonEffAreaIso::PatMuonEffAreaIso(const edm::ParameterSet& iConfig) :
    muonTag_(iConfig.getParameter<edm::InputTag>("src")),
    rho_(iConfig.getParameter<std::string>("rho")),
    deltaR_(iConfig.getParameter<std::string>("deltaR")),
    label_(iConfig.getParameter<std::string>("label"))
{
  produces<pat::MuonCollection>();  
  std::string eaTarget = iConfig.getParameter<std::string>("effectiveAreaTarget");
  if      (eaTarget == "NoCorr") effAreaTarget_ = MuonEffectiveArea::kMuEANoCorr;
  else if (eaTarget == "Data2011") effAreaTarget_ = MuonEffectiveArea::kMuEAData2011;
  else if (eaTarget == "Data2012") effAreaTarget_ = MuonEffectiveArea::kMuEAData2012;
  else if (eaTarget == "Summer11MC") effAreaTarget_ = MuonEffectiveArea::kMuEASummer11MC;
  else if (eaTarget == "Fall11MC") effAreaTarget_ = MuonEffectiveArea::kMuEAFall11MC;
  else throw cms::Exception("Configuration") << "Unknown effective area " << eaTarget << "\n";
  if (deltaR_ == "03") {
      effAreaGamma_   = MuonEffectiveArea::kMuGammaIso03;
      effAreaNeutralHad_ = MuonEffectiveArea::kMuNeutralHadronIso03;
      effAreaTotal_ = MuonEffectiveArea::kMuGammaAndNeutralHadronIso03;
  } else if (deltaR_ == "04") {
      effAreaGamma_   = MuonEffectiveArea::kMuGammaIso04;
      effAreaNeutralHad_ = MuonEffectiveArea::kMuNeutralHadronIso04;
      effAreaTotal_ = MuonEffectiveArea::kMuGammaAndNeutralHadronIso04;
  } else throw cms::Exception("Configuration") << "Unsupported deltaR " << deltaR_ << "\n";
}


void PatMuonEffAreaIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<pat::Muon> > muons;
    iEvent.getByLabel(muonTag_,muons);

    std::auto_ptr<pat::MuonCollection> pOut(new pat::MuonCollection);

    // ----- here is the real loop over the muons ----
    for(edm::View<pat::Muon>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){    
      pat::Muon clone = *mu;
      float rho = clone.userFloat(rho_), abseta = abs(clone.eta());
      float eff_area_ga  = MuonEffectiveArea::GetMuonEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
      float eff_area_nh  = MuonEffectiveArea::GetMuonEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
      float eff_area_tot = MuonEffectiveArea::GetMuonEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
      clone.addUserFloat(label_+"EApho",    eff_area_ga);
      clone.addUserFloat(label_+"EAneuHad", eff_area_nh);
      clone.addUserFloat(label_+"EAtot",    eff_area_tot);
      float iso = 0;
      iso += clone.userFloat("muonPFIsoChHad"+deltaR_);
      iso += max<float>(0., clone.userFloat("muonPFIsoNHad"+deltaR_) + clone.userFloat("muonPFIsoPhoton"+deltaR_) - eff_area_tot*rho);
      clone.addUserFloat(label_, iso);
      pOut->push_back(clone);
    }
    iEvent.put(pOut);
}


//define this as a plug-in
DEFINE_FWK_MODULE(PatElectronEffAreaIso);
DEFINE_FWK_MODULE(PatMuonEffAreaIso);
