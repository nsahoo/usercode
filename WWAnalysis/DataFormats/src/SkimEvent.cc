#include "WWAnalysis/DataFormats/interface/SkimEvent.h"
#include <iostream>
#include <algorithm>
#include "Math/VectorUtil.h"
#include <DataFormats/PatCandidates/interface/JetCorrFactors.h>

std::vector<std::string> reco::SkimEvent::jecFiles_;

struct indexValueStruct {
    indexValueStruct(const float &v, const size_t &i) : value(v), index(i) {}
    float value;
    size_t index;
};

bool lowToHigh(const indexValueStruct &a, const indexValueStruct &b) {
    return a.value < b.value;
}

bool highToLow(const indexValueStruct &a, const indexValueStruct &b) {
    return a.value > b.value;
}

// To be kept in synch with the enumerator definitions in SkimEvent.h file
const std::string reco::SkimEvent::hypoTypeNames[] = { "undefined", "WELNU", "WMUNU", "WWELEL", "WWELMU", "WWMUEL", "WWMUMU"};


std::string reco::SkimEvent::hypoTypeName(reco::SkimEvent::hypoType a){
  if(int(a) < int(hypoTypeSize) && int(a)>0) return hypoTypeNames[int(a)];
  return "undefined";
}



reco::SkimEvent::hypoType reco::SkimEvent::hypoTypeByName(const std::string &name){
  hypoType size = hypoTypeSize;
  int index = std::find(hypoTypeNames, hypoTypeNames+size, name)-hypoTypeNames;
  if(index == size) return undefined; // better this or throw() ?

  // cast
  return hypoType(index);
}



reco::SkimEvent::SkimEvent() : 
        hypo_(-1), sumPts_(0)/*, jec_(0), vtxPoint_(0,0,0) */{ }

reco::SkimEvent::SkimEvent(const reco::SkimEvent::hypoType &h) :
        hypo_(h), sumPts_(0)/*, jec_(0), vtxPoint_(0,0,0) */{ }



void reco::SkimEvent::setLepton(const pat::Electron& ele){
  leps_.push_back(ele);
}

void reco::SkimEvent::setLepton(const pat::Muon& mu){
  leps_.push_back(mu);
}

void reco::SkimEvent::setExtraLepton(const pat::Electron& ele){
  extraLeps_.push_back(ele);
}

void reco::SkimEvent::setExtraLepton(const pat::Muon& mu){
  extraLeps_.push_back(mu);
}

void reco::SkimEvent::setSoftMuon(const pat::Muon& mu){
  softMuons_.push_back(mu);
}


void reco::SkimEvent::setJets(const edm::Handle<pat::JetCollection> & jH) {

    jets_.clear();

    for(size_t i=0;i<jH->size();++i) 
      jets_.push_back(pat::JetRef(jH,i));

    //sortJetsByPt();

}

void reco::SkimEvent::setTagJets(const edm::Handle<pat::JetCollection> & jH) {

    tagJets_.clear();

    for(size_t i=0;i<jH->size();++i) 
      tagJets_.push_back(pat::JetRef(jH,i));

    //sortTagJetsByPt();

}

void reco::SkimEvent::setTCMet(const edm::Handle<reco::METCollection> & mH) {
    tcMet_ = reco::METRef(mH,0);
}

void reco::SkimEvent::setPFMet(const edm::Handle<reco::PFMETCollection> & mH) {
    pfMet_ = reco::PFMETRef(mH,0);
}

void reco::SkimEvent::setVtxSumPts(const edm::Handle<edm::ValueMap<float> > &s) {

    for(size_t i=0;i<vtxs_.size();++i) sumPts_.push_back( (*s)[vtxs_[i]] );

}

void reco::SkimEvent::setVtxSumPt2s(const edm::Handle<edm::ValueMap<float> > &s) {

    for(size_t i=0;i<vtxs_.size();++i) sumPt2s_.push_back( (*s)[vtxs_[i]] );

}

void reco::SkimEvent::setVertex(const edm::Handle<reco::VertexCollection> & vtxH) {

    for(size_t i=0;i<vtxH->size();++i) vtxs_.push_back(reco::VertexRef(vtxH,i));
}


//Lepton variables

const int reco::SkimEvent::nLep(float minPt) const { 
    int count = 0;
    if(minPt < 0) count = leps_.size(); 
    else for(size_t i=0;i<leps_.size();++i) if(leps_[i].pt() > minPt) count++;
    return count;
} 

const int reco::SkimEvent::nExtraLep(float minPt) const { 
    int count = 0;
    if(minPt < 0) count = extraLeps_.size(); 
    else for(size_t i=0;i<extraLeps_.size();++i) if(extraLeps_[i].pt() > minPt) count++;
    return count;
}


const int reco::SkimEvent::nSoftMu(float minPt) const { 
    int count = 0;
    if(minPt < 0) count = softMuons_.size(); 
    else for(size_t i=0;i<softMuons_.size();++i) if(softMuons_[i].pt() > minPt) count++;
    return count;
} 


const int reco::SkimEvent::pdgId(size_t i) const { 
    if(i < leps_.size()) return leps_[i].pdgId();
    else return -9999;
}


/*
const pat::Muon& reco::SkimEvent::mu(size_t a=0) const{
  return static_cast<const pat::Muon&>(leps_[a]);
}

const pat::Electron& reco::SkimEvent::el(size_t a=0) const{
  return static_cast<const pat::Electron&>(leps_[a]);
}
*/

const float reco::SkimEvent::pt(size_t i) const {
  if(i >= leps_.size()) return -9999.0;
  return leps_[i].pt();
}

const float reco::SkimEvent::eta(size_t i) const {
  if(i < leps_.size()) return leps_[i].eta();
  else return -9999.0;
}

const float reco::SkimEvent::phi(size_t i) const {
  if(i < leps_.size()) return leps_[i].phi();
  else return -9999.0;
}


const int reco::SkimEvent::q(size_t i) const{
  if(i < leps_.size()) return leps_[i].charge();
  else return -9999.0;
}


//Jet variables
void reco::SkimEvent::setupJEC(const std::string &l2File, const std::string &l3File, const std::string &residualFile) {

    jecFiles_.clear();
    jecFiles_.push_back(l2File);
    jecFiles_.push_back(l3File);
    if(residualFile.compare("")) jecFiles_.push_back(residualFile);

}


// void reco::SkimEvent::setupJEC(const JetCorrector *c) {
//     jec_ = 

const int reco::SkimEvent::nJets(float minPt,int applyCorrection,int applyID) const {
  return nCentralJets(minPt,99.9,applyCorrection,applyID);
}

const bool reco::SkimEvent::isThisJetALepton(pat::JetRef jet) const {
    bool thisJetIsLepton(false);
    for(size_t j=0; j<leps_.size();++j){
        double dR = fabs(ROOT::Math::VectorUtil::DeltaR(jet->p4(),leps_[j].p4()) );
        if(dR < 0.3){ 
            thisJetIsLepton = true;
            break;
        }
    }

    return thisJetIsLepton;
}

const int reco::SkimEvent::nCentralJets(float minPt,float eta,int applyCorrection,int applyID) const {

    int count = 0;
    for(size_t i=0;i<jets_.size();++i) {
        if( std::fabs(jets_[i]->eta()) >= eta) continue;
        if( jetPt(i,applyCorrection) <= minPt) continue;
        if(isThisJetALepton(jets_[i]))  continue;
	if(applyID && jets_[i]->chargedMultiplicity()==0) continue;
        count++;
    }
    return count;
}

const float reco::SkimEvent::jetPt(size_t i, int applyCorrection) const {
  if(applyCorrection) return jets_[i]->correctedJet("L3Absolute","none","patJetCorrFactorsFastJet").pt();
  else                return jets_[i]->correctedJet("Uncorrected","none","patJetCorrFactorsFastJet").pt();
}

const float reco::SkimEvent::tagJetPt(size_t i, int applyCorrection) const {
  if(applyCorrection) return tagJets_[i]->correctedJet("L3Absolute","none","patJetCorrFactorsFastJet").pt();
  else                return tagJets_[i]->correctedJet("Uncorrected","none","patJetCorrFactorsFastJet").pt();
}



//Event variables
const float reco::SkimEvent::pfMet() const {

    if(pfMet_.isNonnull()) return pfMet_->pt();
    else return -9999.0;
}

const float reco::SkimEvent::tcMet() const {

    if(tcMet_.isNonnull()) return tcMet_->pt();
    else return -9999.0;
}

const float reco::SkimEvent::tcMetX() const {

    if(tcMet_.isNonnull()) return tcMet_->px();
    else return -9999.0;
}

const float reco::SkimEvent::tcMetY() const {

    if(tcMet_.isNonnull()) return tcMet_->py();
    else return -9999.0;
}

const float reco::SkimEvent::mll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0].p4() + leps_[1].p4()).mass();
}

const float reco::SkimEvent::pTll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0].p4() + leps_[1].p4()).pt();
}

const float reco::SkimEvent::dPhill() const {
  if(leps_.size()!=2) return -9999.0;
  return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0].p4(),leps_[1].p4()) );   
}

const float reco::SkimEvent::dRll() const {
  if(leps_.size()!=2) return -9999.0;
  return ROOT::Math::VectorUtil::DeltaR(leps_[0].p4(),leps_[1].p4());
}

const float reco::SkimEvent::dEtall() const {
  if(leps_.size()!=2) return -9999.0;
  return fabs(leps_[0].eta() - leps_[1].eta());
}

const float reco::SkimEvent::etall() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0].p4() + leps_[1].p4()).eta();
}

const float reco::SkimEvent::yll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0].p4() + leps_[1].p4()).y();
}

const float reco::SkimEvent::dPhillPfMet() const {
  if(leps_.size()!=2 || pfMet_.isNonnull()) return -9999.0;
  return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0].p4()+leps_[1].p4(),pfMet_->p4()) );
}

const float reco::SkimEvent::dPhillTcMet() const {
  if(leps_.size()!=2 || tcMet_.isNonnull()) return -9999.0;
  return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[0].p4()+leps_[1].p4(),tcMet_->p4()) );
}

const float reco::SkimEvent::dPhilPfMet(size_t i) const {
  if(i < leps_.size() && pfMet_.isNonnull()) {
    return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[i].p4(),pfMet_->p4()) );
  } else 
    return -9999.0;
}

const float reco::SkimEvent::dPhilTcMet(size_t i) const {
  if(i < leps_.size() && tcMet_.isNonnull()) {
    return fabs(ROOT::Math::VectorUtil::DeltaPhi(leps_[i].p4(),tcMet_->p4()) );
  } else 
    return -9999.0;
}


const float reco::SkimEvent::mTHiggs() const {
    //version 2 from guillelmo's talk
    return sqrt( mll()*mll() + 
                 2*( sqrt(pTll()*pTll()+mTll()*mTll()) * tcMet()  -
                     (pXll()+tcMetX())*(pXll()+tcMetX()) - 
                     (pYll()+tcMetY())*(pYll()+tcMetY()) ) );
}

const float reco::SkimEvent::pXll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0].p4() + leps_[1].p4()).px();
}

const float reco::SkimEvent::pYll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0].p4() + leps_[1].p4()).py();
}

const float reco::SkimEvent::mTll() const {
  if(leps_.size()!=2) return -9999.0;
  return (leps_[0].p4() + leps_[1].p4()).mt();
}

const float reco::SkimEvent::mT(size_t i) const {
  if(i>=leps_.size()) return -9999.0;
  return sqrt(2*pt(i)*tcMet()*(1 - cos(dPhilTcMet(i))));
}

const float reco::SkimEvent::projPfMet() const {
    float dphi = dPhilMet();
    if(dphi < M_PI/2.) return pfMet_->pt()*sin(dphi);
    else               return pfMet_->pt();       
}

const float reco::SkimEvent::projTcMet() const {
    float dphi = dPhilMet();
    if(dphi < M_PI/2.) return tcMet_->pt()*sin(dphi);
    else               return tcMet_->pt();       
}

const float reco::SkimEvent::dPhilMet() const {
    float smallestDphi = 9999.;
    for(size_t l=0; l<leps_.size();++l){
        float dphi = fabs(ROOT::Math::VectorUtil::DeltaPhi(tcMet_->p4(),leps_[l].p4()) );
        if( dphi < smallestDphi) smallestDphi = dphi;
    }
    return smallestDphi;
}



const float reco::SkimEvent::nTracks() const {

    return 0;
}

/*
const float reco::SkimEvent::dPhiJetPfMet() const {

    return 0;
}

const float reco::SkimEvent::dPhiJetTcMet() const {

    return 0;
}
*/

const float reco::SkimEvent::cosThetaStar(size_t i) const {

    return 0;
}

const float reco::SkimEvent::afb(size_t i) const {

    return 0;
}


const bool reco::SkimEvent::leptEtaCut(float maxAbsEtaMu,float maxAbsEtaEl) const{
  bool check0(true);
  bool check1(true);

  if(abs(leps_[0].pdgId())==11 && fabs(leps_[0].eta())>=maxAbsEtaEl) check0=false;
  if(abs(leps_[0].pdgId())==13 && fabs(leps_[0].eta())>=maxAbsEtaMu) check0=false;

  if(abs(leps_[1].pdgId())==11 && fabs(leps_[1].eta())>=maxAbsEtaEl) check1=false;
  if(abs(leps_[1].pdgId())==13 && fabs(leps_[1].eta())>=maxAbsEtaMu) check1=false;

  return (check0 && check1);
}

bool reco::SkimEvent::passTriggerSingleMu(size_t i) const{ 
  bool result(false);

  if( fabs(leps_[i].pdgId()) != 13 ) return false;

  const pat::Muon& mu = static_cast<const pat::Muon&>(leps_[i]);
  if(mu.triggerObjectMatchesByPath("HLT_Mu24*").size()) result=true;

  return result;
}

bool reco::SkimEvent::passTriggerDoubleMu(size_t i) const{
  using namespace std;
  bool result(false);
  
  if( fabs(leps_[i].pdgId()) != 13 ) return false;
  
  const pat::Muon& mu = static_cast<const pat::Muon&>(leps_[i]);
  if(mu.triggerObjectMatchesByPath("HLT_DoubleMu7*").size()) result=true;
  
  return result;
}

bool reco::SkimEvent::passTriggerDoubleEl(size_t i) const{ 
  bool result(false);
  
  if( fabs(leps_[i].pdgId()) != 11 ) return false;
  
  const pat::Electron& el = static_cast<const pat::Electron&>(leps_[i]);
  if(el.triggerObjectMatchesByPath("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL*").size() ||
     el.triggerObjectMatchesByPath("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL*").size() ) result=true;
  
  return result;

}
bool reco::SkimEvent::passTriggerElMu(size_t i) const{ 
  bool result(false);
  using namespace std;
  if( fabs(leps_[i].pdgId()) == 13 ) {
    const pat::Muon& mu = static_cast<const pat::Muon&>(leps_[i]);
    if(mu.triggerObjectMatchesByPath("HLT_Mu8_Ele17_CaloIdL*",false).size() ||
       mu.triggerObjectMatchesByPath("HLT_Mu17_Ele8_CaloIdL*",false).size() ) result=true;
  }

  if( fabs(leps_[i].pdgId()) == 11 ) {
    const pat::Electron& el = static_cast<const pat::Electron&>(leps_[i]);
    if(el.triggerObjectMatchesByPath("HLT_Mu8_Ele17_CaloIdL*").size() ||
       el.triggerObjectMatchesByPath("HLT_Mu17_Ele8_CaloIdL*").size() ) result=true;
  }

  return result;

}


// TO BE FIXED: This guy should take the list of trigger from some 
// sort of configuration file. 
const bool reco::SkimEvent::triggerMatchingCut(SkimEvent::primaryDatasetType pdType) const{
  using namespace std;
  if(pdType==MC)
    return true;

  bool result(false);  
  if(hypo()==WWMUMU){
    if(pdType==DoubleMuon){ //configuration (1)
      result=(passTriggerDoubleMu(0) && passTriggerDoubleMu(1));}
    if(pdType==SingleMuon) {//configuration (2)
      result=(  (passTriggerSingleMu(0) || passTriggerSingleMu(1)) &&
		!(passTriggerDoubleMu(0) && passTriggerDoubleMu(1)) );
    }
  }

  if(hypo()==WWMUEL || hypo()==WWELMU){
    if(pdType==SingleMuon) //configuration (3)
      result=(  passTriggerSingleMu(0) || passTriggerSingleMu(1) );
    if(pdType==MuEG){       //configuration (4)
      result=( (passTriggerElMu(0) && passTriggerElMu(1)) && 
	       !(passTriggerSingleMu(0) || passTriggerSingleMu(1)) );    
    }
  }

  if(hypo()==WWELEL){
    if(pdType==DoubleElectron)//configuration (5)
      result= (passTriggerDoubleEl(0) && passTriggerDoubleEl(1));
  }

  return result;

  /*
  bool result(false);
  for(unsigned int i=0; i<=1; i++){
    if( fabs(leps_[i].pdgId()) == 13 ) {
      const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
      if( !mu.triggerObjectMatchesByPath("HLT_Mu9").empty() || 
	  !mu.triggerObjectMatchesByPath("HLT_Mu15_v1").empty() ) result=true;
    } else if(fabs(leps_[i].pdgId()) == 11 ){
      const pat::Electron &el = static_cast<const pat::Electron&>(leps_[i]);
      if ( !el.triggerObjectMatchesByPath("HLT_Ele10_LW_L1R").empty() || 
	   !el.triggerObjectMatchesByPath("HLT_Ele15_SW_L1R").empty() ||
	   !el.triggerObjectMatchesByPath("HLT_Ele15_SW_CaloEleId_L1R").empty() ||
	   !el.triggerObjectMatchesByPath("HLT_Ele17_SW_CaloEleId_L1R").empty() ||
	   !el.triggerObjectMatchesByPath("HLT_Ele17_SW_TightEleId_L1R").empty() ||
	   !el.triggerObjectMatchesByPath("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2").empty() ||
	   !el.triggerObjectMatchesByPath("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3").empty() ) result=true;
    }
  } 
  return result;
  */
}


//Iso Functions
const float reco::SkimEvent::tkIso(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( fabs(leps_[i].pdgId()) == 11 ) {
        return static_cast<const pat::Electron&>(leps_[i]).dr03TkSumPt();
    } else if ( fabs(leps_[i].pdgId()) == 13 ) {
        return static_cast<const pat::Muon&>(leps_[i]).isolationR03().sumPt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::ecalIso(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( fabs(leps_[i].pdgId()) == 11 ) {
        return static_cast<const pat::Electron&>(leps_[i]).dr03EcalRecHitSumEt();
    } else if ( fabs(leps_[i].pdgId()) == 13 ) {
        return static_cast<const pat::Muon&>(leps_[i]).isolationR03().emEt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::hcalIso(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( fabs(leps_[i].pdgId()) == 11 ) {
        return static_cast<const pat::Electron&>(leps_[i]).dr03HcalTowerSumEt();
    } else if ( fabs(leps_[i].pdgId()) == 13 ) {
        return static_cast<const pat::Muon&>(leps_[i]).isolationR03().hadEt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::getRho(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;
    
    if( fabs(leps_[i].pdgId()) == 11 ) {
        return static_cast<const pat::Electron&>(leps_[i]).userFloat("rhoEl");
    } else if ( fabs(leps_[i].pdgId()) == 13 ) {
        return static_cast<const pat::Muon&>(leps_[i]).userFloat("rhoMu");
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::allIso(size_t i) const {
    
    if( i >= leps_.size() ) return -9999.0;

    if( fabs(leps_[i].pdgId()) == 11 && isEB(i) ) {
        return tkIso(i) + std::max((float)0,ecalIso(i)-1) + hcalIso(i) - getRho(i) * 3.14159265 * 0.3 * 0.3;
    } else if( (fabs(leps_[i].pdgId()) == 11 && isEE(i)) || fabs(leps_[i].pdgId()) == 13 ) {
        return tkIso(i) + ecalIso(i) + hcalIso(i) - getRho(i) * 3.14159265 * 0.3 * 0.3;
    } else {
        std::cout << " Do I ever friggin get here?" << std::endl;
        return 9999.0;
    }
}

const float reco::SkimEvent::tkVeto(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( fabs(leps_[i].pdgId()) == 11 ) {
        return 0;
    } else if ( fabs(leps_[i].pdgId()) == 13 ) {
        return static_cast<const pat::Muon&>(leps_[i]).isolationR03().trackerVetoPt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::ecalVeto(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( fabs(leps_[i].pdgId()) == 11 ) {
        return 0;
    } else if ( fabs(leps_[i].pdgId()) == 13 ) {
        return static_cast<const pat::Muon&>(leps_[i]).isolationR03().emVetoEt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::hcalVeto(size_t i) const {

    if( i >= leps_.size() ) return -9999.0;

    if( fabs(leps_[i].pdgId()) == 11 ) {
        return 0;
    } else if ( fabs(leps_[i].pdgId()) == 13 ) {
        return static_cast<const pat::Muon&>(leps_[i]).isolationR03().hadVetoEt;
    } else {
        return -9999.0;
    }
}

const float reco::SkimEvent::allVeto(size_t i) const {
    
    if( i >= leps_.size() ) return -9999.0;

    if( fabs(leps_[i].pdgId()) == 11 ) {
        return 0;
    } else if( fabs(leps_[i].pdgId()) == 13 ) {
        return tkVeto(i) + ecalVeto(i) + hcalVeto(i);
    } else {
        return -9999.0;
    }
}

const size_t reco::SkimEvent::indexByPt(size_t i) const {

    if( i >= leps_.size() ) return 9999;
    std::vector<indexValueStruct> a;

    for(size_t j=0;j<leps_.size();++j) a.push_back(indexValueStruct(pt(j),j));
    std::sort(a.begin(),a.end(),highToLow);

    return a[i].index;
}

const size_t reco::SkimEvent::indexByIso(size_t i) const {
    
    if( i >= leps_.size() ) return 9999;
    std::vector<indexValueStruct> a;

    for(size_t j=0;j<leps_.size();++j) a.push_back(indexValueStruct(allIso(j),j));
    std::sort(a.begin(),a.end(),lowToHigh);

    return a[i].index;
}

const bool reco::SkimEvent::isEB(size_t i) const {

    if( fabs(leps_[i].pdgId()) == 11 ) {
        return static_cast<const pat::Electron&>(leps_[i]).isEB();
    } else {
        return false;
    }
}

const bool reco::SkimEvent::isEE(size_t i) const {

    if( fabs(leps_[i].pdgId()) == 11 ) {
        return static_cast<const pat::Electron&>(leps_[i]).isEE();
    } else {
        return false;
    }
}

const float reco::SkimEvent::tkPt(size_t i) const {

    if( fabs(leps_[i].pdgId()) == 11 ) {
        return static_cast<const pat::Electron&>(leps_[i]).gsfTrack()->pt();
    } else if( fabs(leps_[i].pdgId()) == 13 && !isSTA(i)) {
        return static_cast<const pat::Muon&>(leps_[i]).track()->pt();
    } else {
        return -9999.9;
    }

}

const bool reco::SkimEvent::passesVtxSel(size_t i) const {
    if(i >= vtxs_.size()) {
        return false;
    } else {
        return ( vtxs_[i]->isValid() && !vtxs_[i]->isFake()         && 
                 vtxs_[i]->ndof()                            >= 4.0 && 
            fabs(vtxs_[i]->position().Rho())                  < 2.0 &&
            fabs(vtxs_[i]->z())                               < 24.0  );
    }
    return false;
}

const int reco::SkimEvent::nGoodVertices() const {

    int count = 0;
    for(size_t i=0;i<vtxs_.size();++i) 
        if( passesVtxSel(i) ) count++;

    return count;
}

const bool reco::SkimEvent::hasGoodVertex() const {

    return (nGoodVertices() > 0);
}

const reco::Vertex reco::SkimEvent::highestPtVtx() const {
    if(vtxs_.size() == 0) return reco::Vertex();
    if(sumPts_.size() == 0) return *vtxs_[0];
    double sum = 0;
    size_t high = 0;
    for(size_t i=0;i<vtxs_.size();++i) {
        if( sumPts_[i] > sum && passesVtxSel(i) ) {
            high = i;
            sum = sumPts_[i];
        }
    }
    return *vtxs_[high];
}
        

const double reco::SkimEvent::d0Reco(size_t i) const {
  double dxyPV = 9999;
  if( fabs(leps_[i].pdgId()) == 11 ) {
    const pat::Electron & e = static_cast<const pat::Electron&>(leps_[i]);
    dxyPV = e.userFloat("dxyPV");
    
  } else if( fabs(leps_[i].pdgId()) == 13 && !isSTA(i) ) {
    const pat::Muon & m = static_cast<const pat::Muon&>(leps_[i]);
    dxyPV = m.userFloat("dxyPV");
  } 
  return dxyPV;
}

const double reco::SkimEvent::dZReco(size_t i) const {
  double dzPV = 9999;
  if( fabs(leps_[i].pdgId()) == 11 ) {
    const pat::Electron & e = static_cast<const pat::Electron&>(leps_[i]);
    dzPV = e.userFloat("dzPV");
  } else if( fabs(leps_[i].pdgId()) == 13 && !isSTA(i) ) {
    const pat::Muon & m = static_cast<const pat::Muon&>(leps_[i]);
    dzPV = m.userFloat("dzPV");
  } 

  return dzPV;
}




const bool reco::SkimEvent::passesConversion(size_t i) const {
  if( fabs(leps_[i].pdgId()) == 11 ) {
    const pat::Electron & e = static_cast<const pat::Electron&>(leps_[i]);
    if( fabs(e.userFloat("convValueMapProd:dist")) < 0.02 &&
	fabs(e.userFloat("convValueMapProd:dcot")) < 0.02 ) {
      return false;
    }

    if( e.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() > 0 ) {
      return false;
    }
    return true;
  } else if( fabs(leps_[i].pdgId()) == 13 ) {
    return true;
  } else {
    return false;
  }

}

const bool reco::SkimEvent::isSTA(size_t i) const {
  if( fabs(leps_[i].pdgId()) == 13 ) {
    const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
    return (mu.type() == 8);
  } else {
    return false;
  }
}


const int reco::SkimEvent::bTaggedJetsUnder(const float& maxPt, const float& cut, std::string discriminator) const {

    int count=0;

    for(size_t i=0;i<tagJets_.size();++i) {      
      if( tagJetPt(i,true) > maxPt ) continue;
      if( tagJets_[i]->bDiscriminator(discriminator) <= cut ) continue;	
      if(isThisJetALepton(tagJets_[i])) continue;
      count++;
    }

    return count;
}

const int reco::SkimEvent::bTaggedJetsOver(const float& maxPt, const float& cut, std::string discriminator) const {

    int count=0;

    for(size_t i=0;i<tagJets_.size();++i) {
	if( tagJetPt(i,true) < maxPt ) continue;
        if( tagJets_[i]->bDiscriminator(discriminator) <= cut ) continue;
	if(isThisJetALepton(tagJets_[i])) continue;
        count++;
    }

    return count;
}

const bool reco::SkimEvent::isMuTriggered(size_t i) const {

    if( fabs(leps_[i].pdgId()) == 13 ) {
        const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
        return ( !mu.triggerObjectMatchesByPath("HLT_Mu9").empty() || !mu.triggerObjectMatchesByPath("HLT_Mu15_v1").empty() );
    } else {
        return false;
    }


}


const bool reco::SkimEvent::isEcalSeeded(size_t i) const {

    if( fabs(leps_[i].pdgId()) == 11 ) {
        const pat::Electron &mu = static_cast<const pat::Electron&>(leps_[i]);
        return mu.isMomentumCorrected();
    } else {
        return true;
    }

}


const int reco::SkimEvent::mitType() const {

    if( abs(pdgIdByPt(0)) == 11 ) {
        if( abs(pdgIdByPt(1)) == 11 ) return 11;
        if( abs(pdgIdByPt(1)) == 13 ) return 12;
    } else if( abs(pdgIdByPt(0)) == 13 ) {
        if( abs(pdgIdByPt(1)) == 11 ) return 13;
        if( abs(pdgIdByPt(1)) == 13 ) return 10;
    } else {
        return -1;
    }

}


const float reco::SkimEvent::nearestJet(int i,float minPt, float eta, bool applyCorrection) const {

    if (i >= (int)leps_.size() || i < -1) return -9999.9;

    float dR = 9999;
    for(size_t j=0;j<jets_.size();++j) {
        if( std::fabs(jets_[j]->eta()) >= eta) continue;
        if( jetPt(j,applyCorrection) <= minPt) continue;
        if(isThisJetALepton(jets_[j]))  continue;
    
        float tempdR;
        if(i != -1) {
            tempdR = fabs(ROOT::Math::VectorUtil::DeltaR(jets_[j]->p4(),leps_[i].p4()) );
            if( tempdR < dR ) dR = tempdR;
        } else {
            for(size_t k=0; k<leps_.size();++k){
               tempdR = fabs(ROOT::Math::VectorUtil::DeltaR(jets_[j]->p4(),leps_[k].p4()) );
               if( tempdR < dR ) dR = tempdR;
            }
        }
    }
    return dR;
}


// ============== Matt's methods ===============================
const int reco::SkimEvent::nExtraLepMatt(float minPt) const { 
    int count = 0;
    for(size_t i=0;i<extraLeps_.size();++i) {
        if( extraLeps_[i].pt() <= minPt) continue;
        if( fabs(extraLeps_[i].pdgId()) == 11 ) {
            pat::Electron e = static_cast<const pat::Electron&>(extraLeps_[i]);
            if( fabs(e.eta()) >= 2.5 ) continue;
            if( fabs(e.userFloat("dxyPV")) >= 0.020 ) continue;
            if( fabs(e.userFloat("dzPV") ) >= 1.0 ) continue;
            if( fabs(e.userFloat("convValueMapProd:dist")) < 0.2 &&
                fabs(e.userFloat("convValueMapProd:dcot")) < 0.2 ) continue;
            if( e.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() > 0 ) continue;
            if(!(( e.isEB() && e.sigmaIetaIeta() < 0.01 &&
                  fabs(e.deltaPhiSuperClusterTrackAtVtx()) < 0.06 &&
                  fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.004 &&
//                   e.hadronicOverEm() < 0.04 && 
                  (e.dr03TkSumPt() + std::max(e.dr03EcalRecHitSumEt()-1,(float)0) + e.dr03HcalTowerSumEt() - e.userFloat("rhoEl")*3.14159*0.3*0.3)/e.pt() < 0.1) ||
                  ( !e.isEB() && e.sigmaIetaIeta() < 0.03  && 
                  fabs(e.deltaPhiSuperClusterTrackAtVtx()) < 0.03 &&
                  fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.007 &&
                  (e.dr03TkSumPt() + e.dr03EcalRecHitSumEt() + e.dr03HcalTowerSumEt() - e.userFloat("rhoEl")*3.14159265*0.3*0.3)/e.pt() < 0.1 //&&
                  /*e.hadronicOverEm() < 0.025 */))) continue;
        } else if ( fabs(extraLeps_[i].pdgId()) == 13 ) {
            pat::Muon m = static_cast<const pat::Muon&>(extraLeps_[i]);
            if( fabs(m.userFloat("dxyPV")) >= 0.020 ) continue;
            if( fabs(m.userFloat("dzPV") ) >= 1.0 ) continue;
            if( !(m.isGlobalMuon() && m.isTrackerMuon() && 
                  m.innerTrack()->found() > 10 &&
                  m.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
                  m.globalTrack()->normalizedChi2() < 10 &&
                  m.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
                  m.numberOfMatches() > 1 && fabs(m.track()->ptError() / m.pt()) < 0.10 )) continue;
            if( (m.isolationR03().emEt  + m.isolationR03().hadEt + m.isolationR03().sumPt - m.userFloat("rhoMu") * 3.14159265 * 0.3 * 0.3) / m.pt() >= 0.15 ) continue;
        }
        count++;
    }
    return count;
} 



const int reco::SkimEvent::nSoftMuMatt(float minPt) const { 
  int count = 0;
  for(size_t i=0;i<softMuons_.size();++i) {
    if(softMuons_[i].pt() <= minPt ) continue;
    if( (softMuons_[i].isolationR03().emEt  +
	 softMuons_[i].isolationR03().hadEt +
	 softMuons_[i].isolationR03().sumPt ) / softMuons_[i].pt() < 0.1 &&
	softMuons_[i].pt() > 20.) continue;
    count++;
  }
  return count;
} 


const bool reco::SkimEvent::passesIDV1(size_t i) const {

    if( fabs(leps_[i].pdgId()) == 11 ) {
        const pat::Electron & e = static_cast<const pat::Electron&>(leps_[i]);
        return (( e.isEB() && e.sigmaIetaIeta() < 0.01 &&
                fabs(e.deltaPhiSuperClusterTrackAtVtx()) < 0.06 &&
                fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.004 //&&
                /*e.hadronicOverEm() < 0.04*/) ||
                ( !e.isEB() && e.sigmaIetaIeta() < 0.03  && 
                fabs(e.deltaPhiSuperClusterTrackAtVtx()) < 0.03 &&
                fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.007 //&&
                /*e.hadronicOverEm() < 0.025*/ ));
    } else if( fabs(leps_[i].pdgId()) == 13 ) {
        const pat::Muon & m = static_cast<const pat::Muon&>(leps_[i]);
        return (m.isGlobalMuon() && m.isTrackerMuon() && 
                m.innerTrack()->found() > 10 &&
                m.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
//                 ( (float)m.innerTrack()->hitPattern().trackerLayersWithMeasurement() / (float)(
//                     m.innerTrack()->hitPattern().trackerLayersWithoutMeasurement() + 
//                     m.innerTrack()->hitPattern().trackerLayersWithMeasurement() + 
//                     m.innerTrack()->trackerExpectedHitsInner().numberOfLostHits() +
//                     m.innerTrack()->trackerExpectedHitsOuter().numberOfLostHits() ) >= 0.75 ) &&
//                 ( (float)m.innerTrack()->hitPattern().numberOfValidHits() /  (float)
//                     (m.innerTrack()->hitPattern().numberOfHits() - 
//                      m.innerTrack()->hitPattern().numberOfInactiveHits()) > 0.92 ) &&
                m.globalTrack()->normalizedChi2() < 10 &&
                m.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
                m.numberOfMatches() > 1 && fabs(m.track()->ptError() / m.pt()) < 0.10 );
    } else {
        return false;
    }

}



// ============== Boris' shit ===============================
const bool reco::SkimEvent::passMuID0() const {
  bool answer=true;
  for(unsigned int i=0; i<=1; i++){
    if( fabs(leps_[i].pdgId()) == 13 ) {
      const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
      answer = answer && (mu.genParticleRef().isNonnull() && fabs(mu.genParticleRef()->pdgId())==13);
    }
  }
  return answer;
}


const bool reco::SkimEvent::passMuID1() const {
  bool answer=true;
  for(unsigned int i=0; i<=1; i++){
    if( fabs(leps_[i].pdgId()) == 13 ) {
      const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
      answer = answer && (mu.isGlobalMuon() || (mu.isTrackerMuon() && mu.muonID("TMLastStationTight")));
    }
  }
  return answer;
}

const bool reco::SkimEvent::passMuID2() const {
  bool answer=true;
  for(unsigned int i=0; i<=1; i++){
    if( fabs(leps_[i].pdgId()) == 13 ) {
      const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
      answer = answer && (mu.isTrackerMuon());
    }
  }
  return answer;
}

const bool reco::SkimEvent::passMuID3() const {
  bool answer=true;
  for(unsigned int i=0; i<=1; i++){
    if( fabs(leps_[i].pdgId()) == 13 ) {
      const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
      answer = answer && (mu.innerTrack()->found()>10);
    }
  }
  return answer;
}

const bool reco::SkimEvent::passMuID4() const {
  bool answer=true;
  for(unsigned int i=0; i<=1; i++){
    if( fabs(leps_[i].pdgId()) == 13 ) {
      const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
      answer = answer && (mu.innerTrack()->hitPattern().numberOfValidPixelHits()>0);
    }
  }
  return answer;
}

const bool reco::SkimEvent::passMuID5() const {
  bool answer=true;
  for(unsigned int i=0; i<=1; i++){
    if( fabs(leps_[i].pdgId()) == 13 ) {
      const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
      if(mu.isGlobalMuon())
	answer = answer && (mu.globalTrack()->normalizedChi2()<10);
    }
  }
  return answer;
}

const bool reco::SkimEvent::passMuID6() const {
  bool answer=true;
  for(unsigned int i=0; i<=1; i++){
    if( fabs(leps_[i].pdgId()) == 13 ) {
      const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
      if(mu.isGlobalMuon())
	answer = answer && (mu.globalTrack()->hitPattern().numberOfValidMuonHits()>0);
    }
  }
  return answer;
}

const bool reco::SkimEvent::passMuID7() const {
  bool answer=true;
  for(unsigned int i=0; i<=1; i++){
    if( fabs(leps_[i].pdgId()) == 13 ) {
      const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
      if(!mu.isGlobalMuon())
	answer = answer && (mu.numberOfMatches()>1);
    }
  }
  return answer;
}

const bool reco::SkimEvent::passMuID8() const {
  bool answer=true;
  for(unsigned int i=0; i<=1; i++){
    if( fabs(leps_[i].pdgId()) == 13 ) {
      const pat::Muon &mu = static_cast<const pat::Muon&>(leps_[i]);
      answer = answer && (mu.track()->ptError()/mu.pt()<0.10);
    }
  }
  return answer;
}
