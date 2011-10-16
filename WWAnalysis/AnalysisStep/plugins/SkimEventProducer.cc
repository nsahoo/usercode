#include "WWAnalysis/AnalysisStep/interface/SkimEventProducer.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include<vector> 
#include "Math/VectorUtil.h"


SkimEventProducer::SkimEventProducer(const edm::ParameterSet& cfg) :
    //   branchAlias_(cfg.getParameter<std::string>("branchAlias")),
    hypoType_(reco::SkimEvent::hypoTypeByName(cfg.getParameter<std::string>("hypoType"))),
    triggerTag_(cfg.getParameter<edm::InputTag>("triggerTag")),
    singleMuData_  ( cfg.getParameter<std::vector<std::string> >("singleMuDataPaths") ),
    singleElData_  ( cfg.getParameter<std::vector<std::string> >("singleElDataPaths") ),
    doubleMuData_  ( cfg.getParameter<std::vector<std::string> >("doubleMuDataPaths") ),
    doubleElData_  ( cfg.getParameter<std::vector<std::string> >("doubleElDataPaths") ),
    muEGData_      ( cfg.getParameter<std::vector<std::string> >("muEGDataPaths") ),
    singleMuMC_    ( cfg.getParameter<std::vector<std::string> >("singleMuMCPaths") ),
    singleElMC_    ( cfg.getParameter<std::vector<std::string> >("singleElMCPaths") ),
    doubleMuMC_    ( cfg.getParameter<std::vector<std::string> >("doubleMuMCPaths") ),
    doubleElMC_    ( cfg.getParameter<std::vector<std::string> >("doubleElMCPaths") ),
    muEGMC_        ( cfg.getParameter<std::vector<std::string> >("muEGMCPaths") )
{
    muTag_          = cfg.getParameter<edm::InputTag>("muTag"     ); 
    elTag_          = cfg.getParameter<edm::InputTag>("elTag"     ); 
    softMuTag_      = cfg.getParameter<edm::InputTag>("softMuTag" ); 
    jetTag_         = cfg.getParameter<edm::InputTag>("jetTag"    ); 
    tagJetTag_      = cfg.getParameter<edm::InputTag>("tagJetTag" ); 
    pfMetTag_       = cfg.getParameter<edm::InputTag>("pfMetTag"  ); 
    tcMetTag_       = cfg.getParameter<edm::InputTag>("tcMetTag"  ); 
    chargedMetTag_  = cfg.getParameter<edm::InputTag>("chargedMetTag" ); 
    vtxTag_         = cfg.getParameter<edm::InputTag>("vtxTag"        ); 
    chCandsTag_     = cfg.getParameter<edm::InputTag>("chCandsTag"    ); 

    if (cfg.exists("sptTag"    )) sptTag_     = cfg.getParameter<edm::InputTag>("sptTag"    ); 
    else                          sptTag_     = edm::InputTag("","","");
    if (cfg.exists("spt2Tag"   )) spt2Tag_    = cfg.getParameter<edm::InputTag>("spt2Tag"   ); 
    else                          spt2Tag_    = edm::InputTag("","","");

    selectorMuTight = new StringCutObjectSelector<reco::Candidate>(cfg.getParameter<std::string>("tightMuSelection"),true);
    selectorMuLoose = new StringCutObjectSelector<reco::Candidate>(cfg.getParameter<std::string>("looseMuSelection"),true);
    selectorEleTight = new StringCutObjectSelector<reco::Candidate>(cfg.getParameter<std::string>("tightEleSelection"),true);
    selectorEleLoose = new StringCutObjectSelector<reco::Candidate>(cfg.getParameter<std::string>("looseEleSelection"),true);


    // ---- Here I initialize the BDT for the mva ElectronID
    eleBDT  = new ElectronIDMVA();
    const char *base=getenv("CMSSW_BASE");
    std::string baseFolder(base);
    baseFolder += "/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/"; 
    eleBDT->Initialize("BDTG method",
		       (baseFolder+"Subdet0LowPt_WithIPInfo_BDTG.weights.xml").c_str(),
		       (baseFolder+"Subdet1LowPt_WithIPInfo_BDTG.weights.xml").c_str(),
		       (baseFolder+"Subdet2LowPt_WithIPInfo_BDTG.weights.xml").c_str(),
		       (baseFolder+"Subdet0HighPt_WithIPInfo_BDTG.weights.xml").c_str(),
		       (baseFolder+"Subdet1HighPt_WithIPInfo_BDTG.weights.xml").c_str(),
		       (baseFolder+"Subdet2HighPt_WithIPInfo_BDTG.weights.xml").c_str(),                
		       ElectronIDMVA::kWithIPInfo);
    // ---------

    produces<std::vector<reco::SkimEvent> >().setBranchAlias(cfg.getParameter<std::string>("@module_label"));
}

void SkimEventProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<reco::SkimEvent> > skimEvent(new std::vector<reco::SkimEvent> );

    //SkimEvent::hypoType type = SkimEvent::hypoTypeByName(hypoType_);

    edm::Handle<pat::JetCollection> jetH;
    iEvent.getByLabel(jetTag_,jetH);

    edm::Handle<pat::JetCollection> tagJetH;
    if(!(tagJetTag_==edm::InputTag(""))) iEvent.getByLabel(tagJetTag_,tagJetH);

    edm::Handle<reco::PFMETCollection> pfMetH;
    iEvent.getByLabel(pfMetTag_,pfMetH);

    edm::Handle<reco::METCollection> tcMetH;
    iEvent.getByLabel(tcMetTag_,tcMetH);

    edm::Handle<edm::ValueMap<reco::PFMET> > chargedMetH;
    iEvent.getByLabel(chargedMetTag_,chargedMetH);

    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vtxTag_,vtxH);

    edm::Handle<reco::CandidateCollection> candsH;
    iEvent.getByLabel(chCandsTag_,candsH);

    edm::Handle<edm::ValueMap<float> > sptH;
    if(!(sptTag_==edm::InputTag(""))) iEvent.getByLabel(sptTag_,sptH);

    edm::Handle<edm::ValueMap<float> > spt2H;
    if(!(spt2Tag_==edm::InputTag(""))) iEvent.getByLabel(spt2Tag_,spt2H);

    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerTag_,triggerResults);

    // May God have mercy on my soul ...
    std::vector<bool> passBits;
    passBits.push_back( singleMuData_.check(iEvent,*triggerResults) );
    passBits.push_back( singleElData_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleMuData_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleElData_.check(iEvent,*triggerResults) );
    passBits.push_back( muEGData_.check(    iEvent,*triggerResults) );
    passBits.push_back( singleMuMC_.check(iEvent,*triggerResults) );
    passBits.push_back( singleElMC_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleMuMC_.check(iEvent,*triggerResults) );
    passBits.push_back( doubleElMC_.check(iEvent,*triggerResults) );
    passBits.push_back( muEGMC_.check(    iEvent,*triggerResults) );


    edm::Handle<edm::View<reco::RecoCandidate> > muons;
    iEvent.getByLabel(muTag_,muons);
    edm::Handle<edm::View<reco::RecoCandidate> > softs;
    iEvent.getByLabel(softMuTag_,softs);
    edm::Handle<edm::View<reco::RecoCandidate> > electrons;
    iEvent.getByLabel(elTag_,electrons);


    //using namespace std;
    //cout << "muon collection ID: " << muons.id() << endl;
    //cout << "ele collection ID: " << electrons.id() << endl;
    //cout << "soft collection ID: " << softs.id() << endl;

    /*
       edm::Handle<pat::ElectronCollection> extraElH;
       iEvent.getByLabel(extraElTag_,extraElH);
     */


    if(hypoType_==reco::SkimEvent::WWELEL){//ELEL
        for(size_t i=0;i<electrons->size();++i) {
            for(size_t j=i+1;j<electrons->size();++j) {
                skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );      
                skimEvent->back().setEventInfo(iEvent);
                // Leptons
		bool passLoose(false);
		bool passTight(false);
		float bdt(-99.);
		passLoose = (*selectorEleLoose)(*(electrons->ptrAt(i)));
		passTight = (*selectorEleTight)(*(electrons->ptrAt(i)));
		bdt = getBdtOutput(electrons,i);
                skimEvent->back().setLepton(electrons,i,passLoose,passTight,bdt);
		passLoose = (*selectorEleLoose)(*(electrons->ptrAt(j)));
		passTight = (*selectorEleTight)(*(electrons->ptrAt(j)));
		bdt = getBdtOutput(electrons,j);
                skimEvent->back().setLepton(electrons,j,passLoose,passTight,bdt);
		

                for(size_t k=0;k<electrons->size();++k) if(k!=i && k!=j) skimEvent->back().setExtraLepton(electrons,k);
                for(size_t k=0;k<muons->size();++k) {
                    float delta1 = ROOT::Math::VectorUtil::DeltaR(muons->at(k).p4(),electrons->at(i).p4());
                    float delta2 = ROOT::Math::VectorUtil::DeltaR(muons->at(k).p4(),electrons->at(j).p4());
                    if(delta1 > 0.1 && delta2 > 0.1) skimEvent->back().setExtraLepton(muons,k);
                }
                for(size_t k=0;k<softs->size();++k) skimEvent->back().setSoftMuon(softs,k);
                // Everything else
                skimEvent->back().setTriggerBits(passBits);
                skimEvent->back().setJets(jetH);
                skimEvent->back().setPFMet(pfMetH);
                skimEvent->back().setTCMet(tcMetH);
                skimEvent->back().setChargedMet(chargedMetH->get(0));
                skimEvent->back().setVertex(vtxH);
                if(sptH.isValid()   ) skimEvent->back().setVtxSumPts(sptH);
                if(spt2H.isValid()  ) skimEvent->back().setVtxSumPt2s(spt2H);
                if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
                else                  skimEvent->back().setTagJets(jetH);
                skimEvent->back().setChargedMetSmurf(doChMET(candsH,&electrons->at(i),&electrons->at(j)));
                //       skimEvent->back().setupJEC(l2File_,l3File_,resFile_);
            }
        }//end loop on main lepton collection
    }else if(hypoType_ == reco::SkimEvent::WWELMU){
        for(size_t i=0;i<electrons->size();++i) {
            for(size_t j=0;j<muons->size();++j) {
                if( muons->at(j).pt() >= electrons->at(i).pt() ) continue;
                skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );      
                skimEvent->back().setEventInfo(iEvent);
                // Leptons
		bool passLoose(false);
		bool passTight(false);
		float bdt(-99.);
		passLoose = (*selectorEleLoose)(*(electrons->ptrAt(i)));
		passTight = (*selectorEleTight)(*(electrons->ptrAt(i)));
		bdt = getBdtOutput(electrons,i);
                skimEvent->back().setLepton(electrons,i,passLoose,passTight,bdt);
		passLoose = (*selectorMuLoose)(*(muons->ptrAt(j)));
		passTight = (*selectorMuTight)(*(muons->ptrAt(j)));
		bdt = getBdtOutput(muons,j);
                skimEvent->back().setLepton(muons,j,passLoose,passTight,bdt);


                for(size_t k=0;k<electrons->size();++k) {
                    float delta1 = ROOT::Math::VectorUtil::DeltaR(electrons->at(k).p4(),muons->at(j).p4());
                    if(k!=i && delta1 > 0.1) 
                        skimEvent->back().setExtraLepton(electrons,k);
                }
                for(size_t k=0;k<muons->size();++k) {
                    float delta1 = ROOT::Math::VectorUtil::DeltaR(muons->at(k).p4(),electrons->at(i).p4());
                    if(delta1 > 0.1 && k!=j) 
                        skimEvent->back().setExtraLepton(muons,k);
                }
                for(size_t k=0;k<softs->size();++k) {
                    if(softs->at(k).pt() != muons->at(j).pt() || softs->at(k).eta() != muons->at(j).eta()) 
                        skimEvent->back().setSoftMuon(softs,k);
                }
                // Everything else
                skimEvent->back().setTriggerBits(passBits);
                skimEvent->back().setJets(jetH);
                skimEvent->back().setPFMet(pfMetH);
                skimEvent->back().setTCMet(tcMetH);
                skimEvent->back().setChargedMet(chargedMetH->get(0));
                skimEvent->back().setVertex(vtxH);
                if(sptH.isValid()   ) skimEvent->back().setVtxSumPts(sptH);
                if(spt2H.isValid()  ) skimEvent->back().setVtxSumPt2s(spt2H);
                if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
                else                  skimEvent->back().setTagJets(jetH);
                skimEvent->back().setChargedMetSmurf(doChMET(candsH,&electrons->at(i),&muons->at(j)));
                //       skimEvent->back().setupJEC(l2File_,l3File_,resFile_);
            }
        }//end loop on main lepton collection
    }else if(hypoType_ == reco::SkimEvent::WWMUEL){
        for(size_t i=0;i<electrons->size();++i) {
            for(size_t j=0;j<muons->size();++j) {
                if( muons->at(j).pt() < electrons->at(i).pt() ) continue;
                skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );      
                skimEvent->back().setEventInfo(iEvent);
                // Leptons
		bool passLoose(false);
		bool passTight(false);
		float bdt(-99.);
		passLoose = (*selectorEleLoose)(*(electrons->ptrAt(i)));
		passTight = (*selectorEleTight)(*(electrons->ptrAt(i)));
		bdt = getBdtOutput(electrons,i);
                skimEvent->back().setLepton(electrons,i,passLoose,passTight,bdt);
		passLoose = (*selectorMuLoose)(*(muons->ptrAt(j)));
		passTight = (*selectorMuTight)(*(muons->ptrAt(j)));
		bdt = getBdtOutput(muons,j);
                skimEvent->back().setLepton(muons,j,passLoose,passTight,bdt);

                for(size_t k=0;k<electrons->size();++k) {
                    float delta1 = ROOT::Math::VectorUtil::DeltaR(electrons->at(k).p4(),muons->at(j).p4());
                    if(k!=i && delta1 > 0.1) 
                        skimEvent->back().setExtraLepton(electrons,k);
                }
                for(size_t k=0;k<muons->size();++k) {
                    float delta1 = ROOT::Math::VectorUtil::DeltaR(muons->at(k).p4(),electrons->at(i).p4());
                    if(delta1 > 0.1 && k!=j) 
                        skimEvent->back().setExtraLepton(muons,k);
                }
                for(size_t k=0;k<softs->size();++k) {
                    if(softs->at(k).pt() != muons->at(j).pt() || softs->at(k).eta() != muons->at(j).eta()) 
                        skimEvent->back().setSoftMuon(softs,k);
                }
                // Everything else
                skimEvent->back().setTriggerBits(passBits);
                skimEvent->back().setJets(jetH);
                skimEvent->back().setPFMet(pfMetH);
                skimEvent->back().setTCMet(tcMetH);
                skimEvent->back().setChargedMet(chargedMetH->get(0));
                skimEvent->back().setVertex(vtxH);
                if(sptH.isValid()   ) skimEvent->back().setVtxSumPts(sptH);
                if(spt2H.isValid()  ) skimEvent->back().setVtxSumPt2s(spt2H);
                if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
                else                  skimEvent->back().setTagJets(jetH);
                skimEvent->back().setChargedMetSmurf(doChMET(candsH,&electrons->at(i),&muons->at(j)));
                //       skimEvent->back().setupJEC(l2File_,l3File_,resFile_);
            }
        }//end loop on main lepton collection
    }else if(hypoType_==reco::SkimEvent::WWMUMU){//MUMU
        for(size_t i=0;i<muons->size();++i) {
            for(size_t j=i+1;j<muons->size();++j) {
                skimEvent->push_back( *(new reco::SkimEvent(hypoType_) ) );      
                skimEvent->back().setEventInfo(iEvent);
                // Leptons
		bool passLoose(false);
		bool passTight(false);
		float bdt(-99.);
		passLoose = (*selectorMuLoose)(*(muons->ptrAt(i)));
		passTight = (*selectorMuTight)(*(muons->ptrAt(i)));
		bdt = getBdtOutput(muons,i);
                skimEvent->back().setLepton(muons,i,passLoose,passTight,bdt);
		passLoose = (*selectorMuLoose)(*(muons->ptrAt(j)));
		passTight = (*selectorMuTight)(*(muons->ptrAt(j)));
		bdt = getBdtOutput(muons,j);
                skimEvent->back().setLepton(muons,j,passLoose,passTight,bdt);

                for(size_t k=0;k<electrons->size();++k) {
                    float delta1 = ROOT::Math::VectorUtil::DeltaR(electrons->at(k).p4(),muons->at(i).p4());
                    float delta2 = ROOT::Math::VectorUtil::DeltaR(electrons->at(k).p4(),muons->at(j).p4());
                    if(delta1 > 0.1 && delta2 > 0.1) 
                        skimEvent->back().setExtraLepton(electrons,k);
                }
                for(size_t k=0;k<muons->size();++k) {
                    if(k!=i && k!=j) 
                        skimEvent->back().setExtraLepton(muons,k);
                }
                for(size_t k=0;k<softs->size();++k) {
                    if( (softs->at(k).pt() != muons->at(i).pt() || softs->at(k).eta() != muons->at(i).eta())  && 
                        (softs->at(k).pt() != muons->at(j).pt() || softs->at(k).eta() != muons->at(j).eta()) ) 
                        skimEvent->back().setSoftMuon(softs,k);
                }
                // Everything else
                skimEvent->back().setTriggerBits(passBits);
                skimEvent->back().setJets(jetH);
                skimEvent->back().setPFMet(pfMetH);
                skimEvent->back().setTCMet(tcMetH);
                skimEvent->back().setChargedMet(chargedMetH->get(0));
                skimEvent->back().setVertex(vtxH);
                if(sptH.isValid()   ) skimEvent->back().setVtxSumPts(sptH);
                if(spt2H.isValid()  ) skimEvent->back().setVtxSumPt2s(spt2H);
                if(tagJetH.isValid()) skimEvent->back().setTagJets(tagJetH);
                else                  skimEvent->back().setTagJets(jetH);
                skimEvent->back().setChargedMetSmurf(doChMET(candsH,&muons->at(i),&muons->at(j)));
                //       skimEvent->back().setupJEC(l2File_,l3File_,resFile_);
            }
        }//end loop on main lepton collection
    }else{
        throw cms::Exception("BadInput") 
            << "ERROR: event type " << hypoType_ << " is not known" ;
    }

    iEvent.put(skimEvent);
}

SkimEventProducer::~SkimEventProducer() { 
  delete selectorMuTight;
  delete selectorMuLoose;
  delete selectorEleTight;
  delete selectorEleLoose;
} 

void SkimEventProducer::beginJob() { }
void SkimEventProducer::endJob() { }



reco::MET SkimEventProducer::doChMET(edm::Handle<reco::CandidateCollection> candsH,
        const reco::Candidate* cand1,const reco::Candidate* cand2){
    using namespace std;
    reco::Candidate::LorentzVector totalP4;
    for(reco::CandidateCollection::const_iterator it= candsH->begin(); it!=candsH->end(); ++it){
        if(fabs(ROOT::Math::VectorUtil::DeltaR(it->p4(),cand1->p4())) <=0.1) continue;
        if(fabs(ROOT::Math::VectorUtil::DeltaR(it->p4(),cand2->p4())) <=0.1) continue;
        totalP4 += it->p4();
    }
    totalP4 +=cand1->p4();
    totalP4 +=cand2->p4();
    reco::Candidate::LorentzVector invertedP4(-totalP4);
    reco::MET met(invertedP4,reco::Candidate::Point(0,0,0));
    return met;
}




float 
SkimEventProducer::getBdtOutput(const edm::Handle<edm::View<reco::RecoCandidate> > &h,size_t i){
  float value(-99.);
  if( fabs(h->ptrAt(i)->pdgId() )== 13){
    //const pat::Muon* muon = static_cast<const pat::Muon*>(h->ptrAt(i).get());

    //Until we don't have a bdg for muons, we always set the bdt value to 1 for the muons; 
    value = 1.0;
  }else if( fabs(h->ptrAt(i)->pdgId() )== 11){
    const pat::Electron* ele = static_cast<const pat::Electron*>(h->ptrAt(i).get());
    double EleOneOverEMinusOneOverP = (1.0/(ele->superCluster()->energy())) - 1.0 / ele->gsfTrack()->p(); 
    double xieSign  = ( (-ele->userFloat("dxyPV")) >=0 )  ? 1: -1;
    value = eleBDT->MVAValue(ele->pt() , 
			     ele->superCluster()->eta(), 
			     ele->sigmaIetaIeta(), 
			     ele->deltaEtaSuperClusterTrackAtVtx(),
			     ele->deltaPhiSuperClusterTrackAtVtx(),
			     ele->hcalOverEcal(),
			     -ele->userFloat("dxyPV"),
			     ele->fbrem(), 
			     ele->eSuperClusterOverP(), 
			     ele->eSeedClusterOverPout(),
			     ele->userFloat("sigmaIphiIphi"),
			     ele->basicClustersSize() - 1,
			     EleOneOverEMinusOneOverP,
			     ele->eSeedClusterOverP(),
			     xieSign*ele->userFloat("ip"),
			     xieSign*ele->userFloat("ip")/ele->userFloat("ipErr"));
  }else{
    std::cout << "LOGIC ERROR: casting to patEle and patMuon failed" << std::endl;
  }
  return value;
}


DEFINE_FWK_MODULE(SkimEventProducer);



