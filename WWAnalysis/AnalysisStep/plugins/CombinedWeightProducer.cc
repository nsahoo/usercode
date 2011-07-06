#include <memory>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"


#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class CombinedWeightProducer : public edm::EDProducer {
    public:
        explicit CombinedWeightProducer(const edm::ParameterSet&);
        ~CombinedWeightProducer();

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);
    
        double baseWeight_;
        bool hasHiggs_;
        edm::InputTag higgsTag_;
        bool hasPU_;
        edm::InputTag puTag_;
        std::vector<double> puWeights_;
        bool hasSrc_;
        edm::InputTag src_;
};

CombinedWeightProducer::CombinedWeightProducer(const edm::ParameterSet& iConfig) :
    baseWeight_(iConfig.existsAs<double>("baseWeight") ? iConfig.getParameter<double>("baseWeight") : 1.0),
    hasHiggs_(iConfig.existsAs<edm::InputTag>("ptWeight")),
    higgsTag_(hasHiggs_ ? iConfig.getParameter<edm::InputTag>("ptWeight") : edm::InputTag()),
    hasPU_(iConfig.existsAs<edm::InputTag>("puLabel")),
    puTag_(hasPU_ ? iConfig.getParameter<edm::InputTag>("puLabel") : edm::InputTag() ),
    puWeights_(hasPU_ ? iConfig.getParameter<std::vector<double> >("puWeight") : std::vector<double>()),
    hasSrc_(iConfig.existsAs<edm::InputTag>("src")),
    src_(hasSrc_ ? iConfig.getParameter<edm::InputTag>("src") : edm::InputTag())
{ 
    if (hasSrc_) produces<edm::ValueMap<float> > ();
    else         produces<double>();
} 


void CombinedWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    std::auto_ptr<double> weight(new double(baseWeight_));

    if(hasPU_) {
        edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
        iEvent.getByLabel(puTag_,puInfoH);
        int nPU = 0;
        for(size_t i=0;i<puInfoH->size();++i) {
            if( puInfoH->at(i).getBunchCrossing()==0 ) {
                nPU = puInfoH->at(i).getPU_NumInteractions();
            }
        }
        *weight *= puWeights_[std::min(nPU,(int)(puWeights_.size()-1))];
    }

    if (hasHiggs_) {
        edm::Handle<double> ptWeight;
        iEvent.getByLabel(higgsTag_,ptWeight);
        *weight *= *ptWeight;
    }

    if(hasSrc_) {
        edm::Handle<edm::View<reco::Candidate> > srcH;
        iEvent.getByLabel(src_ ,srcH);
        std::auto_ptr<edm::ValueMap<float> > floatMap(new edm::ValueMap<float>);
        std::vector<float> floatVec(srcH->size(),*weight);
        edm::ValueMap<float>::Filler filler(*floatMap);
        filler.insert(srcH, floatVec.begin(), floatVec.end());
        filler.fill();
        iEvent.put(floatMap);
    } else {
        iEvent.put(weight);
    }
}

CombinedWeightProducer::~CombinedWeightProducer() { }
DEFINE_FWK_MODULE(CombinedWeightProducer);


