#include <memory>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class CombinedWeightProducer : public edm::EDProducer {
    public:
        explicit CombinedWeightProducer(const edm::ParameterSet&);
        ~CombinedWeightProducer();

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);

        bool hasHiggs_;
        edm::InputTag higgsTag_;
        edm::InputTag puTag_;
        std::vector<double> puWeights_;
};

CombinedWeightProducer::CombinedWeightProducer(const edm::ParameterSet& iConfig) :
    hasHiggs_(iConfig.existsAs<edm::InputTag>("ptWeight")),
    higgsTag_(hasHiggs_ ? iConfig.getParameter<edm::InputTag>("ptWeight") : edm::InputTag()),
    puTag_(iConfig.getParameter<edm::InputTag>("puLabel")),
    puWeights_(iConfig.getParameter<std::vector<double> >("puWeight"))
{ 
    produces<double>();
} 


void CombinedWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    std::auto_ptr<double> weight(new double(1));

    edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
    iEvent.getByLabel(puTag_,puInfoH);
    int nPU = 0;
    for(size_t i=0;i<puInfoH->size();++i) {
        if( puInfoH->at(i).getBunchCrossing()==0 ) {
            nPU = puInfoH->at(i).getPU_NumInteractions();
        }
    }
    *weight *= puWeights_[std::min(nPU,(int)(puWeights_.size()-1))];

    if (hasHiggs_) {
        edm::Handle<double> ptWeight;
        iEvent.getByLabel(higgsTag_,ptWeight);
        *weight *= *ptWeight;
    }

    iEvent.put(weight);
}

CombinedWeightProducer::~CombinedWeightProducer() { }
DEFINE_FWK_MODULE(CombinedWeightProducer);


