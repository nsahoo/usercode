#include "WWAnalysis/AnalysisStep/interface/SkimTupleProducer.h"


SkimTupleProducer::SkimTupleProducer(const edm::ParameterSet& cfg) : 
        skimTag_(cfg.getUntrackedParameter<edm::InputTag>("src")), 
        plotsInfo_(cfg.getUntrackedParameter<std::vector<edm::ParameterSet> >("variables")),
        weight_(cfg.getUntrackedParameter<double>("weight")) {

    fs_ = edm::Service<TFileService>().operator->();
    tree_ = fs_->make<TTree>(cfg.getParameter<std::string>("@module_label").c_str(),cfg.getParameter<std::string>("@module_label").c_str());
    tree_->Branch("weight",&weight_,"weight/D");

}

void SkimTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<std::vector< reco::SkimEvent> > skimH;
    iEvent.getByLabel(skimTag_,skimH);
    

    //Loop Over Each Plot Type
    for(size_t j=0;j<plotsInfo_.size();++j) {
        placeHolder_[j] = objFuncs_[j](skimH->at(0));
    }            

    tree_->Fill();
}

void SkimTupleProducer::beginJob() { 

    placeHolder_.clear();
    objFuncs_.clear();

    placeHolder_.resize(plotsInfo_.size());

    for(size_t j=0;j<plotsInfo_.size();++j) {
        objFuncs_.push_back( StringObjectFunction<reco::SkimEvent>(
                             plotsInfo_[j].getUntrackedParameter<std::string>("quantity")));
        tree_->Branch( plotsInfo_[j].getUntrackedParameter<std::string>("tag").c_str() , 
                       &placeHolder_[j], 
                       (plotsInfo_[j].getUntrackedParameter<std::string>("tag")+"/D").c_str() );
    }

}

SkimTupleProducer::~SkimTupleProducer() { } 
void SkimTupleProducer::endJob() { }
DEFINE_FWK_MODULE(SkimTupleProducer);

