#include "WWAnalysis/AnalysisStep/interface/MCDumper.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"



MCDumper::MCDumper(const edm::ParameterSet& cfg) {

    fs_ = edm::Service<TFileService>().operator->();
    tree_ = fs_->make<TTree>(cfg.getParameter<std::string>("@module_label").c_str(),cfg.getParameter<std::string>("@module_label").c_str());
    tree_->Branch("MHiggs" ,&MHiggs_ ,"MHiggs/F" );
}

void MCDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles",genParticles);


    reco::GenParticleRefVector genParticles_;
    genParticles_.clear();
    for(size_t i=0; i<genParticles->size(); ++i) {
     genParticles_.push_back( reco::GenParticleRef(genParticles,i) );
    }


    float mass = -1;

    const reco::Candidate* mcH = 0;

  // loop over gen particles
    for(size_t gp=0; gp<genParticles_.size();++gp){

     int pdgId  = genParticles_[gp] -> pdgId();
     int status = genParticles_[gp] -> status();

    // Stop {1000006}
     if( (pdgId == 25) && (status == 3) ) {
      mcH = &(*(genParticles_[gp]));
      mass = mcH->mass();
     }
    } // loop over gen particles

    MHiggs_= mass;
    tree_->Fill();

}

void MCDumper::beginJob() { }
MCDumper::~MCDumper() { } 
void MCDumper::endJob() { }
DEFINE_FWK_MODULE(MCDumper);

