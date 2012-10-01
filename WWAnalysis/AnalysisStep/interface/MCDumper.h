#ifndef MYPRODUCERS_MCDumper_MCDumper_H
#define MYPRODUCERS_MCDumper_MCDumper_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "WWAnalysis/AnalysisStep/interface/LatinoReWeighting.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Candidate/interface/Candidate.h"



#include <TTree.h>
#include <vector>
#include <string>


class MCDumper : public edm::EDAnalyzer {
    public:
        explicit MCDumper(const edm::ParameterSet&);
        ~MCDumper();

    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        TFileService *fs_;
        TTree *tree_;
        float MHiggs_;
        reco::GenParticleRefVector genParticles_;
};

#endif
