#ifndef WWANALYSIS_ANALYSISSTEP_LATINOSREWIEGHT_H
#define WWANALYSIS_ANALYSISSTEP_LATINOSREWIEGHT_H

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Provenance/interface/ProcessHistory.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "FWCore/Common/interface/EventBase.h"

namespace edm {
    class LatinoReWeighting : public LumiReWeighting {
        public:
            LatinoReWeighting( std::string generatedFile,
                    std::string dataFile,
                    std::string GenHistName,
                    std::string DataHistName)
                :LumiReWeighting(generatedFile,dataFile,GenHistName,DataHistName) {}
            LatinoReWeighting( std::vector< double > MC_distr, std::vector< double > Lumi_distr);
            LatinoReWeighting( std::vector< float > MC_distr, std::vector< float > Lumi_distr) 
                :LumiReWeighting(MC_distr,Lumi_distr) {}
            LatinoReWeighting ( ) { } 

            double fullWeight( const edm::EventBase &e ) {

                edm::ReleaseVersion TargetRelease = edm::ReleaseVersion("\"CMSSW_4_2_2_patch2\"");
                std::cout << "Target: " << TargetRelease << std::endl;
                edm::ProcessHistory PHist = e.processHistory();
                edm::ProcessHistory::const_iterator PHist_iter = PHist.begin();

                bool isS3 = false;
                for(; PHist_iter<PHist.end() ;++PHist_iter) {
                    edm::ProcessConfiguration PConf = *(PHist_iter);
                    edm::ReleaseVersion Release =  PConf.releaseVersion() ;
                    const std::string Process =  PConf.processName();

                    std::cout << "Process: " << Process << std::endl;
                    std::cout << "Release: " << Release << std::endl;
                    if((Process=="HLT") && (Release==TargetRelease)) {
                        isS3 = true;
                        break;
                    }
                }

                if(isS3) return weightOOT(e);
                else     return weight3BX(e);

            }

    };
    LatinoReWeighting::LatinoReWeighting( std::vector< double > MC_distr, std::vector< double > Lumi_distr) {

        if( MC_distr.size() != Lumi_distr.size() ){   

            std::cerr <<"ERROR: LumiReWeighting: input vectors have different sizes. Quitting... \n";
            return;

        }

        Int_t NBins = MC_distr.size();

        weights_ = boost::shared_ptr<TH1> ( new TH1F("luminumer","luminumer",NBins,-0.5, float(NBins)-0.5) );
        TH1* den = new TH1F("lumidenom","lumidenom",NBins,-0.5, float(NBins)-0.5) ;

        for(int ibin = 1; ibin<NBins+1; ++ibin ) {
            weights_->SetBinContent(ibin, Lumi_distr[ibin-1]);
            den->SetBinContent(ibin,MC_distr[ibin-1]);
        }

        // check integrals, make sure things are normalized

        float deltaH = weights_->Integral();
        if(fabs(1.0 - deltaH) > 0.02 ) { //*OOPS*...
            weights_->Scale( 1.0/ weights_->Integral() );
        }
        float deltaMC = den->Integral();
        if(fabs(1.0 - deltaMC) > 0.02 ) {
            den->Scale(1.0/ den->Integral());
        }

        weights_->Divide( den );  // so now the average weight should be 1.0    

        std::cout << " Lumi/Pileup Reweighting: Computed Weights per In-Time Nint " << std::endl;

        for(int ibin = 1; ibin<NBins+1; ++ibin){
            std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
        }

        weightOOT_init();

        FirstWarning_ = true;
        OldLumiSection_ = -1;

    }

}

#endif
