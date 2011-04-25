#include "WWAnalysis/AnalysisStep/interface/CreateEventHists.h"

#include "DataFormats/Common/interface/Handle.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent.h"


CreateEventHists::CreateEventHists(const edm::ParameterSet& cfg, TFileDirectory& fs): 
        edm::BasicAnalyzer::BasicAnalyzer(cfg, fs),
        eventFiller_(cfg.getParameter<edm::ParameterSet>("histParams"), fs, cfg.getParameter<std::string>("sampleName")),
        selectionParams_(cfg.getParameter<edm::ParameterSet>("selectionParams")),
        myLumiSel_(cfg) {

    hypoTypes_ = selectionParams_.getParameterNamesForType<edm::ParameterSet>();

    std::vector<std::string> cutLabels;
    cutLabels.push_back("10/10");
    cutLabels.push_back("20/10");
//     cutLabels.push_back("Lepton ID");
//     cutLabels.push_back("Lepton ISO");
//     cutLabels.push_back("Conversion");
//     cutLabels.push_back("d0/dZ");
    cutLabels.push_back("Lepton Veto");
    cutLabels.push_back("#slash{E}_{T}");
    cutLabels.push_back("m_{ll} > 12");
    cutLabels.push_back("Z-veto");
    cutLabels.push_back("p#slash{E}_{T}");
    cutLabels.push_back("Jet Veto");
    cutLabels.push_back("Soft #mu Veto");
    cutLabels.push_back("b-tag Veto");
    cutLabels.push_back("m_{ll} < 50");
    cutLabels.push_back("p_{T}^{MAX} > 30");
    cutLabels.push_back("p_{T}^{MIN} > 25");
    cutLabels.push_back("#Delta#phi < 60");

    eventFiller_.setCutLabels(cutLabels);
    eventFiller_.setTotalNumberOfCuts(cutLabels.size());

}


void CreateEventHists::analyze(const edm::EventBase& evt) {

    edm::Handle<std::vector<reco::SkimEvent> > skimH;
    if( myLumiSel_(evt) ) for(size_t hypoI=0;hypoI<hypoTypes_.size();++hypoI) {

        // get the skim evt from the file
        evt.getByLabel(hypoTypes_[hypoI],skimH);

        // get the cut variables
        edm::ParameterSet cuts = selectionParams_.getParameter<edm::ParameterSet>(hypoTypes_[hypoI]);

        //loop on no of hypos in evt
        for(std::vector<reco::SkimEvent>::const_iterator mySkimEvent = skimH->begin(); mySkimEvent != skimH->end(); mySkimEvent++){

            size_t instance = mySkimEvent - skimH->begin();
            eventFiller_(&evt,hypoTypes_[hypoI],instance,SKIMMED,*mySkimEvent);

            if( mySkimEvent->q(0)*mySkimEvent->q(1) < 0 && 
                !mySkimEvent->isSTA(0) && !mySkimEvent->isSTA(1) && 
                mySkimEvent->leptEtaCut(cuts.getParameter<double>("etaMu"), cuts.getParameter<double>("etaEl")) && 
                mySkimEvent->ptMin() > cuts.getParameter<double>("ptMin") && 
                mySkimEvent->ptMax() > cuts.getParameter<double>("ptMax") ) {

                eventFiller_(&evt,hypoTypes_[hypoI],instance,FIDUCIAL,*mySkimEvent);
            }

//             if( mySkimEvent->passesIDV1(0) && mySkimEvent->passesIDV1(1) ) {
//                 eventFiller_(&evt,hypoTypes_[hypoI],instance,ID,*mySkimEvent);
//             }
// 
//             if( !(abs(mySkimEvent->pdgId(0)) == 11 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= cuts.getParameter<double>("isoEl")) && 
//                 !(abs(mySkimEvent->pdgId(1)) == 11 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= cuts.getParameter<double>("isoEl")) && 
//                 !(abs(mySkimEvent->pdgId(0)) == 13 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= cuts.getParameter<double>("isoMu")) && 
//                 !(abs(mySkimEvent->pdgId(1)) == 13 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= cuts.getParameter<double>("isoMu")) ) {
// 
//                 eventFiller_(&evt,hypoTypes_[hypoI],instance,ISO,*mySkimEvent);
// 
//             }
// 
//             if( mySkimEvent->passesConversion(0) && mySkimEvent->passesConversion(1) ) {
//                 eventFiller_(&evt,hypoTypes_[hypoI],instance,CONVERSION,*mySkimEvent);
//             }
// 
//             if( fabs(mySkimEvent->d0Reco(0)) < cuts.getParameter<double>("d0") && fabs(mySkimEvent->d0Reco(1)) < cuts.getParameter<double>("d0") &&
//                 fabs(mySkimEvent->dZReco(0)) < cuts.getParameter<double>("dZ") && fabs(mySkimEvent->dZReco(1)) < cuts.getParameter<double>("dZ") ) {
// 
//                 eventFiller_(&evt,hypoTypes_[hypoI],instance,IP,*mySkimEvent);
//             }

            if( mySkimEvent->nExtraLep(10) <= cuts.getParameter<int>("nExtraLep") ) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,LEPTONVETO,*mySkimEvent);
            }

            if( mySkimEvent->tcMet() > cuts.getParameter<double>("met")) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,MET,*mySkimEvent);
            }

            if( mySkimEvent->mll()>cuts.getParameter<double>("mll")) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,MLL,*mySkimEvent);
            }

            if( fabs(mySkimEvent->mll()-91.1876)>cuts.getParameter<double>("mZ") ) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,MZ,*mySkimEvent);
            }

            if( mySkimEvent->projTcMet() > cuts.getParameter<double>("pMet") ) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,PROJMET,*mySkimEvent);
            }

            if( mySkimEvent->nCentralJets( cuts.getParameter<double>("jetPt"), cuts.getParameter<double>("jetEta"), cuts.getParameter<bool>("useJEC")) <= cuts.getParameter<int>("nCentralJet")) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,JETVETO,*mySkimEvent);
            }

            if( mySkimEvent->nSoftMu(3.) <= cuts.getParameter<int>("nSoftMu") ) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,SOFTMU,*mySkimEvent);
            }

            if( mySkimEvent->bTaggedJetsUnder( cuts.getParameter<double>("jetPt"), cuts.getParameter<double>("bValue")) <= cuts.getParameter<int>("nBtagJets") ) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,BJETS,*mySkimEvent);
            }

            if(mySkimEvent->mll() < cuts.getParameter<double>("mllMaxFinal") ) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,MLL2,*mySkimEvent);
            }

            if(mySkimEvent->ptMax() > cuts.getParameter<double>("ptMaxFinal") ) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,PTMAX,*mySkimEvent);
            }

            if(mySkimEvent->ptMin() > cuts.getParameter<double>("ptMinFinal") ) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,PTMIN,*mySkimEvent);
            }

            if(mySkimEvent->dPhill() < cuts.getParameter<double>("deltaPhiLL") ) {
                eventFiller_(&evt,hypoTypes_[hypoI],instance,DPHI,*mySkimEvent);
            }


        } //end of loop over SkimEvent

    } // end loop over hypothesis

}

void CreateEventHists::endJob() {
    eventFiller_.writeAllYieldHists();
    eventFiller_.writeAllNMinus1Plots();
    eventFiller_.writeAllByCutPlots();
}
