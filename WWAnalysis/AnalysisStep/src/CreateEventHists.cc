#include "WWAnalysis/AnalysisStep/interface/CreateEventHists.h"

#include "DataFormats/Common/interface/Handle.h"


CreateEventHists::CreateEventHists(const edm::ParameterSet& cfg, TFileDirectory& fs): 
        edm::BasicAnalyzer::BasicAnalyzer(cfg, fs),
        eventFiller_(cfg.getParameter<edm::ParameterSet>("histParams"), fs, cfg.getParameter<std::string>("sampleName")),
        myLumiSel_(cfg) {

    // Get a list of all hypotheses
    edm::ParameterSet hypotheses = cfg.getParameter<edm::ParameterSet>("hypotheses");

    // For each PSet in hypotheses, a new hypothesis is made
    hypoNames_ = hypotheses.getParameterNamesForType<edm::ParameterSet>();
    for(size_t hn=0;hn<hypoNames_.size();++hn) {

        // Get the information about this particular hypothesis
        edm::ParameterSet thisPset = hypotheses.getParameter<edm::ParameterSet>(hypoNames_[hn]);

        // Grab the input tag for this hypothesis
        branchTags_.push_back( thisPset.getParameter<edm::InputTag>("src") );

        // Fill the cut labels  and grab the stringselectors for each cut
        std::vector<std::string> cutLabels;
        std::vector<edm::ParameterSet> cuts = thisPset.getParameter<std::vector<edm::ParameterSet> >("cuts");
        for(size_t i=0;i<cuts.size();++i) {
            cutLabels.push_back( cuts[i].getParameter<std::string>("label") );
            stringSelectors_[hypoNames_[hn]].push_back( StringCutObjectSelector<reco::SkimEvent>( cuts[i].getParameter<std::string>("cut") ) );
        }
        eventFiller_.setCutLabels(hypoNames_[hn], cutLabels);
        eventFiller_.setTotalNumberOfCuts(hypoNames_[hn], cutLabels.size());
    }

}


void CreateEventHists::analyze(const edm::EventBase& evt) {

    edm::Handle<std::vector<reco::SkimEvent> > skimH;
    if( myLumiSel_(evt) ) for(size_t hn=0;hn<hypoNames_.size();++hn) {

        // get the skim evt from the file
        evt.getByLabel(branchTags_[hn],skimH);

        //loop on no of hypos in evt
        for(std::vector<reco::SkimEvent>::const_iterator mySkimEvent = skimH->begin(); mySkimEvent != skimH->end(); mySkimEvent++){

            size_t instance = mySkimEvent - skimH->begin();

            // loop over cuts for this hypothesis
            for(size_t i=0;i<stringSelectors_[hypoNames_[hn]].size();++i) {
                if( stringSelectors_[hypoNames_[hn]][i]( (*mySkimEvent) ) ) eventFiller_(&evt,hypoNames_[hn],instance,i,*mySkimEvent);
            }

        } //end of loop over SkimEvent

    } // end loop over hypothesis

}

void CreateEventHists::endJob() {
    eventFiller_.writeAllYieldHists();
    eventFiller_.writeAllNMinus1Plots();
    eventFiller_.writeAllByCutPlots();
}
