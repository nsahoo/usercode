#include "WWAnalysis/AnalysisStep/interface/EventBitsFiller.h"
#include "TFile.h"
#include "TH1F.h"


EventBitsFiller::AdvancedEventBits EventBitsFiller::AdvancedEventBits::operator&(const EventBitsFiller::AdvancedEventBits &rhs) const {
    return (getBits()&rhs.getBits());
}

EventBitsFiller::AdvancedEventBits EventBitsFiller::AdvancedEventBits::operator|(const EventBitsFiller::AdvancedEventBits &rhs) const {
    return (getBits()|rhs.getBits());
}

const size_t EventBitsFiller::AdvancedEventBits::lowestZero() const { 
    for(size_t i = 0;i<MAX;++i) if(passedCuts_[i]==false) return i;
    return MAX;
}

bool EventBitsFiller::AdvancedEventBits::hasSameBits(const EventBitsFiller::AdvancedEventBits& rhs) const {
    return ( getBits() == rhs.getBits() );
}

bool EventBitsFiller::AdvancedEventBits::operator==(const EventBitsFiller::AdvancedEventBits& rhs) const {
    return ( run()==rhs.run() && lumi()==rhs.lumi() && event()==rhs.event());
}

bool EventBitsFiller::AdvancedEventBits::operator>(const EventBitsFiller::AdvancedEventBits& rhs) const {
    return ( lowestZero() > rhs.lowestZero() );
}

bool EventBitsFiller::AdvancedEventBits::operator<(const EventBitsFiller::AdvancedEventBits& rhs) const {
    return ( lowestZero() < rhs.lowestZero() );
}

bool EventBitsFiller::AdvancedEventBits::operator>=(const EventBitsFiller::AdvancedEventBits& rhs) const {
    return ( lowestZero() >= rhs.lowestZero() );
}

bool EventBitsFiller::AdvancedEventBits::operator<=(const EventBitsFiller::AdvancedEventBits& rhs) const {
    return ( lowestZero() <= rhs.lowestZero() );
}



EventBitsFiller::EventBitsFiller(const edm::ParameterSet &p, const TFileDirectory& fs, const std::string &sn) : 
        eventSummaryPopulated_(false), 
        histParams_(p), fd_(fs), sample_(sn) { 

    cutMasks_.push_back(1);
    for(size_t i=1;i<MAX;++i) {
        cutMasks_.push_back( cutMasks_[i-1]|(1<<i) );
    }

    std::vector<std::string> plots = histParams_.getParameterNamesForType<edm::ParameterSet>();
    std::string var;
    for(size_t i=0;i<plots.size();i++) {
        var = histParams_.getParameter<edm::ParameterSet>(plots[i]).getParameter<std::string>("variable");
        stringFunctions_.push_back( StringObjectFunction<reco::SkimEvent>( var ) );
    }
}

void EventBitsFiller::setCutLabels(const std::string &str, const std::vector<std::string> &l) {
    
    for(size_t i=0;i<l.size();++i) {
        cutLabels_[str].push_back(l[i]);
    }

}

void EventBitsFiller::fillHypoVariables(EventBitsFiller::AdvancedEventBits &hypo, const reco::SkimEvent &se) {
    
    hypo.clearEventVariables();
    for(size_t i=0;i<stringFunctions_.size();i++) {
        hypo.addEventVariable( stringFunctions_[i](se) );
    }

}

void EventBitsFiller::operator()(edm::EventBase const * evt,const std::string &str, const size_t &inst, const size_t &cut, const reco::SkimEvent &se ) {

    setEventSummaryPopulated(false);

    EventBitsFiller::AdvancedEventBits temp(
        evt->eventAuxiliary().run(),
        evt->eventAuxiliary().luminosityBlock(),
        evt->eventAuxiliary().event(),
        inst
    );

    HypoList::iterator thisIt = find_if( hypoSummary_[str].begin(), hypoSummary_[str].end(), findHypo(temp) );
    if(thisIt == hypoSummary_[str].end()) {
        hypoSummary_[str].push_back(temp);
        thisIt = hypoSummary_[str].end()-1;
        fillHypoVariables(*thisIt,se);
    }
    thisIt->turnOn(cut);
}

void EventBitsFiller::printHypoSummary() {

    for(HypoSummary::iterator sumIt=hypoSummary_.begin();sumIt!=hypoSummary_.end();++sumIt) {
        std::cout << "Hypotheses Breakdown for " << sumIt->first << std::endl;
        std::vector<int> perCut(numCuts_[sumIt->first],0);
        for(HypoList::iterator evtIt=sumIt->second.begin(); evtIt!=sumIt->second.end();++evtIt) {
            for(size_t i=0;i<numCuts_[sumIt->first];++i) {
                if( ((*evtIt)&cutMasks_[i]).hasSameBits(cutMasks_[i]) ) perCut[i]++;
            }
        }
        for(size_t i=0;i<numCuts_[sumIt->first];++i) {
            std::cout << std::setw(3) << i << std::setw(7) << perCut[i] << std::endl;
        }
    }
    
}

void EventBitsFiller::populateEventSummary() {

    eventSummary_.clear();
    for(HypoSummary::iterator sumIt=hypoSummary_.begin();sumIt!=hypoSummary_.end();++sumIt) {
        for(HypoList::iterator evtIt=sumIt->second.begin(); evtIt!=sumIt->second.end();++evtIt) {
            HypoList::iterator curPos = find_if(eventSummary_[sumIt->first].begin(),eventSummary_[sumIt->first].end(),findEvent(*evtIt));
            if( curPos == eventSummary_[sumIt->first].end() ) {
                eventSummary_[sumIt->first].push_back(*evtIt);
            } else {
                if( (*evtIt) > (*curPos) ) (*curPos) = (*evtIt);
            }
        }
    }
    setEventSummaryPopulated();
}

const std::vector<int> EventBitsFiller::getEventYieldVector(const std::string& str) {
    if( !isEventSummaryPopulated() ) populateEventSummary();
    std::vector<int> perCut(numCuts_[str],0);
    for(HypoList::iterator evtIt=eventSummary_[str].begin(); evtIt!=eventSummary_[str].end();++evtIt) {
        for(size_t i=0;i<numCuts_[str];++i) {
            if( ((*evtIt)&cutMasks_[i]).hasSameBits(cutMasks_[i]) ) perCut[i]++;
        }
    }
    return perCut;
}

void EventBitsFiller::printEventSummary() {

    if( !isEventSummaryPopulated() ) populateEventSummary();

    for(HypoSummary::iterator sumIt=eventSummary_.begin();sumIt!=eventSummary_.end();++sumIt) {
        std::cout << "Event Breakdown for " << sumIt->first << std::endl;
        const std::vector<int> &perCut = getEventYieldVector(sumIt->first);
        for(size_t i=0;i<perCut.size();++i) {
            std::cout << std::setw(3) << i << std::setw(7) << perCut[i] << std::endl;
        }
    }
    
}

void EventBitsFiller::printHypoList(const std::string &str, const size_t &cut) {

    std::cout << "Hypo List for cut " << cut << " and hypothesis " << str << std::endl;
    for(HypoList::iterator evtIt=hypoSummary_[str].begin(); evtIt!=hypoSummary_[str].end();++evtIt) {
        if( ((*evtIt)&cutMasks_[cut]).hasSameBits(cutMasks_[cut]) ) 
            std::cout << "Hypo: " << std::setw(10) << str << ", Cut: " << std::setw(4) << cut << ", Event: " << *evtIt << std::endl;
    }

}

void EventBitsFiller::printEventList(const std::string &str, const size_t &cut) {

    if( !isEventSummaryPopulated() ) populateEventSummary();

    std::cout << "Event List for cut " << cut << " and hypothesis " << str << std::endl;
    for(HypoList::iterator evtIt=eventSummary_[str].begin(); evtIt!=eventSummary_[str].end();++evtIt) {
        if( ((*evtIt)&cutMasks_[cut]).hasSameBits(cutMasks_[cut]) ) 
            std::cout << "Hypo: " << std::setw(10) << str << ", Cut: " << std::setw(4) << cut << ", Event: " << *evtIt << std::endl;
    }

}

void EventBitsFiller::printFuckingEverything() {

    if( !isEventSummaryPopulated() ) populateEventSummary();
    printEventSummary();
    for(HypoSummary::iterator sumIt=eventSummary_.begin();sumIt!=eventSummary_.end();++sumIt) {
        for(size_t i=0;i<numCuts_[sumIt->first];++i) {
            printEventList(sumIt->first,i);
        }
    }

}

void EventBitsFiller::writeYieldHist(const std::string &str) {

    if( !isEventSummaryPopulated() ) populateEventSummary();
    TH1F *h;

    TFileDirectory local = fd_.mkdir("yields/"+str);

    h = local.make<TH1F>(sample_.c_str(),sample_.c_str(),numCuts_[str],0,numCuts_[str]);
    h->Sumw2();

    const std::vector<int> &perCut = getEventYieldVector(str);
    double total = 0;
    for(size_t i=0;i<perCut.size();++i) {
        h->SetBinContent(i+1,perCut[i]);
        h->SetBinError(i+1,TMath::Sqrt(perCut[i]));
        h->GetXaxis()->SetBinLabel(i+1,cutLabels_[str][i].c_str());
        h->LabelsOption("v");
        total+=perCut[i];
    }
    h->SetEntries(total);
    
}

void EventBitsFiller::writeNMinus1Plots(const std::string &str) {

    if( !isEventSummaryPopulated() ) populateEventSummary();


    std::vector<std::string> plots = histParams_.getParameterNamesForType<edm::ParameterSet>();

    std::vector<std::vector<TH1F*> > allHists(plots.size());
    edm::ParameterSet plot;

    std::ostringstream os;
    for(size_t i=0;i<plots.size();++i) {
        plot = histParams_.getParameter<edm::ParameterSet>(plots[i]);
        for(size_t cut=0;cut<numCuts_[str];++cut) {

            os << "nminus1/" << str << "/" << plots[i] << "/" << std::setw(2) << std::setfill('0') << cut;
            TFileDirectory local = fd_.mkdir( os.str() );
            allHists[i].push_back( local.make<TH1F>( 
                sample_.c_str(), 
                sample_.c_str(),
                plot.getParameter<uint>("nbins"),
                plot.getParameter<double>("low"),
                plot.getParameter<double>("high")
            ));
            allHists[i].back()->GetXaxis()->SetTitle( plot.getParameter<std::string>("xtitle").c_str() );
            allHists[i].back()->Sumw2();
            os.clear(); os.str("");
        }
    }

    for(HypoList::iterator evtIt=eventSummary_[str].begin(); evtIt!=eventSummary_[str].end();++evtIt) {
        for(size_t i=0;i<plots.size();++i) {
            for(size_t cut=0;cut<numCuts_[str];++cut) {
                if( ((*evtIt)|(EventBitsFiller::AdvancedEventBits( bits((1<<cut)) )) ).lowestZero() == numCuts_[str]) {
                    allHists[i][cut]->Fill( evtIt->getEventVariable(i) );
                }
            }
        }
    }

}

void EventBitsFiller::writeByCutPlots(const std::string &str) {

    if( !isEventSummaryPopulated() ) populateEventSummary();


    std::vector<std::string> plots = histParams_.getParameterNamesForType<edm::ParameterSet>();

    std::vector<std::vector<TH1F*> > allHists(plots.size());
    edm::ParameterSet plot;

    std::ostringstream os;
    for(size_t i=0;i<plots.size();++i) {
        plot = histParams_.getParameter<edm::ParameterSet>(plots[i]);
        for(size_t cut=0;cut<numCuts_[str];++cut) {

            os << "bycut/" << str << "/" << plots[i] << "/" << std::setw(2) << std::setfill('0') << cut;
            TFileDirectory local = fd_.mkdir( os.str() );
            allHists[i].push_back( local.make<TH1F>( 
                sample_.c_str(), 
                sample_.c_str(),
                plot.getParameter<uint>("nbins"),
                plot.getParameter<double>("low"),
                plot.getParameter<double>("high")
            ));
            allHists[i].back()->GetXaxis()->SetTitle( plot.getParameter<std::string>("xtitle").c_str() );
            allHists[i].back()->Sumw2();
            os.clear(); os.str("");
        }
    }

    for(HypoList::iterator evtIt=eventSummary_[str].begin(); evtIt!=eventSummary_[str].end();++evtIt) {
        for(size_t i=0;i<plots.size();++i) {
            for(size_t cut=0;cut<numCuts_[str];++cut) {
                if( ((*evtIt)&cutMasks_[cut]).hasSameBits(cutMasks_[cut]) ) {
                    allHists[i][cut]->Fill( evtIt->getEventVariable(i) );
                }
            }
        }
    }

}

void EventBitsFiller::writeAllNMinus1Plots() {

    // once for each hypothesis
    if( !isEventSummaryPopulated() ) populateEventSummary();
    for(HypoSummary::iterator sumIt=eventSummary_.begin();sumIt!=eventSummary_.end();++sumIt)  {
        writeNMinus1Plots(sumIt->first);
    }

}

void EventBitsFiller::writeAllByCutPlots() {

    // once for each hypothesis
    if( !isEventSummaryPopulated() ) populateEventSummary();
    for(HypoSummary::iterator sumIt=eventSummary_.begin();sumIt!=eventSummary_.end();++sumIt)  {
        writeByCutPlots(sumIt->first);
    }

}

void EventBitsFiller::writeAllYieldHists() {

    // once for each hypothesis
    if( !isEventSummaryPopulated() ) populateEventSummary();
    for(HypoSummary::iterator sumIt=eventSummary_.begin();sumIt!=eventSummary_.end();++sumIt) 
        writeYieldHist(sumIt->first);

}

