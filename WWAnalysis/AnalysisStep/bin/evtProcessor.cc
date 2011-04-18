#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "PhysicsTools/SelectorUtils/interface/RunLumiSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

//#if !defined(__CINT__) || !defined(__MAKECINT__)
//Headers for the data items
//#include <WWAnalysis/DataFormats/interface/SkimEvent.h>
//#endif

#include "WWAnalysis/DataFormats/interface/SkimEvent.h"

#include <TSystem.h>
#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "Math/VectorUtil.h"

#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>
using namespace std;
const size_t MAX = 30;

class EventFiller {
    public:

        typedef std::bitset<MAX> bits;

        class MyEventStruct {
            public: 
                MyEventStruct(const unsigned int &r, const unsigned int &l, const unsigned int &e) :
                    run_(r), lumi_(l), evt_(e), passedCuts_(0) {}

                void turnOn(const size_t &bit) { passedCuts_[bit] = true; }
                const unsigned int run() const {return run_;}
                const unsigned int event() const {return evt_;}
                const unsigned int lumi() const {return lumi_;}

                friend ostream& operator<<(ostream& out, const MyEventStruct &evt) {
                    out << setw(9) << evt.run() << setw(15) << evt.lumi() << setw(20) << evt.event();
                    return out;
                }
                bits operator&(const bits &rhs) const {
                    return ( passedCuts_&rhs );
                }
                const size_t lowestZero() const { 
                    for(size_t i = 0;i<MAX;++i) 
                        if (passedCuts_[i] == false) return i;
                    return MAX;
                }
                bool operator==(const MyEventStruct& rhs) const {
                    return ( run_==rhs.run() && lumi_==rhs.lumi() && evt_==rhs.event() );
                }
                bool operator>(const MyEventStruct& rhs) const {
                    return ( lowestZero() > rhs.lowestZero() );
                }
                bool operator<(const MyEventStruct& rhs) const {
                    return ( lowestZero() < rhs.lowestZero() );
                }
                bool operator>=(const MyEventStruct& rhs) const {
                    return ( lowestZero() >= rhs.lowestZero() );
                }
                bool operator<=(const MyEventStruct& rhs) const {
                    return ( lowestZero() <= rhs.lowestZero() );
                }
            private:
                unsigned int run_;
                unsigned int lumi_;
                unsigned int evt_;
                bits passedCuts_;
        };

        class MyHypoStruct : public MyEventStruct {
            public:
                MyHypoStruct(const unsigned int &r, const unsigned int &l, const unsigned int &e, const size_t &i) :
                    MyEventStruct(r,l,e), instance_(i) {}

                const unsigned int instance() const {return instance_;}
                friend ostream& operator<<(ostream& out, const MyHypoStruct &evt) {
                    out << setw(9) << evt.run() << setw(15) << evt.lumi() << setw(20) << evt.event() << setw(5) << evt.instance();
                    return out;
                }
                bool operator==(const MyHypoStruct& rhs) const {
                    return ( run()==rhs.run() && lumi()==rhs.lumi() && event()==rhs.event() && instance()==rhs.instance());
                }
            private:
                size_t instance_;
        };

        typedef std::vector<MyEventStruct> EventList;
        typedef std::vector<MyHypoStruct> HypoList;
        typedef std::map<std::string,EventList> EventSummary;
        typedef std::map<std::string,HypoList> HypoSummary;

        EventFiller() : eventSummaryPopulated_(false) { 
            cutMasks_.push_back(bits(1));
            for(size_t i=1;i<MAX;++i) cutMasks_.push_back( cutMasks_[i-1]|(bits(1<<i)) );
        }
        ~EventFiller() { }

        void operator()(edm::EventBase const * evt,const std::string &str, const size_t &inst, const size_t &cut);
        void printHypoSummary();
        void printEventSummary();
        const bool isEventSummaryPopulated() const {return eventSummaryPopulated_;}
        void populateEventSummary();
        HypoList::iterator findInEventSummary(const std::string &str, const HypoList::iterator &hyp);
        void printHypoList(const std::string &str,const size_t &cut);
        void printEventList(const std::string &str,const size_t &cut);
        void printFuckingEverything();

    private:

        void setEventSummaryPopulated(const bool &b=true) { eventSummaryPopulated_ = b; }

        HypoSummary hypoSummary_;
        HypoSummary eventSummary_;
        bool eventSummaryPopulated_;
        vector<bits> cutMasks_;

};

void EventFiller::operator()(edm::EventBase const * evt,const std::string &str, const size_t &inst, const size_t &cut) {

    MyHypoStruct temp(
        evt->eventAuxiliary().run(),
        evt->eventAuxiliary().luminosityBlock(),
        evt->eventAuxiliary().event(),
        inst
    );

    HypoList::iterator thisIt = find(hypoSummary_[str].begin(),hypoSummary_[str].end(),temp);
    if(thisIt == hypoSummary_[str].end()) {
        hypoSummary_[str].push_back(temp);
        thisIt = hypoSummary_[str].end()-1;
    }
    thisIt->turnOn(cut);
}

void EventFiller::printHypoSummary() {

    for(HypoSummary::iterator sumIt=hypoSummary_.begin();sumIt!=hypoSummary_.end();++sumIt) {
        std::cout << "Hypotheses Breakdown for " << sumIt->first << std::endl;
        std::vector<int> perCut(MAX,0);
        for(HypoList::iterator evtIt=sumIt->second.begin(); evtIt!=sumIt->second.end();++evtIt) {
            for(size_t i=0;i<MAX;++i) {
                if( ((*evtIt)&cutMasks_[i]) == cutMasks_[i] ) perCut[i]++;
            }
        }
        for(size_t i=0;i<MAX;++i) {
            std::cout << setw(3) << i << setw(7) << perCut[i] << std::endl;
        }
    }
    
}

EventFiller::HypoList::iterator EventFiller::findInEventSummary(const std::string &str, const HypoList::iterator &hyp) {
    for(HypoList::iterator evtIt=eventSummary_[str].begin();evtIt!=eventSummary_[str].end();++evtIt) {
        if( (*hyp).MyEventStruct::operator==(*evtIt) ) return evtIt;
    }
    return eventSummary_[str].end();
}

void EventFiller::populateEventSummary() {

    eventSummary_.clear();
    for(HypoSummary::iterator sumIt=hypoSummary_.begin();sumIt!=hypoSummary_.end();++sumIt) {
        for(HypoList::iterator evtIt=sumIt->second.begin(); evtIt!=sumIt->second.end();++evtIt) {
            HypoList::iterator curPos = findInEventSummary(sumIt->first,evtIt);
            if( curPos == eventSummary_[sumIt->first].end() ) {
                eventSummary_[sumIt->first].push_back(*evtIt);
            } else {
                if( (*evtIt) > (*curPos) ) (*curPos) = *evtIt;
            }
        }
    }

    setEventSummaryPopulated();
}

void EventFiller::printEventSummary() {

    if( !isEventSummaryPopulated() ) populateEventSummary();

    for(HypoSummary::iterator sumIt=eventSummary_.begin();sumIt!=eventSummary_.end();++sumIt) {
        std::cout << "Event Breakdown for " << sumIt->first << std::endl;
        std::vector<int> perCut(MAX,0);
        for(HypoList::iterator evtIt=sumIt->second.begin(); evtIt!=sumIt->second.end();++evtIt) {
            for(size_t i=0;i<MAX;++i) {
                if( ((*evtIt)&cutMasks_[i]) == cutMasks_[i] ) perCut[i]++;
            }
        }
        for(size_t i=0;i<MAX;++i) {
            std::cout << setw(3) << i << setw(7) << perCut[i] << std::endl;
        }
    }
    
}

void EventFiller::printHypoList(const std::string &str, const size_t &cut) {

//     std::cout << "Hypo List for cut " << cut << " and hypothesis " << str << std::endl;
    for(HypoList::iterator evtIt=hypoSummary_[str].begin(); evtIt!=hypoSummary_[str].end();++evtIt) {
        if( ((*evtIt)&cutMasks_[cut]) == cutMasks_[cut] ) 
            std::cout << "Hypo: " << setw(10) << str << ", Cut: " << setw(4) << cut << ", Event: " << *evtIt << std::endl;
    }

}

void EventFiller::printEventList(const std::string &str, const size_t &cut) {

    if( !isEventSummaryPopulated() ) populateEventSummary();

//     std::cout << "Event List for cut " << cut << " and hypothesis " << str << std::endl;
    for(HypoList::iterator evtIt=eventSummary_[str].begin(); evtIt!=eventSummary_[str].end();++evtIt) {
        if( ((*evtIt)&cutMasks_[cut]) == cutMasks_[cut] ) 
            std::cout << "Hypo: " << setw(10) << str << ", Cut: " << setw(4) << cut << ", Event: " << *evtIt << std::endl;
    }

}

void EventFiller::printFuckingEverything() {

    printEventSummary();
    for(HypoSummary::iterator sumIt=eventSummary_.begin();sumIt!=eventSummary_.end();++sumIt) {
        for(size_t i=0;i<MAX;++i) {
            printEventList(sumIt->first,i);
        }
    }

}


//cut bits
enum cutBits {
    SKIMMED = 0,
    FIDUCIAL,
    ID,
    ISO,
    CONVERSION,
    IP,
    LEPTONVETO,
    MET,
    MLL,
    MZ,
    PROJMET,
    JETVETO,
    SOFTMU,
    BJETS,
    MLL2,
    PTMAX,
    PTMIN,
    DPHI
};


int main(int argc,char* argv[]) {

    gSystem->Load( "libFWCoreFWLite" );
    AutoLibraryLoader::enable();

    if ( argc < 2 ) {
        std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
        return 0;
    }

    // Get the python configuration
    PythonProcessDesc builder(argv[1]);
    edm::ParameterSet const& allPars = 
        builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("yieldParams");

    // output name, to be handled here or by the Event Describer thingy
    string outputFileName(  allPars.getParameter<string>("fileOutName") );

    // if allPars has the right parameter, it will filter events by lumi
    RunLumiSelector myLumiSel(allPars);

    // selection variables
    double etaMu = 0;
    double etaEl = 0;
    double ptMin = 0;
    double ptMax = 0;
    double d0 = 0;
    double dZ = 0;
    double isoEl = 0;
    double isoMu = 0;
    double met = 0;
    double mll = 0;
    double mZ = 0;
    double pMet = 0;
    double jetPt = 0;
    double jetEta = 0;
    bool useJEC = false;
    int nCentralJet = 0;
    int nSoftMu = 0;
    int nExtraLep = 0;
    double bValue = 0;
    int nBtagJets = 0;
    double mllMaxFinal = 9999999;
    double ptMinFinal  = 0;
    double ptMaxFinal  = 0;
    double deltaPhiLL  = 0; 

    EventFiller eventFiller;

    // vector of strings, each name is a SkimEvent stored in the file
    edm::ParameterSet selectionParams = allPars.getParameter<edm::ParameterSet>("selectionParams");
    vector<string> hypoTypes = selectionParams.getParameterNamesForType<edm::ParameterSet>();
    vector<fwlite::Handle<std::vector<reco::SkimEvent> >* > hypoHandles;
    for(vector<string>::const_iterator itHypo=hypoTypes.begin();itHypo!=hypoTypes.end();++itHypo) {
        hypoHandles.push_back( new fwlite::Handle<std::vector<reco::SkimEvent> >() );
    }

    // vector of strings, each describes an input file
    edm::ParameterSet sampleInputParams = allPars.getParameter<edm::ParameterSet>("inputParams");
    vector<string> sampleInputs = sampleInputParams.getParameterNamesForType<edm::ParameterSet>();


    //loop on datasets
    for(vector<string>::const_iterator itSample=sampleInputs.begin();itSample!=sampleInputs.end();++itSample) {

        // parameter set describing the input file
        edm::ParameterSet input = sampleInputParams.getParameter<edm::ParameterSet>(*itSample);
        fwlite::ChainEvent ev(input.getParameter<vector<string> >("files"));

        // variable just for letting user know wtf is going on
        int evtCount = 0;
        int maxCount = input.getParameter<int>("maxEvents");
        if( maxCount == -1 ) maxCount = numeric_limits<int>::max();

        //loop on events                                     // and fail if not part of lumi list
        if(ev.size()) for( ev.toBegin(); ! ev.atEnd(); ++ev) if(myLumiSel(ev)) {

            //stop if the user is tired o this stuff
            if(evtCount > maxCount) break;
            // let the user know, ya know?
            if( evtCount++%1000 == 0 ) cerr << "Processing " << *itSample << ": " << setw(10) << evtCount << endl;

            //loop on hypothesis
            for(size_t hypoI=0;hypoI<hypoHandles.size();++hypoI) {

                // get the skim event from the file
                hypoHandles[hypoI]->getByLabel(ev,hypoTypes[hypoI].c_str());

                // update the cut variables
                edm::ParameterSet thisSelection = selectionParams.getParameter<edm::ParameterSet>(hypoTypes[hypoI]);
                etaMu = thisSelection.getParameter<double>("etaMu");
                etaEl = thisSelection.getParameter<double>("etaEl");
                ptMax = thisSelection.getParameter<double>("ptMax");
                ptMin = thisSelection.getParameter<double>("ptMin");
                d0 = thisSelection.getParameter<double>("d0");
                dZ = thisSelection.getParameter<double>("dZ");
                isoEl = thisSelection.getParameter<double>("isoEl");
                isoMu = thisSelection.getParameter<double>("isoMu");
                met = thisSelection.getParameter<double>("met");
                mll = thisSelection.getParameter<double>("mll");
                mZ = thisSelection.getParameter<double>("mZ");
                pMet = thisSelection.getParameter<double>("pMet");
                jetPt = thisSelection.getParameter<double>("jetPt");
                jetEta = thisSelection.getParameter<double>("jetEta");
                useJEC = thisSelection.getParameter<bool>("useJEC");
                nCentralJet = thisSelection.getParameter<int>("nCentralJet");
                nSoftMu = thisSelection.getParameter<int>("nSoftMu");
                nExtraLep = thisSelection.getParameter<int>("nExtraLep");
                bValue = thisSelection.getParameter<double>("bValue");
                nBtagJets = thisSelection.getParameter<int>("nBtagJets");
                mllMaxFinal = thisSelection.getParameter<double>("mllMaxFinal");
                ptMinFinal = thisSelection.getParameter<double>("ptMinFinal");
                ptMaxFinal = thisSelection.getParameter<double>("ptMaxFinal");
                deltaPhiLL = thisSelection.getParameter<double>("deltaPhiLL");

                //loop on no of hypos in event
                for(vector<reco::SkimEvent>::const_iterator mySkimEvent = hypoHandles[hypoI]->ptr()->begin();
                        mySkimEvent != hypoHandles[hypoI]->ptr()->end(); mySkimEvent++){

                       size_t instance = mySkimEvent - hypoHandles[hypoI]->ptr()->begin();
                       eventFiller(ev.event(),hypoTypes[hypoI],instance,SKIMMED);

                    if( mySkimEvent->q(0)*mySkimEvent->q(1) < 0 && 
                       !mySkimEvent->isSTA(0) && !mySkimEvent->isSTA(1) && 
                        mySkimEvent->leptEtaCut(etaMu, etaEl) && 
                        mySkimEvent->ptMin() > ptMin && 
                        mySkimEvent->ptMax() > ptMax ) {

                        eventFiller(ev.event(),hypoTypes[hypoI],instance,FIDUCIAL);
                    }

                    if( mySkimEvent->passesIDV1(0) && mySkimEvent->passesIDV1(1) ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,ID);
                    }

                    if( !(abs(mySkimEvent->pdgId(0)) == 11 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= isoEl) && 
                        !(abs(mySkimEvent->pdgId(1)) == 11 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= isoEl) && 
                        !(abs(mySkimEvent->pdgId(0)) == 13 && mySkimEvent->allIso(0)/mySkimEvent->pt(0) >= isoMu) && 
                        !(abs(mySkimEvent->pdgId(1)) == 13 && mySkimEvent->allIso(1)/mySkimEvent->pt(1) >= isoMu) ) {

                        eventFiller(ev.event(),hypoTypes[hypoI],instance,ISO);

                    }

                    if( mySkimEvent->passesConversion(0) && mySkimEvent->passesConversion(1) ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,CONVERSION);
                    }

                    if( fabs(mySkimEvent->d0Reco(0)) < d0 && fabs(mySkimEvent->d0Reco(1)) < d0 &&
                        fabs(mySkimEvent->dZReco(0)) < dZ && fabs(mySkimEvent->dZReco(1)) < dZ ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,IP);
                    }

                    if( mySkimEvent->nExtraLep(10) <= nExtraLep ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,LEPTONVETO);
                    }

                    if( mySkimEvent->tcMet() > met) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,MET);
                    }

                    if( mySkimEvent->mll()>mll) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,MLL);
                    }

                    if( fabs(mySkimEvent->mll()-91.1876)>mZ ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,MZ);
                    }

                    if( mySkimEvent->projTcMet() > pMet ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,PROJMET);
                    }

                    if( mySkimEvent->nCentralJets( jetPt, jetEta, useJEC) <= nCentralJet) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,JETVETO);
                    }

                    if( mySkimEvent->nSoftMu(3.) <= nSoftMu ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,SOFTMU);
                    }

                    if( mySkimEvent->bTaggedJetsUnder( jetPt, bValue) <= nBtagJets ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,BJETS);
                    }

                    if(mySkimEvent->mll() < mllMaxFinal ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,MLL2);
                    }

                    if(mySkimEvent->ptMax() > ptMaxFinal ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,PTMAX);
                    }

                    if(mySkimEvent->ptMin() > ptMinFinal ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,PTMIN);
                    }

                    if(mySkimEvent->dPhill() < deltaPhiLL ) {
                        eventFiller(ev.event(),hypoTypes[hypoI],instance,DPHI);
                    }


                } //end of loop over SkimEvent

            } // end loop over hypothesis

        } // end loop over edm::events 

//         double scaleTo1pb1 = input.getParameter<double>("scale");
//         cout << "Counts: " << *itSample << " == " << "all" << endl;
//         eventFiller.printHypoSummary();
//         eventFiller.printEventSummary();
        eventFiller.printFuckingEverything();
//         eventFiller.printHypoList("wwmumu0",ISO);

    } //end loop over input datasets


    return 0;
}


