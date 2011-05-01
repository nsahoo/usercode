/** Debug macro to print trigger paths, objects and matches for events
    Must be compiled:
        root.exe -b -l -q file.root printTriggers.cxx+ 

    You can specify the process name, and optionally a trigger path name (if you do, only events passing the path are dumped)

    When printing trigger paths associated to an object, it will put a (!) after objects that are from the last filter of the
    path (if the path is succeeding), and a (-) for intermediate objects (or objects in a path that failed)
*/
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include <TFile.h>
#include <iostream>

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include <algorithm>
#endif


#if !defined(__CINT__) && !defined(__MAKECINT__)
void dumpTOS( const pat::TriggerObjectStandAloneCollection & objs ) {
    for (unsigned i = 0; i < objs.size(); ++i) {
        const pat::TriggerObjectStandAlone &o = objs[i];
        std::cout << "\tTrigger object " << i << ", pt " << o.pt() << ", eta " << o.eta() << ", phi " << o.phi() << std::endl;
        std::cout << "\t   Collection: " << o.collection() << std::endl;
        std::cout << "\t   Type IDs:   ";
        for (unsigned h = 0; h < o.filterIds().size(); ++h) std::cout << " " << o.filterIds()[h] ;
        std::cout << std::endl;
        std::cout << "\t   Filters:    ";
        for (unsigned h = 0; h < o.filterLabels().size(); ++h) std::cout << " " << o.filterLabels()[h];
        std::cout << std::endl;
        std::vector<std::string> pathNamesAll  = o.pathNames(false);
        std::vector<std::string> pathNamesLast = o.pathNames(true);
        std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
            bool isLast = std::count(pathNamesLast.begin(), pathNamesLast.end(), pathNamesAll[h]);
            std::cout << " " << pathNamesAll[h] << (isLast?"(!) ":"(-) ");
        }
        std::cout << std::endl;
    }
}
#endif

void printTriggers(const char *process = "HLT", const char *path="") {
#if !defined(__CINT__) && !defined(__MAKECINT__)
    fwlite::Event ev(gFile);
    fwlite::Handle<edm::TriggerResults> hTriggerResults;
    fwlite::Handle<std::vector<pat::Electron> > hEle;
    fwlite::Handle<std::vector<pat::Muon> > hMu;
    fwlite::Handle<pat::TriggerObjectStandAloneCollection> hHLTObjs;

    int iEvent = 0;
    for (ev.toBegin(); ! ev.atEnd(); ++ev) {
        ++iEvent;
        hTriggerResults.getByLabel(ev,"TriggerResults","",process);
        edm::TriggerNames const&  triggerNames = ev.triggerNames(*hTriggerResults);
        bool pathPass = 0;
        if (path != NULL && strlen(path) != 0) {
            int index = triggerNames.triggerIndex(path);
            if (index == triggerNames.size()) {  std::cerr << "Trigger " << path << " not found." << std::endl; continue; }
            if (!hTriggerResults->accept(index)) continue;
        }
        std::cout << "\n\n\n === EVENT " << iEvent << ": HLT RESULTS ====" << std::endl;
        for (unsigned i = 0; i < triggerNames.size(); ++i) {
            if (strstr(triggerNames.triggerName(i).c_str(), "Mu")) {
                std::cout << "\t" << triggerNames.triggerName(i) << ": " << (hTriggerResults->accept(i) ? "PASS" : "FAIL" ) << std::endl;
            }
        }

        hHLTObjs.getByLabel(ev,"patTrigger");
        if (!hHLTObjs.failedToGet()) {
            std::cout << " === EVENT " << iEvent << ": ALL HLT Objects ====" << std::endl;
            const pat::TriggerObjectStandAloneCollection & objs = *hHLTObjs;
            dumpTOS(objs);
        }

        std::cout << " === EVENT " << iEvent << ": Muons ====" << std::endl;
        hMu.getByLabel(ev,"boostedMuons");
        const std::vector<pat::Muon>     & mus = *hMu;
        for (unsigned i = 0; i < mus.size(); ++i) {
            const pat::Muon &mu = mus[i];
            std::cout << "  mu " << i << ", pt " << mu.pt() << ", eta " << mu.eta() << std::endl;
            dumpTOS(mu.triggerObjectMatches());
        } // muon  

        std::cout << " === EVENT " << iEvent << ": Electrons ====" << std::endl;
        hEle.getByLabel(ev,"boostedElectrons");
        const std::vector<pat::Electron> & els = *hEle;
        for (unsigned i = 0; i < els.size(); ++i) {
            const pat::Electron &el = els[i];
            std::cout << "  el " << i << ", pt " << el.pt() << ", eta " << el.eta() << std::endl;
            dumpTOS(el.triggerObjectMatches());
        } // el
        if (iEvent >= 20) break;
    } // event
#endif
} // macro
