#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <string>
#include <cstdlib>

void dumpTrees(TString dir = "zz4lTree", const char *cut="", int maxEvents = 200, const char *process=0) {
    TTree *Events = (TTree *) gFile->Get(dir+"/probe_tree");
    Events->SetScanField(0);
    if (dir.Contains("4l") || dir.Contains("4L")) {
        Events->Scan("event:mass:z1mass:z2mass:l1pdgId:l1pt:l1eta:l2pt:l2eta:l3pdgId:l3pt:l3eta:l4pt:l4eta:massErr:melaLD",cut,"",maxEvents);
    //} else if (dir.Contains("anyZllTree")) {
    //    Events->Scan("event:zmass:l1pdgId:l1pt:l1eta:l2pt:l2eta:l1pfIsoComb04EACorr/l1pt:l2pfIsoComb04EACorr/l2pt:l1sip:l2sip",cut,"",maxEvents);
    } else if (dir.Contains("muonTree")) {
        Events->Scan("event:pt:eta:phi:looseID:newID:prlID:mvaISO:sip:pfIsoComb04EACorr:bdtIso:mvaISO:pfIsoChHad04:pfIsoNHad04_NoEA:pfIsoPhot04_NoEA:pfIsoEAtot:rho:rhoAA",cut,"",maxEvents);
    } else if (dir.Contains("electronTree")) {
        Events->Scan("event:pt:eta:phi:looseID:newID:prlID:sip:pfIsoComb04EACorr:bdtIso:bdtID:mvaISO:mvaID:pfIsoChHad04:pfIsoNHad04_NoEA:pfIsoPhot04_NoEA:pfIsoEAtot:rho:rhoAA",cut,"",maxEvents);
    }
}
