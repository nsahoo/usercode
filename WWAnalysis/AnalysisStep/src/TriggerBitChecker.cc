#include "WWAnalysis/AnalysisStep/interface/TriggerBitChecker.h"

#include "FWCore/Common/interface/TriggerNames.h"

bool TriggerBitChecker::check(const edm::EventBase &event, const edm::TriggerResults &result) const {
    if (event.id().run() != lastRun_) { syncIndices(event, result); lastRun_ = event.id().run(); }
    for (std::vector<unsigned int>::const_iterator it = indices_.begin(), ed = indices_.end(); it != ed; ++it) {
        if (result.accept(*it)) return true;
    }
    return false;
}

void TriggerBitChecker::syncIndices(const edm::EventBase &event, const edm::TriggerResults &result) const {
    indices_.clear();
    const edm::TriggerNames &names = event.triggerNames(result);
    std::vector<std::string>::const_iterator itp, bgp = paths_.begin(), edp = paths_.end();
    for (size_t i = 0, n = names.size(); i < n; ++i) {
        const std::string &thispath = names.triggerName(i);
        for (itp = bgp; itp != edp; ++itp) {
            if (thispath.find(*itp) == 0) indices_.push_back(i);
        }
    }
}
