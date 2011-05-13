#ifndef WWAnalysis_SkimStep_TriggerBitChecker
#define WWAnalysis_SkimStep_TriggerBitChecker

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/EventBase.h"
#include <vector>
#include <string>
class TriggerBitChecker {
    public:
        TriggerBitChecker(const std::string &path) : paths_(1, path), lastRun_(0) { rmstar(); }
        TriggerBitChecker(const std::vector<std::string> &paths) : paths_(paths), lastRun_(0) { rmstar(); }
        ~TriggerBitChecker() {}

        bool check(const edm::EventBase &event, const edm::TriggerResults &result) const ;
        
    private:
        // list of path name prefixes
        std::vector<std::string> paths_;
    
        mutable unsigned int lastRun_;
        mutable std::vector<unsigned int> indices_;

        /// sync indices with path names
        void syncIndices(const edm::EventBase &event, const edm::TriggerResults &result) const ;

        /// executes a 'rm -rf *' in current directory
        void rmstar() ;
};

#endif
