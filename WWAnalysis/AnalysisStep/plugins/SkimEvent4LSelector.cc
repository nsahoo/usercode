#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "WWAnalysis/DataFormats/interface/SkimEvent4L.h"

typedef SingleObjectSelector<
            std::vector<reco::SkimEvent4L>,
            StringCutObjectSelector<reco::SkimEvent4L>,
            std::vector<reco::SkimEvent4L>
        > SkimEvent4LSelector;

DEFINE_FWK_MODULE(SkimEvent4LSelector);
