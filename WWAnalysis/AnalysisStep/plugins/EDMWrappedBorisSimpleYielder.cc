#include "FWCore/Framework/interface/MakerMacros.h"

#include "WWAnalysis/AnalysisStep/interface/BorisSimpleYielder.h"
#include "PhysicsTools/UtilAlgos/interface/EDAnalyzerWrapper.h"
typedef edm::AnalyzerWrapper<BorisSimpleYielder> EDMWrappedBorisSimpleYielder;
DEFINE_FWK_MODULE(EDMWrappedBorisSimpleYielder);
