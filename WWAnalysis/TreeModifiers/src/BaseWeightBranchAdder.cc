#include "WWAnalysis/TreeModifiers/interface/BaseWeightBranchAdder.h"


float BaseWeightBranchAdder::calculateValue(TTree& tree,int k, float& value)  {
  using namespace std;

  tree.GetEntry(k);

  //this is just a test. Real calculation has still to be implemented
  value = 1./k;

  return 0;
}
