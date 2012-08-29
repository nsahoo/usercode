#include "WWAnalysis/TreeModifiers/interface/BranchAdder.h"
#include <iostream>

int BranchAdder::addBranch(TTree& tree){
  using namespace std;
  
  TObjArray* listBranches = tree.GetListOfBranches(); 
  int nBranches = listBranches->GetEntries();
  for(int i=0; i< nBranches; ++i){
    string nameThisBranch = ((TBranch*) (*listBranches)[i])->GetName();
    if(nameThisBranch == branchName()){
      cout << ">>> ERROR: the branch \"" << branchName() << "\" that you are trying to add is already present in the tree." << endl;
      cout << ">>> The addBranch() operation will be skipped. Please try again with a different input file." << endl << endl;      
      return -1;
    }
  }

  float value;
  TBranch* newBranch = tree.Branch(branchName().c_str(), &value);

  int  nentries = tree.GetEntries();
  initialize(tree);
  for(int k=0; k<nentries; k++){
    calculateValue(tree,k,value);  
    newBranch->Fill(); 
  }

  tree.Write("",TObject::kOverwrite);
    
  return 0;
}
