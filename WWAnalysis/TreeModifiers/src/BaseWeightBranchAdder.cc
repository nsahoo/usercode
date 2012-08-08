#include "WWAnalysis/TreeModifiers/interface/BaseWeightBranchAdder.h"


BaseWeightBranchAdder::BaseWeightBranchAdder(const edm::ParameterSet& pset) : BranchAdder(pset){
  std::string datasetConf_     = pset.getParameter<std::string>("datasetConf"); 
  if(datasetConf_ == "7TeV") set7TeVMap();
  else if(datasetConf_ == "8TeV") set8TeVMap();
  else{
    std::cout << "WARNING: you didn't specify the \"datasetConf\" parameter. Using default: 8TeV" 
	      << std::endl;
    set8TeVMap();    
  }


}

float BaseWeightBranchAdder::calculateValue(TTree& tree,int k,  int id,
					    float& value)  {
  using namespace std;

  tree.GetEntry(k);

  map<int,double>::iterator it =  baseWeightMap_.find(id);
  if(it != baseWeightMap_.end())
    value = it->second;
  else{
    cout << "WARNING: baseWeight is not available for id=" << id << " and so it is set to 1." << endl;
  }

  return 0;
}


void BaseWeightBranchAdder::set8TeVMap(){
  using namespace std;
  /* --- To parse Adish vector
    cat fromAdish.txt |awk '{if($1) print $1$3$4$5$6$7$8$9$10}'|awk -F[ '{print $2}'| \
    awk -F] '{print "pair<int,double>("$1","$2");"}' |sed 's/;//g' | awk '{ print "baseWeightMap_.insert("$1");" }';
  */

  baseWeightMap_.insert(pair<int,double>(100,2.24/550862.0));
  baseWeightMap_.insert(pair<int,double>(101,4.47/400973.0));
  baseWeightMap_.insert(pair<int,double>(102,76.91/1489093.0));
  baseWeightMap_.insert(pair<int,double>(103,76.91/1476065.0));
  baseWeightMap_.insert(pair<int,double>(104,76.91/824466.0));
  baseWeightMap_.insert(pair<int,double>(105,176.7/1617951.0));
  baseWeightMap_.insert(pair<int,double>(106,176.7/823911.0));
  baseWeightMap_.insert(pair<int,double>(107,176.7/823922.0));
  baseWeightMap_.insert(pair<int,double>(111,3503700.0/17782545.0));
  baseWeightMap_.insert(pair<int,double>(112,23640.0/13948878.0));
  baseWeightMap_.insert(pair<int,double>(113,3890.0/497658.0));
  baseWeightMap_.insert(pair<int,double>(114,55530.0/493460.0));
  baseWeightMap_.insert(pair<int,double>(115,11770.0/23777.0));
  baseWeightMap_.insert(pair<int,double>(116,1760.0/1925072.0));
  baseWeightMap_.insert(pair<int,double>(117,30000.0/259961.0));
  baseWeightMap_.insert(pair<int,double>(118,117700.0/119974.0));
  baseWeightMap_.insert(pair<int,double>(119,1057.0/2170000.0));
  baseWeightMap_.insert(pair<int,double>(120,5995.0/1753235.0));
  //
  baseWeightMap_.insert(pair<int,double>(1115,22.96*1000.0*9.22e-5/275378.));
  baseWeightMap_.insert(pair<int,double>(1117,22.20*1000.0*1.18e-4/299982.));
  baseWeightMap_.insert(pair<int,double>(1119,21.48*1000.0*1.50e-4/274180.));
  baseWeightMap_.insert(pair<int,double>(1120,21.13*1000.0*1.67e-4/999921.));
  baseWeightMap_.insert(pair<int,double>(1121,20.80*1000.0*1.86e-4/297566.));
  baseWeightMap_.insert(pair<int,double>(1123,20.15*1000.0*2.29e-4/299977.));
  baseWeightMap_.insert(pair<int,double>(1124,19.83*1000.0*2.52e-4/299982.));
  baseWeightMap_.insert(pair<int,double>(1125,19.52*1000.0*2.76e-4/995117.));
  baseWeightMap_.insert(pair<int,double>(1126,19.22*1000.0*3.02e-4/299974.));
  baseWeightMap_.insert(pair<int,double>(1127,18.92*1000.0*3.28e-4/299971.));
  baseWeightMap_.insert(pair<int,double>(1145,14.59*1000.0*8.14e-4/299971.));
  baseWeightMap_.insert(pair<int,double>(1150,13.65*1000.0*8.50e-4/999314.));
  baseWeightMap_.insert(pair<int,double>(1180,8.874*1000.0*6.12e-4/297273.));
  baseWeightMap_.insert(pair<int,double>(1200,7.127*1000.0*2.59e-3/894023.));
  baseWeightMap_.insert(pair<int,double>(1220,6.038*1000.0*2.87e-3/299968.));
  //
  baseWeightMap_.insert(pair<int,double>(2115,1.709*1000.0*9.22e-5/49937.));
  baseWeightMap_.insert(pair<int,double>(2117,1.678*1000.0*1.18e-4/49958.));
  baseWeightMap_.insert(pair<int,double>(2119,1.647*1000.0*1.50e-4/49940.));
  baseWeightMap_.insert(pair<int,double>(2120,1.632*1000.0*1.67e-4/49936.));
  baseWeightMap_.insert(pair<int,double>(2121,1.617*1000.0*1.86e-4/49941.));
  baseWeightMap_.insert(pair<int,double>(2123,1.588*1000.0*2.29e-4/49936.));
  baseWeightMap_.insert(pair<int,double>(2124,1.573*1000.0*2.52e-4/49941.));
  baseWeightMap_.insert(pair<int,double>(2125,1.559*1000.0*2.76e-4/49939.));
  baseWeightMap_.insert(pair<int,double>(2126,1.544*1000.0*3.02e-4/49945.));
  baseWeightMap_.insert(pair<int,double>(2127,1.531*1000.0*3.28e-4/49953.));
  baseWeightMap_.insert(pair<int,double>(2130,1.490*1000.0*4.14e-4/49946.));
  baseWeightMap_.insert(pair<int,double>(2150,1.251*1000.0*8.50e-4/49945.));
  baseWeightMap_.insert(pair<int,double>(2180,0.9820*1000.0*6.12e-4/49951.));
  baseWeightMap_.insert(pair<int,double>(2200,0.8441*1000.0*2.59e-3/49946.));
  baseWeightMap_.insert(pair<int,double>(2220,0.7229*1000.0*2.87e-3/49945.));

}


void BaseWeightBranchAdder::set7TeVMap(){
  using namespace std;
  /* --- To parse Adish vector
    cat fromAdish.txt |awk '{if($1) print $1$3$4$5$6$7$8$9$10}'|awk -F[ '{print $2}'| \
    awk -F] '{print "pair<int,double>("$1","$2");"}' |sed 's/;//g' | awk '{ print "baseWeightMap_.insert("$1");" }';
  */


  baseWeightMap_.insert(pair<int,double>(100,1.74/524301.0));
  baseWeightMap_.insert(pair<int,double>(101,3.48/486993.0));
  baseWeightMap_.insert(pair<int,double>(102,15.34/499924.0));
  baseWeightMap_.insert(pair<int,double>(103,15.34/377430.0));
  baseWeightMap_.insert(pair<int,double>(104,15.34/499026.0));
  baseWeightMap_.insert(pair<int,double>(105,30.67/224340.0));
  baseWeightMap_.insert(pair<int,double>(106,30.67/493618.0));
  baseWeightMap_.insert(pair<int,double>(107,30.67/499927.0));
  baseWeightMap_.insert(pair<int,double>(111,3048000.0/36209629.0));
  baseWeightMap_.insert(pair<int,double>(112,17320.0/10189374.0));
  baseWeightMap_.insert(pair<int,double>(113,3190.0/259971.0));
  baseWeightMap_.insert(pair<int,double>(114,41920.0/3900171.0));
  baseWeightMap_.insert(pair<int,double>(115,7870.0/814390.0));
  baseWeightMap_.insert(pair<int,double>(116,1440.0/137980.0));
  baseWeightMap_.insert(pair<int,double>(117,22650.0/1944826.0));
  baseWeightMap_.insert(pair<int,double>(118,7870.0/323401.0));
  baseWeightMap_.insert(pair<int,double>(119,868.0/997759.0));
  baseWeightMap_.insert(pair<int,double>(120,4880.0/210667.0));
  baseWeightMap_.insert(pair<int,double>(121,0.1525*1000.0*0.3283/410728.));
  baseWeightMap_.insert(pair<int,double>(122,0.0659*1000.0*0.4131/397435.));
  baseWeightMap_.insert(pair<int,double>(123,0.1523*1000.0*0.3413/318885.));
  baseWeightMap_.insert(pair<int,double>(124,0.0664*1000.0*0.3764/478033.));
  baseWeightMap_.insert(pair<int,double>(125,0.1517*1000.0*0.3410/418014.));
  baseWeightMap_.insert(pair<int,double>(126,0.0661*1000.0*0.3765/479652.));
  baseWeightMap_.insert(pair<int,double>(200,18.14*1000.0*9.22e-5/296850.));
  baseWeightMap_.insert(pair<int,double>(201,16.65*1000.0*1.67e-4/238879.));
  baseWeightMap_.insert(pair<int,double>(202,14.16*1000.0*4.14e-4/295356.));
  baseWeightMap_.insert(pair<int,double>(203,12.18*1000.0*7.11e-4/299973.));
  baseWeightMap_.insert(pair<int,double>(204,10.58*1000.0*8.50e-4/285262.));
  baseWeightMap_.insert(pair<int,double>(205,9.202*1000.0*4.25e-4/237830.));
  baseWeightMap_.insert(pair<int,double>(206,7.786*1000.0*2.41e-4/292891.));
  baseWeightMap_.insert(pair<int,double>(207,6.869*1000.0*6.12e-4/285581.));
  baseWeightMap_.insert(pair<int,double>(208,6.017*1000.0*2.12e-3/243498.));
  baseWeightMap_.insert(pair<int,double>(209,5.377*1000.0*2.59e-3/214799.));
  baseWeightMap_.insert(pair<int,double>(210,4.922*1000.0*2.78e-3/199978.));
  baseWeightMap_.insert(pair<int,double>(211,4.528*1000.0*2.87e-3/299075.));
  baseWeightMap_.insert(pair<int,double>(250,1.332*1000.0*9.22e-5/49939.));
  baseWeightMap_.insert(pair<int,double>(251,1.269*1000.0*1.67e-4/49935.));
  baseWeightMap_.insert(pair<int,double>(252,1.154*1000.0*4.14e-4/49948.));
  baseWeightMap_.insert(pair<int,double>(253,1.052*1000.0*7.11e-4/49950.));
  baseWeightMap_.insert(pair<int,double>(254,0.9617*1000.0*8.50e-4/45443.));
  baseWeightMap_.insert(pair<int,double>(255,0.8787*1000.0*4.25e-4/49955.));
  baseWeightMap_.insert(pair<int,double>(256,0.8163*1000.0*2.41e-4/49944.));
  baseWeightMap_.insert(pair<int,double>(257,0.7480*1000.0*6.12e-4/49357.));
  baseWeightMap_.insert(pair<int,double>(258,0.6925*1000.0*2.12e-3/49955.));
  baseWeightMap_.insert(pair<int,double>(259,0.6371*1000.0*2.59e-3/49938.));
  baseWeightMap_.insert(pair<int,double>(260,0.5869*1000.0*2.78e-3/49416.));
  baseWeightMap_.insert(pair<int,double>(261,0.5420*1000.0*2.87e-3/49936.));

}
