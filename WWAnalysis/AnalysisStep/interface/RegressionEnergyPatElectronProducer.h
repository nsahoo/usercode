
#ifndef RegressionEnergyPatElectronProducer_h
#define RegressionEnergyPatElectronProducer_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEnergyRegressionEvaluate.h"

class RegressionEnergyPatElectronProducer: public edm::EDProducer 
 {
  public:

    explicit RegressionEnergyPatElectronProducer( const edm::ParameterSet & ) ;
    virtual ~RegressionEnergyPatElectronProducer();
    virtual void produce( edm::Event &, const edm::EventSetup & ) ;

  private:

    edm::InputTag inputPatElectrons ;
    uint32_t energyRegressionType ;
    std::string regressionInputFile;
    bool debug ;
    ElectronEnergyRegressionEvaluate *regressionEvaluator;

 } ;

#endif
