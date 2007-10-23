#include "HiggsAnalysis/HZZ4lAnalyzers/interface/HiggsAnalyzer.h"

using namespace edm;
using namespace std;


HiggsAnalyzer::HiggsAnalyzer(const edm::ParameterSet& iConfig)
{
    LogDebug("RWK") << "HiggsAnalzyer::HiggsAnalyzer()" << "\n";
		
    outputRootFileName  = iConfig.getUntrackedParameter<std::string>("outputRootFileName", "higgsStudy.root");
    hepMCLabel          = iConfig.getUntrackedParameter("hepMCLabel", std::string("source") );

}


HiggsAnalyzer::~HiggsAnalyzer()
{
    LogDebug("RWK") << "HiggsAnalzyer::~HiggsAnalyzer()" << "\n";
    // nothing to do here
}


// ------------ method called once each job just before starting event loop  ------------
void HiggsAnalyzer::beginJob(const edm::EventSetup&)
{
  LogDebug("RWK") << "HiggsAnalzyer::BeginJob()" << "\n";
	
  // create a new root file
  outputRootFile = new TFile(outputRootFileName.c_str(), "RECREATE");
  
  // create a new event class
  higgsEvent = new HiggsEvent();

  // create the tree
  higgsEventTree = new TTree("T", "Tree of Higgs Event data");
  higgsEventTree->Branch("T", "HiggsEvent", &higgsEvent);

  // create event counting histo
  //hNumEvents = new TH1I("hNE", "Number of Events Processed", 3, 0.0, 2.0);
  numEventsProcessed = 0;

}


// ------------ method called to for each event  ------------
void HiggsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  LogDebug("RWK") << "HiggsAnalzyer::analyze()" << "\n";
  
  // increment the event counting histo  
  std::cout << "The event being processed:  " <<  ++numEventsProcessed << std::endl;
  
  // create a new event class
  higgsEvent = new HiggsEvent();
  
  // get the event from the record
  edm::Handle<HepMCProduct> hepMCEvent;
  iEvent.getByLabel(hepMCLabel, hepMCEvent);
  const HepMC::GenEvent * generatedEvent = new HepMC::GenEvent( *(hepMCEvent->GetEvent()) );
	
  
  for(HepMC::GenEvent::particle_const_iterator MCParticle = generatedEvent->particles_begin(); 
      MCParticle != generatedEvent->particles_end(); ++MCParticle)
    {
      if( ( (abs( (*MCParticle)->pdg_id() ) == 11) ||
	    (abs( (*MCParticle)->pdg_id() ) == 13) ||
	    (abs( (*MCParticle)->pdg_id() ) == 15)    )	  
	  && (*MCParticle)->status() == 1)
	{
	  
	  HepMC::GenParticle* mc 	= (*MCParticle);
	  HepMC::GenParticle* mcMom 	= GenParticleMother(*MCParticle);
	  HepMC::GenParticle* mcGMom 	= GenParticleGrandMother(*MCParticle);
	  
	  

	  //cout << "=== found lepton type, pt: " << mc->pdg_id() << " , " << mc->momentum().perp() << endl;
	  
	  /*
	  Lepton* e = new Lepton( mc->momentum().px(), mc->momentum().py(), mc->momentum().pz(), 
				  mc->momentum().e(),  mc->generatedMass(), 
				  mc->pdg_id(), mcMom->pdg_id(),  mcGMom->pdg_id(), 
				  mc->barcode(), mcMom->barcode() );
	  */

	  //higgsEvent->e.push_back(e);
	      
	  // if e's mother is a Z, add it the vector (do for e+ to ensure only one copy in the vector)
	  /*
	  if(e->isMomZ() && e->charge() > 0)
	    {
	      ZBoson* Z = new ZBoson( mcMom->momentum().px(), mcMom->momentum().py(), mcMom->momentum().pz(),
				      mcMom->momentum().e(), mcMom->generatedMass(),
				      mcMom->pdg_id(), mcGMom->pdg_id(), GenParticleMother(mcGMom)->pdg_id(),
				      mcMom->barcode(), mcGMom->barcode() );
	      //cout << "=== found Z type, pt: " 
	      //   << mcMom->pdg_id() << " , " 
	      //   << mcMom->momentum().perp() << endl;
	      
	      if( higgsEvent->Z.empty() ) higgsEvent->Z.push_back(Z);		
	      else
		{
		  if(higgsEvent->Z[0]->M_Gen() > Z->M_Gen()) higgsEvent->Z.push_back(Z);
		  else higgsEvent->Z.insert(higgsEvent->Z.begin(), Z);
		}	      
	    }
	  
	  // if e's grand mother is a H, add it the vector 
	  if(e->isGMomHiggs() && higgsEvent->H.empty())
	    {	      
	      HBoson* H = new HBoson( mcGMom->momentum().px(), mcGMom->momentum().py(), mcGMom->momentum().pz(),
				      mcGMom->momentum().e(), mcGMom->generatedMass(),
				      mcGMom->pdg_id(), GenParticleMother(mcGMom)->pdg_id(), 
				      GenParticleGrandMother(mcGMom)->pdg_id(), 
				      mcGMom->barcode(), GenParticleMother(mcGMom)->barcode() );
	      
	      //cout << "=== found H type, pt: " 
	      //   << mcGMom->pdg_id() << " , " 
	      //   << mcGMom->momentum().perp() << endl;

	      higgsEvent->H.push_back(H);			     
	    }
	  */
	}

      if(   (*MCParticle)->pdg_id()  == 25)
	{
	  HepMC::GenParticle* mc 	= (*MCParticle);
	  cout << "=== found higgs type,status, pt: " 
	       << mc->pdg_id() << " , " 
	       << mc->status() << " , " 
	       << mc->momentum().perp() << endl;
	  TLorentzVector vector(mc->momentum().px(),
				mc->momentum().py(),
				mc->momentum().pz(),
				mc->momentum().e());

	  higgsEvent->mc().setHiggs(vector);

	}
      
      if(   abs( (*MCParticle)->pdg_id()  == 23) )
	{
	  HepMC::GenParticle* mc 	= (*MCParticle);
	  cout << "=== found Z type, status, pt: " 
	       << mc->pdg_id() << " , " 
	       << mc->status() << " , " 
	       << mc->momentum().perp() << endl;
	}
    }

   

  
  // Associate the true Z with the their true electrons's
  /*
  for(unsigned short i = 0; i < higgsEvent->Z.size(); i++)
    {
      for(unsigned short j = 0; j < higgsEvent->e.size(); j++)
	{
	  if(higgsEvent->Z[i]->barcode() == higgsEvent->e[j]->mbarcode() && higgsEvent->e[j]->charge() > 0)
	    {
	      higgsEvent->Z[i]->AddElectron(higgsEvent->e[j]);
	    }
	}
      for(unsigned short j = 0; j < higgsEvent->e.size(); j++)
	{
	  if(higgsEvent->Z[i]->barcode() == higgsEvent->e[j]->mbarcode() && higgsEvent->e[j]->charge() < 0)
	    {
	      higgsEvent->Z[i]->AddElectron(higgsEvent->e[j]);
	    }
	}
    }
  */

  // ================== to be adapted. It makes it crash 
  /* 
  // Set the associated e's Zflag (Z1 higher mass, Z1 lower mass)
  if(!higgsEvent->Z.empty())
    {
      higgsEvent->Z[0]->SetZFlag(0,Z1);
      higgsEvent->Z[0]->SetZFlag(1,Z1);
      higgsEvent->Z[1]->SetZFlag(0,Z2);
      higgsEvent->Z[1]->SetZFlag(1,Z2);
    }
  // Set the remaining e's Zflags to noZ's
 
  for(unsigned short i = 0; i < higgsEvent->e.size(); i++)
    {
      if(higgsEvent->e[i]->WhichZ() == notSet) higgsEvent->e[i]->SetZFlag(noZ);
    }
  
	
  
  // Associate the true H with the their true Z's (assumption: only one H per event)
  for(unsigned short j = 0; j < higgsEvent->Z.size(); j++)
    {
      if(higgsEvent->Z[j]->isMomHiggs())
	{
	  higgsEvent->H[0]->AddZ(higgsEvent->Z[j]);
	}
    }
  */
  // =============================================  
  
  // Fill the Tree @ end of loop
  higgsEventTree->Fill();
  

}


// ------------ method called once each job just after ending the event loop  ------------
void HiggsAnalyzer::endJob() 
{
    LogDebug("RWK") << "HiggsAnalzyer::endJob()" << "\n";
	std::cout << "\n\nTotal number of events proecessed:  " << numEventsProcessed << std::endl;

	// write the tree and histo

	higgsEventTree->Write();
	//hNumEvents->Write();
	
	// write and close the output root file
	outputRootFile->Write();
	outputRootFile->Close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Analyzer specific functions

// Return the mother of a generator level particle
HepMC::GenParticle* HiggsAnalyzer::GenParticleMother(HepMC::GenParticle* MCParticle)
{
	HepMC::GenParticle* gp = MCParticle;
	
	// this gets the first genParticle mother that is different.
	while( (*gp->production_vertex()->particles_in_const_begin())->pdg_id() == gp->pdg_id() )
	{
		gp = *gp->production_vertex()->particles_in_const_begin();
	}

	// return the value;
	if( *gp->production_vertex()->particles_in_const_begin() != NULL )
	{
		return *gp->production_vertex()->particles_in_const_begin();
	}
	else
	{
		return new HepMC::GenParticle();
		std::cout << "no mother found" << std::endl;
	}
}


// Return the grandmother of a generator level particle
HepMC::GenParticle* HiggsAnalyzer::GenParticleGrandMother(HepMC::GenParticle* MCParticle)
{
	HepMC::GenParticle* grandmother = GenParticleMother( GenParticleMother(MCParticle) );
	return grandmother;
}


/*
// Set the leptons values
void HiggsAnalyzer::SetLeptonData(HepMC::GenParticle* mcLepton, Lepton *lepton)
{
    lepton->Set4Momentum( mcLepton->momentum().px(), mcLepton->momentum().py(),
                          mcLepton->momentum().pz(), mcLepton->momentum().e() );
    lepton->SetGenMass( mcLepton->generatedMass() );
    lepton->SetPdgId( mcLepton->pdg_id() );
    lepton->SetMother( GenParticleMother(mcLepton)->pdg_id() );
    lepton->SetGrandMother( GenParticleGrandMother(mcLepton)->pdg_id() );
	lepton->SetBarcode( mcLepton->barcode() );
    return;
}
*/
