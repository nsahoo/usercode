#include "FWCore/MessageLogger/interface/AbstractMLscribe.h"
#include "FWCore/MessageLogger/interface/ErrorObj.h"
#include "FWCore/MessageLogger/interface/MessageLoggerQ.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"



class SilentMLscribe : public edm::service::AbstractMLscribe {
  
 public:
  SilentMLscribe() {}
  
  // ---------- member functions ---------------------------
  virtual
    void  runCommand(edm::MessageLoggerQ::OpCode  opcode, void * operand);
  
 private:
  SilentMLscribe(const SilentMLscribe&); // stop default
  
  const SilentMLscribe& operator=(const SilentMLscribe&); // stop default
  
  // ---------- member data --------------------------------
  
};      

