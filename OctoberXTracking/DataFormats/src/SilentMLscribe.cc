#include "OctoberXTracking/DataFormats/interface/SilentMLscribe.h"

void  SilentMLscribe::runCommand(edm::MessageLoggerQ::OpCode  opcode, void * operand) {
  //even though we don't print, have to clean up memory
  switch (opcode) {
  case edm::MessageLoggerQ::LOG_A_MESSAGE: {
    edm::ErrorObj *  errorobj_p = static_cast<edm::ErrorObj *>(operand);
    //std::cerr<<errorobj_p->xid().severity.getInputStr()<<" "<<errorobj_p->xid().id<<" -----------------------"<<std::endl;
    //std::cerr <<errorobj_p->fullText()<<std::endl;
    delete errorobj_p;
    break;
  }
  case edm::MessageLoggerQ::JOBREPORT:
  case edm::MessageLoggerQ::JOBMODE:
  case edm::MessageLoggerQ::GROUP_STATS:
    std::string* string_p = static_cast<std::string*> (operand);
    delete string_p;
    break;
  default:
    break;
  }
}   

