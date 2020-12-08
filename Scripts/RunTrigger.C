#include <TROOT.h>
#include "TriggerEstimate.C"
#include <TString.h>

void RunTrigger(){
   //TString era[7]={"B","C","D","E","F","G","H"};
   //for(unsigned i=3; i<7; ++i){
      //TriggerEstimate T(era[i]);
      TriggerEstimate T("B");    
      T.Loop();
  //}
}
