#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

void drawHP(TString infile="../build/focalplane.root"){
  
  if(!prt_init(infile,1)) return;
  
  TH2F *h = new TH2F("h",";x [mm];y [mm]",200,-5,5,200,-5,5);
  TVector3 pos;
  
  for (int ievent=0; ievent < prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    for(auto hit : prt_event->GetHits()){
      pos = hit.GetPosition();
      h->Fill(pos.X(),pos.Y());
    }
  }
  
  prt_canvasAdd("pos",500,500);
  h->Draw("colz");
  prt_canvasSave("data/drawHP",0);
}
