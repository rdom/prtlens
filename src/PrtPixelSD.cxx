#include "PrtPixelSD.h"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include <TVector3.h>
#include "G4TransportationManager.hh"

#include "PrtEvent.h"
#include "PrtPrizmHit.h"

#include "PrtRunAction.h"
#include "PrtManager.h"

PrtPixelSD::PrtPixelSD( const G4String& name, 
			const G4String& hitsCollectionName,
			G4int nofCells)
  : G4VSensitiveDetector(name){
  collectionName.insert(hitsCollectionName);
}

PrtPixelSD::~PrtPixelSD(){ 

}

void PrtPixelSD::Initialize(G4HCofThisEvent* hce){

}

G4bool PrtPixelSD::ProcessHits(G4Step* step, G4TouchableHistory* hist){
  
  if(step == 0) return false; 
  G4TouchableHistory* touchable = (G4TouchableHistory*)(step->GetPostStepPoint()->GetTouchable());

  G4Track* track = step->GetTrack();
  G4ThreeVector globalpos = step->GetPostStepPoint()->GetPosition();
  G4ThreeVector localpos = touchable->GetHistory()->GetTopTransform().TransformPoint(globalpos);

  G4ThreeVector globalvec = track->GetVertexMomentumDirection();
  G4ThreeVector localvec = touchable->GetHistory()->GetTopTransform().TransformAxis(globalvec);

  G4ThreeVector g4mom = track->GetVertexMomentumDirection(); // track->GetMomentum(); 
  G4ThreeVector g4pos = track->GetVertexPosition();
 
  TVector3 globalPos(globalpos.x(),globalpos.y(),globalpos.z());
  TVector3 localPos(localpos.x(),localpos.y(),localpos.z());
    
  //focal plane scan
  globalPos = TVector3(globalpos.x(),globalpos.y(),globalpos.z());
  localPos = TVector3(g4pos.x(),g4pos.y(),g4pos.z());
  
  TVector3 vertexPos(g4pos.x(),g4pos.y(),g4pos.z());
  G4ThreeVector lp = touchable->GetHistory()->GetTransform(1).TransformPoint(g4pos); //pos in wDirc
  TVector3 position(lp.x(),lp.y(),lp.z());
    
  PrtHit hit;
  hit.SetPosition(globalPos);
  hit.SetLeadTime(step->GetPreStepPoint()->GetLocalTime());
  double wavelength = 1.2398/(track->GetMomentum().mag()*1E6)*1000;
  hit.SetTotTime(wavelength);

  PrtManager::Instance()->AddHit(hit,localPos,localPos,position,vertexPos);

  return true; 
}

void PrtPixelSD::EndOfEvent(G4HCofThisEvent*){
  memset(fMultHit, 0, sizeof(fMultHit[0][0])*960);
  G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  if(eventNumber%1==0 && PrtManager::Instance()->GetRunType()==0) std::cout<<"Event # "<<eventNumber <<std::endl;
  if(eventNumber%1000==0 && PrtManager::Instance()->GetRunType()!=0) std::cout<<"Event # "<<eventNumber <<std::endl;
  PrtManager::Instance()->Fill();
}

