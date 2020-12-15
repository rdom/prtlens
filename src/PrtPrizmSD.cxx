// -----------------------------------------
// PrtPrizmSD class
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#include "PrtPrizmSD.h"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"
#include <TVector3.h>

#include "PrtEvent.h"

#include "PrtRunAction.h"
#include "PrtManager.h"

PrtPrizmSD::PrtPrizmSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName,
                            G4int nofCells)
  : G4VSensitiveDetector(name), fHitsCollection(NULL){
  collectionName.insert(hitsCollectionName);
}

PrtPrizmSD::~PrtPrizmSD(){ 
}

void PrtPrizmSD::Initialize(G4HCofThisEvent* hce){ 
  // Create hits collection
  fHitsCollection = new PrtPrizmHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
}

G4bool PrtPrizmSD::ProcessHits(G4Step* aStep, G4TouchableHistory* hist){
  G4StepPoint* pPostStepPoint = aStep->GetPostStepPoint(); 
  if (pPostStepPoint->GetStepStatus() == fGeomBoundary){
    PrtPrizmHit* newHit = new PrtPrizmHit();
    
    newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
    newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());

    // // store normal to the closest boundary
    G4Navigator* theNavigator 
      = G4TransportationManager::GetTransportationManager()
      ->GetNavigatorForTracking();

    int nid = 0;
    G4bool valid;
    G4ThreeVector normal = theNavigator->GetLocalExitNormal(&valid);
    if (valid ){
      // G4ThreeVector gnormal = theNavigator->GetLocalToGlobalTransform().TransformAxis(normal);
      if(normal.y()> 0.99) nid = 1; // right
      if(normal.y()<-0.99) nid = 2; // left
      if(normal.x()<-0.99) nid = 3; // bottom
      if(fabs(normal.x()-0.838671)<0.4 ) nid = 4; // top
    }
    newHit->SetNormalId(nid);
    fHitsCollection->insert(newHit);
  }

  return true;
}

void PrtPrizmSD::EndOfEvent(G4HCofThisEvent*){ 
  if ( verboseLevel>1 ) { 
    G4int nofHits = fHitsCollection->entries();
    G4cout << "\n-------->Prizm Hits Collection: in this event they are " << nofHits 
	   << " hits in the tracker chambers: " << G4endl;
    for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}
