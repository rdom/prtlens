#include "PrtSteppingAction.h"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"

#include "G4TrackingManager.hh"
#include "PrtManager.h"

PrtSteppingAction::PrtSteppingAction()
: G4UserSteppingAction()
{ 
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fEventNumber = -1;
}

PrtSteppingAction::~PrtSteppingAction(){ }


void PrtSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* track = step->GetTrack();

  if(track->GetCurrentStepNumber()>50000 || track->GetTrackLength() > 10000) {
    track->SetTrackStatus(fStopAndKill);
  }

  if(step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wPixel") {
    track->SetTrackStatus(fStopAndKill);
  }  
}

