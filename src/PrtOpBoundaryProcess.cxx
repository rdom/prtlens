#include "PrtOpBoundaryProcess.h"

#include "PrtManager.h"

PrtOpBoundaryProcess::PrtOpBoundaryProcess()
  : G4OpBoundaryProcess(){
  fLensId = PrtManager::Instance()->GetLens();
}

G4VParticleChange* PrtOpBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep){
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  G4VParticleChange* particleChange = G4OpBoundaryProcess::PostStepDoIt(aTrack, aStep); 
  
  // int parentId = aTrack.GetParentID();
  // std::cout<<"parentId   "<<parentId <<std::endl;
  // if(parentId==1) particleChange->ProposeTrackStatus(fStopAndKill);
  

  // ideal focusing
  if(PrtManager::Instance()->GetLens() == 10){
    G4String ParticleName = aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName();  
    if (ParticleName == "opticalphoton"){
      double endofbar = 1250/2.;
      G4ThreeVector theGlobalPoint1 = pPostStepPoint->GetPosition();
      G4TouchableHistory* touchable = (G4TouchableHistory*)(pPostStepPoint->GetTouchable());
      G4ThreeVector lpoint =  touchable->GetHistory()->GetTransform( 1 ).TransformPoint(theGlobalPoint1);
     
      if(lpoint.getZ() < endofbar+0.0001 && lpoint.getZ() > endofbar-0.0001){
	G4ThreeVector ww  = pPreStepPoint->GetTouchableHandle()->GetHistory()->
	  GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0,0,endofbar));
	if(aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()!="wBar") 
	  particleChange->ProposeTrackStatus(fStopAndKill);
	else{
	  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
	  theNavigator->LocateGlobalPointWithinVolume(ww);
	  aParticleChange.ProposePosition(ww.getX(), ww.getY(),ww.getZ());
	}
	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
    }
  }

  // // bar surface scattering
  // if(1){
  //   G4String ParticleName = aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName();  
  //   if (ParticleName == "opticalphoton" && aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wBar"){
  //     G4double z = 0.5*PrtManager::Instance()->GetRadiatorL();
  //     G4double w = 0.5*PrtManager::Instance()->GetRadiatorW();
  //     G4double h = 0.5*PrtManager::Instance()->GetRadiatorH();

  //     G4ThreeVector theGlobalPoint1 = pPostStepPoint->GetPosition();
  //     G4TouchableHistory* touchable = (G4TouchableHistory*)(pPostStepPoint->GetTouchable());
  //     G4ThreeVector lpoint =  touchable->GetHistory()->GetTransform(1).TransformPoint(theGlobalPoint1);     
  // 	{
  // 	  std::cout<<w<<" lpoint "<<lpoint<<std::endl;
	  
  // 	  if(lpoint.getY() > w-0.01)
  // 	  {
  // 	  std::cout<<h<<" lpoint.getX() "<<lpoint.getX()<<" "<<lpoint.getY()<<std::endl;
	  
  // 	  G4ThreeVector ww  = pPreStepPoint->GetTouchableHandle()->GetHistory()->
  // 	    GetTopTransform().Inverse().TransformPoint(G4ThreeVector(lpoint.getX(),lpoint.getY(),lpoint.getZ()));
	  
  // 	  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  // 	  theNavigator->LocateGlobalPointWithinVolume(ww);
  // 	  aParticleChange.ProposePosition(ww.getX(), ww.getY(),ww.getZ());
  // 	  // aParticleChange.ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz)

  // 	  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  // 	}
  //     }
  //   }
  // }

  if(PrtManager::Instance()->GetRunType() == 1 && pPostStepPoint->GetPosition().z()<pPreStepPoint->GetPosition().z()){    
    particleChange->ProposeTrackStatus(fStopAndKill);
  }
  
  // kill reflections from FP
  if(PrtManager::Instance()->GetRunType() == 0 && pPreStepPoint->GetPhysicalVolume()->GetName()=="wPrizm"){
    auto touchable = (G4TouchableHistory*)(pPreStepPoint->GetTouchable());
    auto pos1 = touchable->GetHistory()->GetTopTransform().TransformPoint(pPreStepPoint->GetPosition());
    auto pos2 = touchable->GetHistory()->GetTopTransform().TransformPoint(pPostStepPoint->GetPosition());
    if(pos2.y()>pos1.y()) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if(PrtManager::Instance()->GetRunType() == 5 &&  pPreStepPoint->GetPhysicalVolume()->GetName()=="wDirc" && pPostStepPoint->GetPhysicalVolume()->GetName()=="wPrizm" && GetStatus() == FresnelRefraction){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }
  
  // kill photons outside bar and prizm
  if(GetStatus() == FresnelRefraction 
     && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    //rd for air gap
    if(PrtManager::Instance()->GetLens()!=4) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens1" 
      || aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens2"
      || aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3") 
     &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    //if(PrtManager::Instance()->GetLens()!=4) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wOpticalGreased"
      || aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wOpticalGrease")
     &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    if(PrtManager::Instance()->GetLens()==2) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  return particleChange;

}
