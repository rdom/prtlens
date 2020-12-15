#include "PrtDetectorConstructionMessenger.h"
#include "PrtPrimaryGeneratorAction.h"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"

#include "PrtManager.h"

PrtDetectorConstructionMessenger::PrtDetectorConstructionMessenger
(PrtDetectorConstruction *PrtGeom): G4UImessenger(),
    fPrtGeom(PrtGeom)
{
  fGeomDir = new G4UIdirectory("/prt");
  fGeomDir->SetGuidance("Geometry control");

  fAngleCmd = new G4UIcmdWithADoubleAndUnit("/prt/angle",this);
  fAngleCmd->SetGuidance("Rotation angle of the lens");
  fAngleCmd->SetParameterName("angle",true);
  fAngleCmd->SetRange("angle >= -180 && angle <= 180");
  fAngleCmd->SetDefaultValue(0.);
  fAngleCmd->SetDefaultUnit("deg");

  fLensIdCmd = new G4UIcmdWithAnInteger("/prt/lensId",this);
  fLensIdCmd->SetGuidance("Lens Id");
  fLensIdCmd->SetParameterName("lenseId",true);
  fLensIdCmd->SetRange("lenseId>=0");
  fLensIdCmd->SetDefaultValue(1);
}

PrtDetectorConstructionMessenger::~PrtDetectorConstructionMessenger()
{
  delete fAngleCmd;
  delete fGeomDir;
}

void PrtDetectorConstructionMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{ 
  if( command == fAngleCmd ) {
    G4double angle = fAngleCmd->GetNewDoubleValue(newValue);
    fPrtGeom->SetRotation(angle);
    PrtManager::Instance()->SetAngle(angle);
    std::cout<<"new angle "<<angle/deg <<std::endl;    
  } 
  
  if( command == fLensIdCmd ) {
    G4int id = fLensIdCmd->GetNewIntValue(newValue);
    PrtManager::Instance()->SetLens(id);
    fPrtGeom->SetLens(id);
  }

  if( command == fDetEffType ) {
    G4int id = fDetEffType->GetNewIntValue(newValue);
    //PrtManager::Instance()->SetQuantumEfficiency(id);
    //PrtManager::Instance()->SetDetEffType();
  }
}
