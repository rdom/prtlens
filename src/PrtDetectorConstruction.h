// -----------------------------------------
// PrtDetectorConstruction class
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtDetectorConstruction_h
#define PrtDetectorConstruction_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"

#include "PrtDetectorConstructionMessenger.h"

class PrtDetectorConstructionMessenger;

class PrtDetectorConstruction : public G4VUserDetectorConstruction
{ 
public:
  PrtDetectorConstruction();
  virtual ~PrtDetectorConstruction();

public:
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  void DefineMaterials();
  void SetVisualization();
  void SetRotation(G4double angle);
  void SetLens(G4int id);
  void SetQuantumEfficiency(G4int id);
  

private:
  G4LogicalVolume* lExpHall;
  G4LogicalVolume* lLens1;
  G4LogicalVolume* lLens2;
  G4LogicalVolume* lLens3;
  G4LogicalVolume* lTank;
  G4LogicalVolume* lPixel;

  G4Material*        defaultMaterial; // material for bars
  G4Material*        BarMaterial; // material for bars
  G4Material*        OilMaterial;
  G4Material*        MirrorMaterial; // material of mirror
  G4Material*        opticalGreaseMaterial; // material of mirror
  G4Material*        opticalCookieMaterial;
  G4Material*        epotekMaterial;  
  G4Material*        Nlak33aMaterial;
  G4Material*        PbF2Material;  
  G4Material*        frontMaterial;
  
  G4int fGeomId;
  G4int fLensId;
  G4double fTank[3];
  G4double fLens[4];
  G4double fPixelZ;
  G4ThreeVector fCenterShift;
  
  G4double fRotAngle;
  G4RotationMatrix *fPrtRot;
  PrtDetectorConstructionMessenger* fGeomMessenger;
  G4double *fQuantumEfficiency;
};

#endif
