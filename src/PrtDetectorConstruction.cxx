
#include "PrtDetectorConstruction.h"

#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "PrtManager.h"
#include "PrtTriggerSD.h"
#include "PrtBarSD.h"
#include "PrtPrizmSD.h"
#include "PrtPixelSD.h"

PrtDetectorConstruction::PrtDetectorConstruction()
  : G4VUserDetectorConstruction(){

  fGeomId = PrtManager::Instance()->GetGeometry();
  fLensId = PrtManager::Instance()->GetLens();
  fPixelZ = PrtManager::Instance()->GetBeamZ();
  
  fTank[0] = 200; fTank[1] = 200, fTank[2] = 300;
  
  fLens[0] = fLens[1] = 40; fLens[2]=10;
  if(fLensId == 2){
    fLens[0] = 50; fLens[1] = 175; fLens[2]=14.4;
  }else if(fLensId == 3){
    fLens[0] = 60; fLens[1] = 60; fLens[2]=15;
  }else if(fLensId == 4 ||fLensId == 5){
    fLens[0] = 50; fLens[1] = 50; fLens[2]=5.7;
  }else if(fLensId == 6 || fLensId==7  || fLensId==8){
    fLens[0] = 50; fLens[1] = 175; fLens[2]=12;
  }
     
  fCenterShift =  G4ThreeVector(0., 0., 0.);
    
  fPrtRot = new G4RotationMatrix();
  
  fGeomMessenger = new PrtDetectorConstructionMessenger(this);
  PrtManager::Instance()->AddInfo("Initialization done");
}

PrtDetectorConstruction::~PrtDetectorConstruction(){}

G4VPhysicalVolume* PrtDetectorConstruction::Construct(){
  DefineMaterials();

  auto gExpHall = new G4Box("gExpHall",200,200,500);
  lExpHall = new G4LogicalVolume(gExpHall,defaultMaterial,"lExpHall",0,0,0);    
  auto wExpHall  = new G4PVPlacement(0,G4ThreeVector(),lExpHall,"gExpHall",0,false,0);
    
  auto gTank = new G4Box("gTank",0.5*fTank[0],0.5*fTank[1],0.5*fTank[2]);

  G4Material* tankMaterial = defaultMaterial; // air
  if(fGeomId==2) tankMaterial = OilMaterial;
  lTank = new G4LogicalVolume(gTank,tankMaterial,"lTank",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(),lTank,"wTank",lExpHall,false,0);
    
  // The Lens
  G4Box* gfbox = new G4Box("Fbox",fLens[0]/2.,fLens[1]/2.,fLens[2]/2.);
  
  if(fLensId == 1){ // 2-layer spherical lens
    G4double r1 = 0; // PrtManager::Instance()->GetTest1(); 
    G4double lensrad1 = (r1==0)? 73.58: r1;
    G4double lensMinThikness = 2; 

    G4ThreeVector zTrans1(0, 0, -lensrad1+fLens[2]/2.-lensMinThikness);
    G4Tubs* gftub = new G4Tubs("Ftub",0,fLens[0]/2.,fLens[2]/2.,0.*deg,360.*deg);
    G4Sphere* gsphere = new G4Sphere("Sphere",0,lensrad1,0,360*deg,0,360*deg);
    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Ftub*Sphere", gftub, gsphere,new G4RotationMatrix(),zTrans1); 
    G4SubtractionSolid* gLens2 = new G4SubtractionSolid("Ftub-Sphere", gftub, gsphere,new G4RotationMatrix(),zTrans1);
    lLens1 = new G4LogicalVolume(gLens1,Nlak33aMaterial,"lLens1",0,0,0); //Nlak33aMaterial  
    lLens2 = new G4LogicalVolume(gLens2,BarMaterial,"lLens2",0,0,0);
  }

  if(fLensId == 2){ // 2-layer cylindrical lens
    G4double lensrad = 73.58;
    G4double lensMinThikness = 8;
    G4ThreeVector zTrans(0, 0, -lensrad+fLens[2]/2.-lensMinThikness);
    G4Tubs* gcylinder = new G4Tubs("Cylinder",0,lensrad,fLens[1]/2.+1.1,0,360*deg);
    G4RotationMatrix* xRot = new G4RotationMatrix();
    xRot->rotateX(M_PI/2.*rad);
    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Fbox*Cylinder", gfbox, gcylinder,xRot,zTrans); 
    G4SubtractionSolid* gLens2 = new G4SubtractionSolid("Fbox-Cylinder", gfbox, gcylinder,xRot,zTrans);
    lLens1 = new G4LogicalVolume(gLens1,Nlak33aMaterial,"lLens1",0,0,0); //Nlak33aMaterial  
    lLens2 = new G4LogicalVolume(gLens2,BarMaterial,"lLens2",0,0,0);
  }

  if(fLensId == 334){ // 3-component spherical lens
    G4double lensMinThikness = 2; 
  
    G4double r1 = 0; //PrtManager::Instance()->GetTest1();
    G4double r2 = 0; //PrtManager::Instance()->GetTest2();
  
    G4double lensrad1 = (r1==0)? 47.8: r1;
    G4double lensrad2 = (r2==0)? 29.1: r2;
    
 
    G4ThreeVector zTrans1(0, 0, -lensrad1-fLens[2]/2.+lensrad1-sqrt(lensrad1*lensrad1-fLens[0]/2.*fLens[0]/2.)+lensMinThikness);
    G4ThreeVector zTrans2(0, 0, -lensrad2+fLens[2]/2.-lensMinThikness);
    
    G4Sphere* gsphere1 = new G4Sphere("Sphere1",0,lensrad1,0,360*deg,0,360*deg);
    G4Sphere* gsphere2 = new G4Sphere("Sphere2",0,lensrad2,0,360*deg,0,360*deg);


    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Fbox*Sphere1", gfbox, gsphere1,new G4RotationMatrix(),zTrans1); 
    G4SubtractionSolid* gLenst = new G4SubtractionSolid("Fbox-Sphere1", gfbox, gsphere1, new G4RotationMatrix(),zTrans1);

    G4IntersectionSolid* gLens2 = new G4IntersectionSolid("gLenst*Sphere2", gLenst, gsphere2, new G4RotationMatrix(),zTrans2);
    G4SubtractionSolid* gLens3 = new G4SubtractionSolid("gLenst-Sphere2", gLenst, gsphere2,new G4RotationMatrix(),zTrans2);
    
    lLens1 = new G4LogicalVolume(gLens1,BarMaterial,"lLens1",0,0,0);
    lLens2 = new G4LogicalVolume(gLens2,Nlak33aMaterial,"lLens2",0,0,0);
    lLens3 = new G4LogicalVolume(gLens3,BarMaterial,"lLens3",0,0,0);
  }

  if(fLensId == 3){ // 3-component spherical lens
    G4double lensMinThikness = 2.0; 
  
    G4double r1 = PrtManager::Instance()->GetTest1();
    G4double r2 = PrtManager::Instance()->GetTest2();

    if(PrtManager::Instance()->GetRunType() == 6){ //focal plane scan
       r1 = PrtManager::Instance()->GetTest1();
       r2 = PrtManager::Instance()->GetTest2();
    }
    
    r1 = (r1==0)? 47.80: r1;
    r2 = (r2==0)? 29.12: r2;
    G4double shight = 40;
    G4double bwidth = fLens[2]-lensMinThikness*2;

    G4ThreeVector zTrans1(0, 0, -r1-fLens[2]/2.+r1-sqrt(r1*r1-shight/2.*shight/2.) +lensMinThikness);
    G4ThreeVector zTrans2(0, 0, -r2-fLens[2]/2.+r2-sqrt(r2*r2-shight/2.*shight/2.) +lensMinThikness*2);

    G4Box* gfbox0 = new G4Box("Fbox0",0.5*fLens[0]+1,0.5*fLens[1]+1,0.5*fLens[2]);
    G4Tubs* gftub = new G4Tubs("Ftub",0,0.5*fLens[0],0.5*fLens[2],0.*deg,360.*deg);
    G4Box* gfsbox = new G4Box("Fsbox",0.5*shight,0.5*fLens[1],0.5*fLens[2]);
    G4Tubs* gfstube = new G4Tubs("ftube",0,0.5*shight,0.5*fLens[2],0.*deg,360.*deg);

    G4Sphere* gsphere1 = new G4Sphere("Sphere1",0,r1,0,360*deg,0,360*deg);
    G4Sphere* gsphere2 = new G4Sphere("Sphere2",0,r2,0,360*deg,0,360*deg);

    G4IntersectionSolid* gbbox = new G4IntersectionSolid("bbox", gftub, gfbox0,new G4RotationMatrix(),G4ThreeVector(0,0,lensMinThikness*2)); 
    G4IntersectionSolid* gsbox = new G4IntersectionSolid("sbox", gfstube, gfbox0,new G4RotationMatrix(),G4ThreeVector(0,0,-lensMinThikness*2)); 

    G4UnionSolid* gubox = new G4UnionSolid("unionbox", gbbox, gsbox,new G4RotationMatrix(),G4ThreeVector(0,0,0)); 

    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Lens1", gubox, gsphere1,new G4RotationMatrix(),zTrans1); 
    G4SubtractionSolid*  gLenst = new G4SubtractionSolid("temp", gubox, gsphere1, new G4RotationMatrix(),zTrans1);

    G4IntersectionSolid* gLens2 = new G4IntersectionSolid("Lens2", gLenst, gsphere2, new G4RotationMatrix(),zTrans2);
    G4SubtractionSolid*  gLens3 = new G4SubtractionSolid("Lens3", gLenst, gsphere2,new G4RotationMatrix(),zTrans2);
    
    lLens1 = new G4LogicalVolume(gLens1,BarMaterial,"lLens1",0,0,0);
    if(PrtManager::Instance()->GetTest3()>1) lLens2 = new G4LogicalVolume(gLens2,PbF2Material,"lLens2",0,0,0);
    else lLens2 = new G4LogicalVolume(gLens2,Nlak33aMaterial,"lLens2",0,0,0);
    lLens3 = new G4LogicalVolume(gLens3,BarMaterial,"lLens3",0,0,0);
  }

  if(fLensId == 4){ // Spherical lens with air gap // f =250 , d = , w = 5.7
    G4double r1 = 0; // PrtManager::Instance()->GetTest1(); 
    G4double lensrad1 = (r1==0)? 250: r1;
    G4double lensMinThikness = 2; 

    G4ThreeVector zTrans1(0, 0, -lensrad1+fLens[2]/2.-lensMinThikness);
    G4Tubs* gftub = new G4Tubs("Ftub",0,fLens[0]/2.,fLens[2]/2.,0.*deg,360.*deg);
    G4Sphere* gsphere = new G4Sphere("Sphere",0,lensrad1,0,360*deg,0,360*deg);
    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Ftub*Sphere", gftub, gsphere,new G4RotationMatrix(),zTrans1); 
    lLens1 = new G4LogicalVolume(gLens1,BarMaterial,"lLens1",0,0,0); //Nlak33aMaterial
  }

  if(fLensId == 5){ // Spherical lens with air gap // f =250 , d = , w = 5.7 // black edges
    G4double r1 = 0; // PrtManager::Instance()->GetTest1(); 
    G4double lensrad1 = (r1==0)? 250: r1;
    G4double lensMinThikness = 2; 

    G4ThreeVector zTrans1(0, 0, -lensrad1+fLens[2]/2.-lensMinThikness);
    G4Tubs* gftub = new G4Tubs("Ftub",0,fLens[0]/2.,fLens[2]/2.,0.*deg,360.*deg);
    G4Sphere* gsphere = new G4Sphere("Sphere",0,lensrad1,0,360*deg,0,360*deg);
    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Ftub*Sphere", gftub, gsphere,new G4RotationMatrix(),zTrans1);
    G4SubtractionSolid* gLens2 = new G4SubtractionSolid("Ftub-Sphere", gftub, gsphere,new G4RotationMatrix(),zTrans1);
    lLens1 = new G4LogicalVolume(gLens1,BarMaterial,"lLens1",0,0,0); //Nlak33aMaterial
    lLens2 = new G4LogicalVolume(gLens2,defaultMaterial,"lLens2",0,0,0);
  }

  if(fLensId == 6){ // 3-component cylindrical lens
    G4double lensMinThikness = 2.0; 
  
    G4double r1 = 0; //PrtManager::Instance()->GetTest1();
    G4double r2 = 0; //PrtManager::Instance()->GetTest2();
    
    if(PrtManager::Instance()->GetRunType() == 6){ //focal plane scan
      r1 = PrtManager::Instance()->GetTest1();
      r2 = PrtManager::Instance()->GetTest2();
    }
    
    //RMI lens
    fLens[2]=13.12;
    lensMinThikness = 2.51;
    G4double layer12 = lensMinThikness+ 3.525; //lensMinThikness*2;
    
    // r1 = (r1==0)? 27.45: r1;
    // r2 = (r2==0)? 20.02: r2;

    r1 = (r1==0)? 33: r1;
    r2 = (r2==0)? 24: r2;
    G4double shight = 20;

    G4ThreeVector zTrans1(0, 0, -r1-fLens[2]/2.+r1-sqrt(r1*r1-shight/2.*shight/2.) +lensMinThikness);
    G4ThreeVector zTrans2(0, 0, -r2-fLens[2]/2.+r2-sqrt(r2*r2-shight/2.*shight/2.) +layer12);

    G4Box* gftub = new G4Box("ftub",0.5*fLens[0],0.5*fLens[1],0.5*fLens[2]);
    G4Box* gcbox = new G4Box("cbox",0.5*fLens[0],0.5*fLens[1]+1,0.5*fLens[2]);
    G4ThreeVector tTrans1( 0.5*(fLens[0]+shight),0,-fLens[2]+layer12);
    G4ThreeVector tTrans0(-0.5*(fLens[0]+shight),0,-fLens[2]+layer12);
    G4SubtractionSolid*  tubox = new G4SubtractionSolid("tubox", gftub, gcbox,new G4RotationMatrix(),tTrans1);
    G4SubtractionSolid*  gubox = new G4SubtractionSolid("gubox", tubox, gcbox,new G4RotationMatrix(),tTrans0);

    G4Tubs* gcylinder1  = new G4Tubs("Cylinder1",0,r1,0.5*fLens[1],0*deg,360*deg);
    G4Tubs* gcylinder2  = new G4Tubs("Cylinder2",0,r2,0.5*fLens[1]-0.5,0*deg,360*deg);
    G4Tubs* gcylinder1c = new G4Tubs("Cylinder1c",0,r1,0.5*fLens[1]+0.5,0*deg,360*deg);
    G4Tubs* gcylinder2c = new G4Tubs("Cylinder2c",0,r2,0.5*fLens[1]+0.5,0*deg,360*deg);
    G4RotationMatrix* xRot = new G4RotationMatrix();
    xRot->rotateX(M_PI/2.*rad);
    
    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Lens1", gubox, gcylinder1,xRot,zTrans1); 
    G4SubtractionSolid*  gLenst = new G4SubtractionSolid("temp", gubox, gcylinder1c, xRot,zTrans1);

    G4IntersectionSolid* gLens2 = new G4IntersectionSolid("Lens2", gLenst, gcylinder2, xRot,zTrans2);
    G4SubtractionSolid*  gLens3 = new G4SubtractionSolid("Lens3", gLenst, gcylinder2c, xRot,zTrans2);
     
    lLens1 = new G4LogicalVolume(gLens1,BarMaterial,"lLens1",0,0,0);
    lLens2 = new G4LogicalVolume(gLens2,Nlak33aMaterial,"lLens2",0,0,0);
    lLens3 = new G4LogicalVolume(gLens3,BarMaterial,"lLens3",0,0,0);
  }
   
  if(fLensId == 7){ // 3-component cylindrical lens
    G4double lensMinThikness = 2.0; 
  
    G4double r1 = 0; //PrtManager::Instance()->GetTest1();
    G4double r2 = 0; //PrtManager::Instance()->GetTest2();

    if(PrtManager::Instance()->GetRunType() == 6){ //focal plane scan
      r1 = PrtManager::Instance()->GetTest1();
      r2 = PrtManager::Instance()->GetTest2();
    }
    
    // r1 = (r1==0)? 27.45: r1;
    // r2 = (r2==0)? 20.02: r2;

    r1 = (r1==0)? 33: r1;
    r2 = (r2==0)? 25: r2;
    G4double shight = 20;

    G4ThreeVector zTrans1(0, 0, -r1-fLens[2]/2.+r1-sqrt(r1*r1-shight/2.*shight/2.) +lensMinThikness);
    G4ThreeVector zTrans2(0, 0, -r2-fLens[2]/2.+r2-sqrt(r2*r2-shight/2.*shight/2.) +lensMinThikness*2);

    G4Box* gfboxx = new G4Box("fboxx",0.5*fLens[0],0.5*fLens[1],0.5*fLens[2]);
    G4Box* gcbox = new G4Box("cbox",0.5*fLens[0],0.5*fLens[1]+1,0.5*fLens[2]);
    G4ThreeVector tTrans1( 0.5*(fLens[0]+shight),0,-fLens[2]+lensMinThikness);
    G4ThreeVector tTrans0(-0.5*(fLens[0]+shight),0,-fLens[2]+lensMinThikness);
    G4SubtractionSolid*  tubox = new G4SubtractionSolid("tubox", gfboxx, gcbox,new G4RotationMatrix(),tTrans1);
    G4SubtractionSolid*  gubox = new G4SubtractionSolid("gubox", tubox, gcbox,new G4RotationMatrix(),tTrans0);

    G4Tubs* gcylinder1  = new G4Tubs("Cylinder1",0,r1,0.5*fLens[1],0*deg,360*deg);
    G4Tubs* gcylinder2  = new G4Tubs("Cylinder2",0,r2,0.5*fLens[1]-0.5,0*deg,360*deg);
    G4Tubs* gcylinder1c = new G4Tubs("Cylinder1c",0,r1,0.5*fLens[1]+0.5,0*deg,360*deg);
    G4Tubs* gcylinder2c = new G4Tubs("Cylinder2c",0,r2,0.5*fLens[1]+0.5,0*deg,360*deg);
    G4RotationMatrix* xRot = new G4RotationMatrix();
    xRot->rotateX(M_PI/2.*rad);
    
    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Lens1", gubox, gcylinder1,xRot,zTrans1); 
    G4SubtractionSolid*  gLenst = new G4SubtractionSolid("temp", gubox, gcylinder1c, xRot,zTrans1);

    G4IntersectionSolid* gLens2 = new G4IntersectionSolid("Lens2", gLenst, gcylinder2, xRot,zTrans2);
    G4SubtractionSolid*  gLens3 = new G4SubtractionSolid("Lens3", gLenst, gcylinder2c, xRot,zTrans2);
     
    lLens1 = new G4LogicalVolume(gLens1,BarMaterial,"lLens1",0,0,0);
    lLens2 = new G4LogicalVolume(gLens2,Nlak33aMaterial,"lLens2",0,0,0);
    lLens3 = new G4LogicalVolume(gLens3,BarMaterial,"lLens3",0,0,0);
  }

  if(fLensId == 8){ // 3-component cylindrical lens
    G4double lensMinThikness = 2.0; 
    G4double r1 = 0; //PrtManager::Instance()->GetTest1();
    G4double r2 = 0; //PrtManager::Instance()->GetTest2();

    if(PrtManager::Instance()->GetRunType() == 6){ //focal plane scan
      r1 = PrtManager::Instance()->GetTest1();
      r2 = PrtManager::Instance()->GetTest2();
    }
    
    // // thickness scan
    // G4double d = PrtManager::Instance()->GetTest1();
    // d = (d==0)? 3: d;    
    
   // r1 = (r1==0)? 27.45: r1;
   // r2 = (r2==0)? 20.02: r2;

    r1 = (r1==0)? 33: r1;
    r2 = (r2==0)? 25: r2;
    G4double shight = 0; //19

    G4ThreeVector zTrans1(0, 0, -r1-fLens[2]/2.+r1-sqrt(r1*r1-shight/2.*shight/2.) +3.0); //1.5
    G4ThreeVector zTrans2(0, 0, -r2-fLens[2]/2.+r2-sqrt(r2*r2-shight/2.*shight/2.) +3.0+5);// 3.5

    G4Box* gfboxx = new G4Box("fboxx",0.5*fLens[0],0.5*fLens[1],0.5*fLens[2]);
    G4Box* gcbox = new G4Box("cbox",0.5*fLens[0],0.5*fLens[1]+1,0.5*fLens[2]);
    G4ThreeVector tTrans1( 0.5*(fLens[0]+shight),0,-fLens[2]);
    G4ThreeVector tTrans0(-0.5*(fLens[0]+shight),0,-fLens[2]);
    G4SubtractionSolid*  tubox = new G4SubtractionSolid("tubox", gfboxx, gcbox,new G4RotationMatrix(),tTrans1);
    G4SubtractionSolid*  gubox = new G4SubtractionSolid("gubox", tubox, gcbox,new G4RotationMatrix(),tTrans0);

    G4Tubs* gcylinder1  = new G4Tubs("Cylinder1",0,r1,0.5*fLens[1],0*deg,360*deg);
    G4Tubs* gcylinder2  = new G4Tubs("Cylinder2",0,r2,0.5*fLens[1]-0.5,0*deg,360*deg);
    G4Tubs* gcylinder1c = new G4Tubs("Cylinder1c",0,r1,0.5*fLens[1]+0.5,0*deg,360*deg);
    G4Tubs* gcylinder2c = new G4Tubs("Cylinder2c",0,r2,0.5*fLens[1]+0.5,0*deg,360*deg);
    G4RotationMatrix* xRot = new G4RotationMatrix();
    xRot->rotateX(M_PI/2.*rad);
    
    G4IntersectionSolid* gLens1 = new G4IntersectionSolid("Lens1", gubox, gcylinder1,xRot,zTrans1); 
    G4SubtractionSolid*  gLenst = new G4SubtractionSolid("temp", gubox, gcylinder1c, xRot,zTrans1);

    G4IntersectionSolid* gLens2 = new G4IntersectionSolid("Lens2", gLenst, gcylinder2, xRot,zTrans2);
    G4SubtractionSolid*  gLens3 = new G4SubtractionSolid("Lens3", gLenst, gcylinder2c, xRot,zTrans2);
     
    lLens1 = new G4LogicalVolume(gLens1,BarMaterial,"lLens1",0,0,0);
    lLens2 = new G4LogicalVolume(gLens2,Nlak33aMaterial,"lLens2",0,0,0);
    lLens3 = new G4LogicalVolume(gLens3,BarMaterial,"lLens3",0,0,0);
  }

  if (fLensId == 9) { // TL ACR spherical lens

    double r1 = 50.08;
    double r2 = 41.7;
    double r3 = 247.7;
    double diameter = 50.8;
    double ht = 0.5*23;

    double lensMinThikness = 2;

    G4Tubs *gtub = new G4Tubs("tub", 0, 0.5 * diameter, 50, 0. * deg, 360. * deg);
    G4Sphere *gsphere1 = new G4Sphere("sphere1", 0, r1, 0, 360 * deg, 0, 360 * deg);
    G4Sphere *gsphere2 = new G4Sphere("sphere2", 0, r2, 0, 360 * deg, 0, 360 * deg);
    G4Sphere *gsphere3 = new G4Sphere("sphere3", 0, r3, 0, 360 * deg, 0, 360 * deg);
    auto r = new G4RotationMatrix();

    G4ThreeVector zTrans1(0, 0, r1 -ht + 3);
    G4ThreeVector zTrans2(0, 0, -r2 + ht );
    G4ThreeVector zTrans3(0, 0, -r3 + ht + 3);
    
    auto gLens01 = new G4IntersectionSolid("p01", gtub, gsphere2, r, zTrans2);
    auto gLens1 = new G4IntersectionSolid("p1", gLens01, gsphere1, r, zTrans1);
    auto gLens02 = new G4SubtractionSolid("p02", gtub, gsphere2, r, zTrans2);
    auto gLens2 = new G4IntersectionSolid("p2", gLens02, gsphere3, r, zTrans3);

    lLens1 = new G4LogicalVolume(gLens1, Nlak33aMaterial, "lLens1", 0, 0, 0);
    lLens2 = new G4LogicalVolume(gLens2, BarMaterial, "lLens2", 0, 0, 0);
  }

  fRotAngle = PrtManager::Instance()->GetAngle()*deg; 
  fPrtRot->rotateY(fRotAngle);
  fPrtRot->rotateX(PrtManager::Instance()->GetPhi()*deg);

  new G4PVPlacement(fPrtRot,G4ThreeVector(0,0,0.5*fLens[2]),lLens1,"wLens1", lTank,false,0);
  if(fLensId != 4) new G4PVPlacement(fPrtRot,G4ThreeVector(0,0,0.5*fLens[2]),lLens2,"wLens2", lTank,false,0);
  if(fLensId == 3 || fLensId==6 || fLensId==7 || fLensId==8)  new G4PVPlacement(fPrtRot,G4ThreeVector(0,0,0.5*fLens[2]),lLens3,"wLens3", lTank,false,0);

  auto gPixel = new G4Box("gPixel",100,100,1);
  lPixel = new G4LogicalVolume(gPixel,BarMaterial,"lPixel",0,0,0);
  if(fPixelZ < 0.5*fTank[2]) new G4PVPlacement(0,G4ThreeVector(0,0,fPixelZ),lPixel,"wPixel", lTank,false,0); 
  else new G4PVPlacement(0,G4ThreeVector(0,0,fPixelZ),lPixel,"wPixel", lExpHall,false,0);
  
  SetVisualization();

  return wExpHall;
}

void PrtDetectorConstruction::DefineMaterials(){
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;  

  G4int ncomponents, natoms;
  G4double fractionmass;

  // define Elements
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* Si = new G4Element("Silicon" ,symbol="Si", z= 14., a= 28.09*g/mole);

  G4Element* Al = new G4Element("Aluminum",symbol="Al",z=13.,a=26.98*g/mole);

  // quartz material = SiO2
  G4Material* SiO2 = new G4Material("quartz",density= 2.200*g/cm3, ncomponents=2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O , natoms=2);

  Nlak33aMaterial  = new G4Material("Nlak33a",density= 4.220*g/cm3, ncomponents=2);
  Nlak33aMaterial->AddElement(Si, natoms=1);
  Nlak33aMaterial->AddElement(O , natoms=2);

  PbF2Material  = new G4Material("PbF2",density= 4.220*g/cm3, ncomponents=2);
  PbF2Material->AddElement(Si, natoms=1);
  PbF2Material->AddElement(O , natoms=2);

  G4Material* Vacuum = new G4Material("interGalactic", 1., 1.008*g/mole, 
				      1.e-25*g/cm3, kStateGas, 
				      2.73*kelvin, 3.e-18*pascal);
  G4Material* Air = new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material* Aluminum = new G4Material("Aluminum",density=2.7*g/cm3,ncomponents=1);
  Aluminum->AddElement(Al,fractionmass=1.0);

  G4Material* KamLandOil = new G4Material("KamLandOil",density=0.914*g/cm3,ncomponents=2);
  KamLandOil->AddElement(C,natoms=12);
  KamLandOil->AddElement(H,natoms=26);

  G4Material* CarbonFiber = new G4Material("CarbonFiber", density=0.145*g/cm3, ncomponents=1);
  CarbonFiber->AddElement(C,fractionmass=1.0);
			
  G4Material* Epotek = new G4Material("Epotek",density=1.2*g/cm3,ncomponents=3);

  Epotek->AddElement(C,natoms=3);
  Epotek->AddElement(H,natoms=5);
  Epotek->AddElement(O,natoms=2);


  // assign main materials
  defaultMaterial = Air; //Vacuum // material of world
  frontMaterial = CarbonFiber; 
  BarMaterial = SiO2; // material of all Bars, Quartz and Window
  OilMaterial = KamLandOil; // material of volume 1,2,3,4
  MirrorMaterial = Aluminum; // mirror material
  epotekMaterial = Epotek; // Epotek material - glue between bars


  // ------------ Generate & Add Material Properties Table ------------

  static const G4double LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;
  const G4int num = 36;
  G4double WaveLength[num];
  G4double AirAbsorption[num]; // absorption value for air
  G4double AirRefractiveIndex[num]; // air refractive index
  G4double PhotonEnergy[num]; // energy of photons which correspond to the given 
  // refractive or absoprtion values

  G4double PhotonEnergyNlak33a[76] = {1,1.2511,1.26386,1.27687,1.29016,1.30372,1.31758,1.33173,1.34619,1.36097,1.37607,1.39152,1.40731,1.42347,1.44,1.45692,1.47425,1.49199,1.51016,1.52878,1.54787,1.56744,1.58751,1.6081,1.62923,1.65092,1.6732,1.69609,1.71961,1.7438,1.76868,1.79427,1.82062,1.84775,1.87571,1.90452,1.93423,1.96488,1.99652,2.0292,2.06296,2.09787,2.13398,2.17135,2.21006,2.25017,2.29176,2.33492,2.37973,2.42631,2.47473,2.52514,2.57763,2.63236,2.68946,2.7491,2.81143,2.87666,2.94499,3.01665,3.09187,3.17095,3.25418,3.34189,3.43446,3.53231,3.6359,3.74575,3.86244,3.98663,4.11908,4.26062,4.41225,4.57506,4.75035,4.93961};

  /*************************** ABSORPTION COEFFICIENTS *****************************/

  // absorption of KamLandOil per 50 cm - from jjv
  G4double KamLandOilAbsorption[num]=
    {0.97469022,0.976603956,0.978511548,0.980400538,0.982258449,0.984072792,
     0.985831062,0.987520743,0.989129303,0.990644203,0.992052894,
     0.993342822,0.994501428,0.995516151,0.996374433,0.997063719,
     0.997571464,0.997885132,0.997992205,0.997880183,0.997536591,
     0.99,0.98,0.97,0.96,0.94,0.93,0.924507,0.89982,0.883299,
     0.85657,0.842637,0.77020213,0.65727,0.324022,0.019192};

  // absorption of quartz per 1 m - from jjv
  G4double QuartzAbsorption[num] = 
    {0.999572036,0.999544661,0.999515062,0.999483019,0.999448285,
     0.999410586,0.999369611,0.999325013,0.999276402,0.999223336,
     0.999165317,0.999101778,0.999032079,0.998955488,0.998871172,
     0.998778177,0.99867541,0.998561611,0.998435332,0.998294892,0.998138345,
     0.997963425,0.997767484,0.997547418,
     0.99729958,0.99701966,0.99670255,0.996342167,0.995931242,0.995461041,
     0.994921022,0.994298396,0.993577567,0.992739402,0.991760297,0.990610945};
  
  // absorption of epotek per one layer - thicknes 0.001'' - from jjv
  G4double EpotekAbsorption[num] = 
    {0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.99999999,0.99999999,0.99999999,0.99999999,0.99999999,
     0.9999,0.9998,0.9995,0.999,0.998,0.997,0.996,0.9955,0.993,
     0.9871,0.9745};

  //N-Lak 33a
  G4double Nlak33aAbsorption[76]={371813,352095,331021,310814,291458,272937,255238,238342,222234,206897,192313,178463,165331,152896,141140,130043,119585,109747,100507,91846.3,83743.1,76176.7,69126.1,62570.2,56488,50858.3,45660.1,40872.4,36474.6,32445.8,28765.9,25414.6,22372.2,19619.3,17136.9,14906.5,12910.2,11130.3,9550.13,8153.3,6924.25,5848.04,4910.46,4098.04,3398.06,2798.54,2288.32,1856.99,1494.92,1193.28,943.973,739.657,573.715,440.228,333.94,250.229,185.064,134.967,96.9664,68.5529,47.6343,32.4882,21.7174,14.2056,9.07612,5.65267,3.4241,2.01226,1.14403,0.62722,0.330414,0.166558,0.0799649,0.0363677,0.0155708,0.00623089};

  // NLak33b from refractiveindex for 1 cm
  G4double Nlak33bEn[25]={0.4959, 0.5332 ,0.6293 ,0.8103 ,1.1696 ,1.7712 ,1.8785,1.9997,2.1376,2.2707,2.4796,2.6953,2.8436,2.9520,3.0613,3.0996,3.1790,3.2627,3.3509,3.3968,3.5424,3.7121,3.8745,3.9994,4.1328};
  G4double Nlak33bAb[25]={0.398114,0.679068,0.937060,0.985032,0.997996,0.997996,0.997595,0.997194,0.997595,0.997996,0.997194,0.994376,0.991546,0.988297,0.982161,0.979691,0.971388,0.954455,0.928177,0.910019,0.820600,0.657099,0.455454,0.245954,0.158490};

  G4int n_PbF2=56;
  G4double en_PbF2[] = {1.55 ,1.569,1.59 ,1.61 ,1.631,1.653,1.675,1.698,1.722,1.746,1.771,1.797,1.823,1.851,1.879,1.907,1.937,1.968,2    ,2.033,2.066,2.101,2.138,2.175,2.214,2.254,2.296,2.339,2.384,2.431,2.48 ,2.53 ,2.583,2.638,2.695,2.755,2.818,2.883,2.952,3.024,3.1  ,3.179,3.263,3.351,3.444,3.542,3.647,3.757,3.875,3.999,4.133,4.275,4.428,4.592,4.769,4.959};

  G4double ab_PbF2[]= {407  ,403.3,379.1,406.3,409.7,408.9,406.7,404.7,391.7,397.7,409.6,403.7,403.8,409.7,404.9,404.2,407.1,411.1,403.1,406.1,415.4,399.1,405.8,408.2,385.7,405.6,405.2,401.6,402.6,407.1,417.7,401.1,389.9,411.9,400.9,398.3,402.1,408.7,384.8,415.8,413.1,385.7,353.7,319.1,293.6,261.9,233.6,204.4,178.3,147.6,118.2,78.7 ,51.6 ,41.5 ,24.3 ,8.8};
  G4double ref_PbF2[]= {1.749,1.749,1.75 ,1.75 ,1.751,1.752,1.752,1.753,1.754,1.754,1.755,1.756,1.757,1.757,1.758,1.759,1.76 ,1.761,1.762,1.764,1.765,1.766,1.768,1.769,1.771,1.772,1.774,1.776,1.778,1.78 ,1.782,1.785,1.787,1.79 ,1.793,1.796,1.8  ,1.804,1.808,1.813,1.818,1.824,1.83 ,1.837,1.845,1.854,1.865,1.877,1.892,1.91 ,1.937,1.991,1.38 ,1.915,1.971,2.019};
  
  for(int i=0;i<num;i++){
    WaveLength[i]= (300 +i*10)*nanometer;
    //    AirAbsorption[i] = 4*cm; // if photon in the air -> kill it immediately
    AirAbsorption[i] = 400000000*cm; //rd for air gap 
    AirRefractiveIndex[i] = 1; 
    PhotonEnergy[num-(i+1)]= LambdaE/WaveLength[i];

    /* as the absorption is given per length and G4 needs 
       mean free path length, calculate it here
       mean free path length - taken as probability equal 1/e
       that the photon will be absorbed */
      
    EpotekAbsorption[i] = (-1)/log(EpotekAbsorption[i])*0.001*2.54*cm;
    QuartzAbsorption[i] = (-1)/log(QuartzAbsorption[i])*100*cm;
    KamLandOilAbsorption[i] = (-1)/log(KamLandOilAbsorption[i])*50*cm;
  }
  for(Int_t i=0; i<25; i++){
    Nlak33bEn[i] *= eV;
    Nlak33bAb[i] = (-0.8)/log(Nlak33bAb[i])*1*cm; // account for glue in lens  
  }

  
  /**************************** REFRACTIVE INDEXES ****************************/
  
  // only phase refractive indexes are necessary -> g4 calculates group itself !!
  
  G4double QuartzRefractiveIndex[num]={
    1.456535,1.456812,1.4571,1.457399,1.457712,1.458038,1.458378,
    1.458735,1.459108,1.4595,1.459911,1.460344,1.460799,1.46128,
    1.461789,1.462326,1.462897,1.463502,1.464146,1.464833,
    1.465566,1.46635,1.46719,1.468094,1.469066,1.470116,1.471252,1.472485,
    1.473826,1.475289,1.476891,1.478651,1.480592,1.482739,1.485127,1.487793};

  G4double EpotekRefractiveIndex[num]={
    1.554034,1.555575,1.55698,1.558266,1.559454,1.56056,1.561604,
    1.562604,1.563579,1.564547,1.565526,1.566536,1.567595,
    1.568721,1.569933,1.57125,1.57269,1.574271,1.576012,
    1.577932,1.580049,1.582381,1.584948,1.587768,1.590859,
    1.59424,1.597929,1.601946,1.606307,1.611033,1.616141,1.621651,1.62758,
    1.633947,1.640771,1.64807};

  G4double KamLandOilRefractiveIndex[num]={
    1.433055,1.433369,1.433698,1.434045,1.434409,1.434793,1.435198,
    1.435626,1.436077,1.436555,1.4371,1.4376,1.4382,1.4388,1.4395,
    1.4402,1.4409,1.4415,1.4425,1.4434,1.4444,1.4455,1.4464,1.4479,1.4501,
    1.450428,1.451976,1.453666,1.455513,1.45754,1.45977,1.462231,1.464958,
    1.467991,1.471377,1.475174};

  G4double Nlak33aRefractiveIndex[76]={1.73816,1.73836,1.73858,1.73881,1.73904,1.73928,1.73952,1.73976,1.74001,1.74026,1.74052,1.74078,1.74105,1.74132,1.7416,1.74189,1.74218,1.74249,1.74279,1.74311,1.74344,1.74378,1.74412,1.74448,1.74485,1.74522,1.74562,1.74602,1.74644,1.74687,1.74732,1.74779,1.74827,1.74878,1.7493,1.74985,1.75042,1.75101,1.75163,1.75228,1.75296,1.75368,1.75443,1.75521,1.75604,1.75692,1.75784,1.75882,1.75985,1.76095,1.76211,1.76335,1.76467,1.76608,1.76758,1.7692,1.77093,1.77279,1.7748,1.77698,1.77934,1.7819,1.7847,1.78775,1.79111,1.79481,1.79889,1.80343,1.8085,1.81419,1.82061,1.8279,1.83625,1.84589,1.85713,1.87039};

  /* ASSIGNING REFRACTIVE AND ABSORPTION PROPERTIES TO THE GIVEN MATERIALS */

  // Quartz material => Si02
  G4MaterialPropertiesTable* QuartzMPT = new G4MaterialPropertiesTable();
  QuartzMPT->AddProperty("RINDEX",       PhotonEnergy, QuartzRefractiveIndex,num);
  QuartzMPT->AddProperty("ABSLENGTH",    PhotonEnergy, QuartzAbsorption,           num);
  BarMaterial->SetMaterialPropertiesTable(QuartzMPT);
  
  // Air
  G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
  AirMPT->AddProperty("RINDEX",    PhotonEnergy, AirRefractiveIndex, num);
  AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAbsorption,      num);
  //  assign this parameter table to the air 
  defaultMaterial->SetMaterialPropertiesTable(AirMPT);

  // KamLandOil                                                
  G4MaterialPropertiesTable* KamLandOilMPT = new G4MaterialPropertiesTable();
  KamLandOilMPT->AddProperty("RINDEX", PhotonEnergy, KamLandOilRefractiveIndex, num);
  KamLandOilMPT->AddProperty("ABSLENGTH", PhotonEnergy, KamLandOilAbsorption, num);
  OilMaterial->SetMaterialPropertiesTable(KamLandOilMPT);  

  // N-Lak 33a                                                
  for(Int_t i=0; i<76; i++){
    PhotonEnergyNlak33a[i]*=eV;
    Nlak33aAbsorption[i]*=cm; // cm to mm
  }
  G4MaterialPropertiesTable* Nlak33aMPT = new G4MaterialPropertiesTable();
  Nlak33aMPT->AddProperty("RINDEX", PhotonEnergyNlak33a, Nlak33aRefractiveIndex, 76);
  //Nlak33aMPT->AddProperty("ABSLENGTH",PhotonEnergyNlak33a, Nlak33aAbsorption, 76);
  Nlak33aMPT->AddProperty("ABSLENGTH",Nlak33bEn, Nlak33bAb, 25);
  Nlak33aMaterial->SetMaterialPropertiesTable(Nlak33aMPT);

  // PbF2
  G4MaterialPropertiesTable* PbF2MPT = new G4MaterialPropertiesTable();
  PbF2MPT->AddProperty("RINDEX", en_PbF2, ref_PbF2, n_PbF2);
  PbF2MPT->AddProperty("ABSLENGTH",en_PbF2, ab_PbF2, n_PbF2);
  PbF2Material->SetMaterialPropertiesTable(PbF2MPT);

  // Optical grease                                                
  G4MaterialPropertiesTable* opticalGreaseMPT = new G4MaterialPropertiesTable();
  G4double og_en[6]={1*eV,2*eV,3*eV,4*eV,4.2*eV,10*eV};
  G4double og_ab[6]={0.66*cm,0.66*cm,0.47*cm,0.14*cm,0.06*cm,0.02*cm};
  G4double og_re[6]={1.55,1.56,1.59,1.64,1.64,1.64};
  opticalGreaseMPT->AddProperty("RINDEX", og_en, og_re, 6);
  opticalGreaseMPT->AddProperty("ABSLENGTH",og_en, og_ab, 6);

  opticalGreaseMaterial = new G4Material("opticalGreaseMaterial",density= 2.200*g/cm3, ncomponents=2);
  opticalGreaseMaterial->AddElement(Si, natoms=1);
  opticalGreaseMaterial->AddElement(O , natoms=2);
  opticalGreaseMaterial->SetMaterialPropertiesTable(opticalGreaseMPT);  
  
  // Optical cookie (RTV615)                                                
  G4MaterialPropertiesTable* opticalCookieMPT = new G4MaterialPropertiesTable();
  G4double oc_en[9]={1.50*eV,2.00*eV,2.50*eV,3.00*eV,3.50*eV,4.00*eV,4.10*eV,4.50*eV,5.00*eV};
  G4double oc_ab[9]={14.2*cm,14.2*cm,14.2*cm,11.54*cm,5.29*cm,2.98*cm,2.43*cm,2.43*cm,2.43*cm};
  G4double oc_re[9]={1.406,1.406,1.406,1.406,1.406,1.406,1.406,1.406,1.406};
  opticalCookieMPT->AddProperty("RINDEX", oc_en, oc_re, 9);
  opticalCookieMPT->AddProperty("ABSLENGTH",oc_en, oc_ab, 9);

  opticalCookieMaterial = new G4Material("opticalCookieMaterial",density= 2.200*g/cm3, ncomponents=2);
  opticalCookieMaterial->AddElement(Si, natoms=1);
  opticalCookieMaterial->AddElement(O , natoms=2);
  opticalCookieMaterial->SetMaterialPropertiesTable(opticalCookieMPT);  

  // Epotek Glue                                        
  G4MaterialPropertiesTable* EpotekMPT = new G4MaterialPropertiesTable();
  EpotekMPT->AddProperty("RINDEX", PhotonEnergy, EpotekRefractiveIndex, num);
  EpotekMPT->AddProperty("ABSLENGTH", PhotonEnergy, EpotekAbsorption, num);
  // assign this parameter table to the epotek
  epotekMaterial->SetMaterialPropertiesTable(EpotekMPT);
}

void PrtDetectorConstruction::SetVisualization(){

  G4Colour blue = G4Colour(0.0,0.0,1.0);
  G4Colour green = G4Colour(0.0,1.0,.0);
  G4Colour red = G4Colour(1.0,0.0,.0); 
  G4Colour DircColour = G4Colour(1.,1.0,0.);
  G4Colour Dark = G4Colour(0.,0.05,0.05,0.15);

  G4VisAttributes *vaExpHall = new G4VisAttributes(G4Colour(0.,1.,1.,0));
  vaExpHall->SetForceWireframe(true);
  //vaExpHall->SetVisibility(false);
  lExpHall->SetVisAttributes(vaExpHall);

  G4VisAttributes *vaTank = new G4VisAttributes(blue);
  //vaTank->SetVisibility(false);
  vaTank->SetForceWireframe(true);
  lTank->SetVisAttributes(vaTank);
 
  G4double transp = 0.15;
  if(fLensId!=0 && fLensId!=10){
    G4VisAttributes * vaLens = new G4VisAttributes(G4Colour(0.,1.,1.,transp));
    //vaLens->SetForceWireframe(true);
    //vaLens->SetForceAuxEdgeVisible(true);
    lLens1->SetVisAttributes(vaLens);
    G4VisAttributes * vaLens1 = new G4VisAttributes(G4Colour(0.,0.5,1.0,transp));
    if(fLensId!=4) lLens2->SetVisAttributes(vaLens1);
    if(fLensId==2) {
      lLens1->SetVisAttributes(vaLens1);
      lLens2->SetVisAttributes(vaLens);
    }
    if(fLensId==3 || fLensId==6  || fLensId==7  || fLensId==8) lLens3->SetVisAttributes(vaLens);
  }

  G4VisAttributes *waPixel = new G4VisAttributes(red);
  //waPixel->SetForceWireframe(true);
  lPixel->SetVisAttributes(waPixel);
}

void PrtDetectorConstruction::ConstructSDandField(){

  PrtPixelSD* pixelSD = new PrtPixelSD("PixelSD", "PixelHitsCollection", 0);
  G4SDManager::GetSDMpointer()->AddNewDetector(pixelSD);
  SetSensitiveDetector("lPixel",pixelSD);
}

void PrtDetectorConstruction::SetRotation(G4double angle){

  fPrtRot->rotateY(-fRotAngle);
  fPrtRot->rotateY(angle);
  fRotAngle=angle;
  
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void PrtDetectorConstruction::SetLens(G4int id){
 
}

void PrtDetectorConstruction::SetQuantumEfficiency(G4int id){

}
