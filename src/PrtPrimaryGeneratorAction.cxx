#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "PrtPrimaryGeneratorAction.h"
#include "PrtPrimaryGeneratorMessenger.h"
#include "PrtManager.h"

PrtPrimaryGeneratorAction::PrtPrimaryGeneratorAction():G4VUserPrimaryGeneratorAction(),fParticleGun(0){

  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  fGunMessenger = new PrtPrimaryGeneratorMessenger(this);
}

PrtPrimaryGeneratorAction::~PrtPrimaryGeneratorAction(){
  delete fParticleGun;
  delete fGunMessenger;
}

void PrtPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){

  PrtManager::Instance()->AddEvent(PrtEvent());
  double sigma = PrtManager::Instance()->GetBeamDimension();
  double divergence = PrtManager::Instance()->GetBeamDivergence();
  double events = PrtManager::Instance()->GetEvents();
  double beamx = PrtManager::Instance()->GetBeamX();
  double beamy = PrtManager::Instance()->GetBeamY();

  for(int i=0; i<events; i++){
    double x = G4RandGauss::shoot(beamx,sigma);
    double y = G4RandGauss::shoot(beamy,sigma);
    G4ThreeVector dir = G4ThreeVector(0,0,1); //fParticleGun->GetParticleMomentumDirection();

    dir *= fParticleGun->GetParticleMomentum();
    dir.setTheta(G4RandGauss::shoot(0, divergence)); // beam divergence
    dir.setPhi(2 * M_PI * G4UniformRand());

    fParticleGun->SetParticlePosition(G4ThreeVector(x, y,-200));
    fParticleGun->SetParticleMomentumDirection(dir);
    fParticleGun->GeneratePrimaryVertex(anEvent);

  }
  G4ThreeVector dir = G4ThreeVector(0,0,1);
  PrtManager::Instance()->SetMomentum(TVector3(dir.x(),dir.y(),dir.z()));

}

void PrtPrimaryGeneratorAction::SetOptPhotonPolar(){
 G4double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}

void PrtPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle){
 if (fParticleGun->GetParticleDefinition()->GetParticleName()!="opticalphoton"){
     G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
       "the particleGun is not an opticalphoton " << 
       fParticleGun->GetParticleDefinition()->GetParticleName()<< G4endl;
     return;
   }

 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton);
 G4double modul2       = product*product;
 
 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
 
 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 fParticleGun->SetParticlePolarization(polar);
}


