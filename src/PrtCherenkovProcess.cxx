#include "PrtCherenkovProcess.h"

#include "PrtManager.h"

#include "G4Poisson.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

PrtCherenkovProcess::PrtCherenkovProcess(const G4String& processName, G4ProcessType type)
  : G4Cerenkov(processName, type)
{
 
}

G4VParticleChange* PrtCherenkovProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine is called for each tracking Step of a charged particle
// in a radiator. A Poisson-distributed number of photons is generated
// according to the Cerenkov formula, distributed evenly along the track
// segment and uniformly azimuth w.r.t. the particle direction. The 
// parameters are then transformed into the Master Reference System, and 
// they are added to the particle change. 

{
	//////////////////////////////////////////////////////
	// Should we ensure that the material is dispersive?
	//////////////////////////////////////////////////////

        aParticleChange.Initialize(aTrack);

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        const G4Material* aMaterial = aTrack.GetMaterial();

	G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
	G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

	G4ThreeVector x0 = pPreStepPoint->GetPosition();
        G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
	G4double t0 = pPreStepPoint->GetGlobalTime();

        G4MaterialPropertiesTable* aMaterialPropertiesTable =
                               aMaterial->GetMaterialPropertiesTable();
        if (!aMaterialPropertiesTable) return pParticleChange;

	G4MaterialPropertyVector* Rindex = 
                aMaterialPropertiesTable->GetProperty("RINDEX"); 
        if (!Rindex) return pParticleChange;

        // particle charge
        const G4double charge = aParticle->GetDefinition()->GetPDGCharge();

        // particle beta
        const G4double beta = (pPreStepPoint ->GetBeta() +
                               pPostStepPoint->GetBeta())/2.;

	G4double MeanNumberOfPhotons = 
                 GetAverageNumberOfPhotons(charge,beta,aMaterial,Rindex);

        if (MeanNumberOfPhotons <= 0.0) {

                // return unchanged particle and no secondaries

                aParticleChange.SetNumberOfSecondaries(0);
 
                return pParticleChange;

        }

        G4double step_length;
        step_length = aStep.GetStepLength();

	MeanNumberOfPhotons = MeanNumberOfPhotons * step_length;

	G4int NumPhotons = (G4int) G4Poisson(MeanNumberOfPhotons);

	if (NumPhotons <= 0) {

		// return unchanged particle and no secondaries  

		aParticleChange.SetNumberOfSecondaries(0);
		
                return pParticleChange;
	}

	////////////////////////////////////////////////////////////////

	aParticleChange.SetNumberOfSecondaries(NumPhotons);

        if (20) { //fTrackSecondariesFirst
           if (aTrack.GetTrackStatus() == fAlive )
                   aParticleChange.ProposeTrackStatus(fSuspend);
        }
	
	////////////////////////////////////////////////////////////////

	// G4double Pmin = Rindex->GetMinLowEdgeEnergy();
	// G4double Pmax = Rindex->GetMaxLowEdgeEnergy();
	
	// //d/ monochromatic photons
	G4double Pmin = 3.18*1E-6;//(Rindex->GetMinLowEdgeEnergy()+Rindex->GetMaxLowEdgeEnergy())/2.;
	G4double Pmax = Pmin;
 
	G4double dp = Pmax - Pmin;

	G4double nMax = Rindex->GetMaxValue();

        G4double BetaInverse = 1./beta;

	G4double maxCos = BetaInverse / nMax; 
	G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);

        const G4double beta1 = pPreStepPoint ->GetBeta();
        const G4double beta2 = pPostStepPoint->GetBeta();

        G4double MeanNumberOfPhotons1 =
                     GetAverageNumberOfPhotons(charge,beta1,aMaterial,Rindex);
        G4double MeanNumberOfPhotons2 =
                     GetAverageNumberOfPhotons(charge,beta2,aMaterial,Rindex);

	for (G4int i = 0; i < NumPhotons; i++) {

		// Determine photon energy

		G4double rand;
		G4double sampledEnergy, sampledRI; 
		G4double cosTheta, sin2Theta;
		
		// sample an energy

		do {
			rand = G4UniformRand();	
			sampledEnergy = Pmin + rand * dp; 
			sampledRI = Rindex->Value(sampledEnergy);
			cosTheta = BetaInverse / sampledRI;  
			if(cosTheta>1) {
			  std::cout<<"Warning - PrtCherenkovProcess:  cosTheta "<< cosTheta <<std::endl;
			   return pParticleChange;
			}

			sin2Theta = (1.0 - cosTheta)*(1.0 + cosTheta);
			rand = G4UniformRand();	
 
		} while (rand*maxSin2 > sin2Theta);

		// Generate random position of photon on cone surface 
		// defined by Theta 

		rand = G4UniformRand();

		G4double phi = twopi*rand; //twopi*PrtManager::Instance()->GetTest()/360.;//twopi*rand;
		G4double sinPhi = std::sin(phi);
		G4double cosPhi = std::cos(phi);

		// calculate x,y, and z components of photon energy
		// (in coord system with primary particle direction 
		//  aligned with the z axis)

		G4double sinTheta = std::sqrt(sin2Theta); 
		G4double px = sinTheta*cosPhi;
		G4double py = sinTheta*sinPhi;
		G4double pz = cosTheta;

		// Create photon momentum direction vector 
		// The momentum direction is still with respect
	 	// to the coordinate system where the primary
		// particle direction is aligned with the z axis  

		G4ParticleMomentum photonMomentum(px, py, pz);

		// Rotate momentum direction back to global reference
		// system 

                photonMomentum.rotateUz(p0); //test

		// Determine polarization of new photon 

		G4double sx = cosTheta*cosPhi;
		G4double sy = cosTheta*sinPhi; 
		G4double sz = -sinTheta;

		G4ThreeVector photonPolarization(sx, sy, sz);

		// Rotate back to original coord system 

                photonPolarization.rotateUz(p0);
		
                // Generate a new photon:

                G4DynamicParticle* aCerenkovPhoton =
                  new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), 
  					                 photonMomentum);
		aCerenkovPhoton->SetPolarization
				     (photonPolarization.x(),
				      photonPolarization.y(),
				      photonPolarization.z());

		aCerenkovPhoton->SetKineticEnergy(sampledEnergy);

                // Generate new G4Track object:

                G4double delta, NumberOfPhotons, N;

                do {
		  rand = G4UniformRand();//PrtManager::Instance()->GetTest();//G4UniformRand();
		  delta = rand * aStep.GetStepLength();
		  NumberOfPhotons = MeanNumberOfPhotons1 - delta *
		    (MeanNumberOfPhotons1-MeanNumberOfPhotons2)/
		    aStep.GetStepLength();
		  N = G4UniformRand() *
		    std::max(MeanNumberOfPhotons1,MeanNumberOfPhotons2);
                } while (N > NumberOfPhotons);

		G4double deltaTime = delta /
                       ((pPreStepPoint->GetVelocity()+
                         pPostStepPoint->GetVelocity())/2.);

                G4double aSecondaryTime = t0 + deltaTime;

                G4ThreeVector aSecondaryPosition =
                                    x0 + rand * aStep.GetDeltaPosition();

		G4Track* aSecondaryTrack = 
		new G4Track(aCerenkovPhoton,aSecondaryTime,aSecondaryPosition);

                aSecondaryTrack->SetTouchableHandle(
                                 aStep.GetPreStepPoint()->GetTouchableHandle());

                aSecondaryTrack->SetParentID(aTrack.GetTrackID());

		aParticleChange.AddSecondary(aSecondaryTrack);
	}

	if (verboseLevel>0) {
	   G4cout <<"\n Exiting from G4Cerenkov::DoIt -- NumberOfSecondaries = "
	          << aParticleChange.GetNumberOfSecondaries() << G4endl;
	}

        return pParticleChange;
}

G4double PrtCherenkovProcess::GetAverageNumberOfPhotons(const G4double charge,
                              const G4double beta, 
			      const G4Material* aMaterial,
			      G4MaterialPropertyVector* Rindex) const
{
	const G4double Rfact = 369.81/(eV * cm);

        if(beta <= 0.0)return 0.0;

        G4double BetaInverse = 1./beta;

	// Vectors used in computation of Cerenkov Angle Integral:
	// 	- Refraction Indices for the current material
	//	- new G4PhysicsOrderedFreeVector allocated to hold CAI's
 
	G4int materialIndex = aMaterial->GetIndex();

	// Retrieve the Cerenkov Angle Integrals for this material  

	G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals =
	(G4PhysicsOrderedFreeVector*)((*thePhysicsTable)(materialIndex));

        if(!(CerenkovAngleIntegrals->IsFilledVectorExist()))return 0.0;

	// Min and Max photon energies 
	G4double Pmin = Rindex->GetMinLowEdgeEnergy();
	G4double Pmax = Rindex->GetMaxLowEdgeEnergy();

	// Min and Max Refraction Indices 
	G4double nMin = Rindex->GetMinValue();	
	G4double nMax = Rindex->GetMaxValue();

	// Max Cerenkov Angle Integral 
	G4double CAImax = CerenkovAngleIntegrals->GetMaxValue();

	G4double dp, ge;

	// If n(Pmax) < 1/Beta -- no photons generated 

	if (nMax < BetaInverse) {
		dp = 0;
		ge = 0;
	} 

	// otherwise if n(Pmin) >= 1/Beta -- photons generated  

	else if (nMin > BetaInverse) {
		dp = Pmax - Pmin;	
		ge = CAImax; 
	} 

	// If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then
	// we need to find a P such that the value of n(P) == 1/Beta.
	// Interpolation is performed by the GetEnergy() and
	// Value() methods of the G4MaterialPropertiesTable and
	// the GetValue() method of G4PhysicsVector.  

	else {
		Pmin = Rindex->GetEnergy(BetaInverse);
		dp = Pmax - Pmin;

		// need boolean for current implementation of G4PhysicsVector
		// ==> being phased out
		G4bool isOutRange;
		G4double CAImin = CerenkovAngleIntegrals->
                                  GetValue(Pmin, isOutRange);
		ge = CAImax - CAImin;

		if (verboseLevel>0) {
			G4cout << "CAImin = " << CAImin << G4endl;
			G4cout << "ge = " << ge << G4endl;
		}
	}
	
	// Calculate number of photons 
	G4double NumPhotons = Rfact * charge/eplus * charge/eplus *
                                 (dp - ge * BetaInverse*BetaInverse);

	return NumPhotons;		
}
