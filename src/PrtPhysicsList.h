#ifndef PrtPhysicsList_h
#define PrtPhysicsList_h

#include "globals.hh"
#include "G4VUserPhysicsList.hh"

class G4Cerenkov;
class PrtCherenkovProcess;
class G4Scintillation;
class G4OpAbsorption;
class G4OpRayleigh;
class G4OpMieHG;
class G4OpBoundaryProcess;

class PrtPhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrtPhysicsList : public G4VUserPhysicsList
{ 
  public:
    PrtPhysicsList();
    virtual ~PrtPhysicsList();

  public:
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    virtual void SetCuts();

    //these methods Construct physics processes and register them
    void ConstructDecay();
    void ConstructEM();
    void ConstructOp();

    //for the Messenger 
    void SetVerbose(G4int);
    void SetNbOfPhotonsCerenkov(G4int);
 
  private:
    G4int                fPhysList;
    G4Cerenkov*          fCerenkovProcess;
    PrtCherenkovProcess* fCerenkovProcess0;
    G4Scintillation*     fScintillationProcess;
    G4OpAbsorption*      fAbsorptionProcess;
    G4OpRayleigh*        fRayleighScatteringProcess;
    G4OpMieHG*           fMieHGScatteringProcess;
    G4OpBoundaryProcess* fBoundaryProcess;
 
    PrtPhysicsListMessenger* fMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* PrtPhysicsList_h */
