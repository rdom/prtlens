// -----------------------------------------
// PrtOpBoundaryProcess.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtOpBoundaryProcess_h
#define PrtOpBoundaryProcess_h

#include "globals.hh"
#include "G4OpBoundaryProcess.hh"

class PrtOpBoundaryProcess : public G4OpBoundaryProcess
{
public:
  PrtOpBoundaryProcess();
  ~PrtOpBoundaryProcess(){};

public:
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

private:
  int fLensId;
};


#endif /*PrtOpBoundaryProcess_h*/
