// simulation software for the Panda Barrel DIRC prototype
// contact: r.dzhygadlo@gsi.de

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "TROOT.h"

#include "PrtPhysicsList.h"
#include "PrtDetectorConstruction.h"

#include "PrtActionInitialization.h"
#include "time.h"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "TApplication.h"

#include "PrtManager.h"
#include "PrtLutReco.h"

namespace {
  void PrintUsage() {
    G4cerr<<"prtlens [OPTION] [ARGUMENT] ... [OPTION] [ARGUMENT]"<<G4endl;
    G4cerr<<"see README.md"<<G4endl;
  }
}

int main(int argc,char** argv)
{
  for ( G4int i=1; i<argc; i=i+2 ) std::cout<< argv[i] << "  "<<argv[i+1] <<std::endl;
  
  // Evaluate arguments
  if ( argc > 100 ) {
    PrintUsage();
    return 1;
  }
 
#ifdef G4MULTITHREADED
  G4int nThreads = 1;
#endif
  TApplication theApp("App", 0, 0);

  G4String macro,radiator, physlist,session,testVal1,testVal2,testVal3,
    prismStepX,prismStepY,beamX,timeRes,lutfile,mcpLayout
    ,geometry("1")
    ,outfile("focalplane.root")
    ,beamZ("200")
    ,beamDimension("1")
    ,lensId("3")
    ,particle("opticalphoton")
    ,geomAng("0")
    ,geomPhi("0")
    ,momentum("2 eV")
    ,infile("focalplane.root");
  G4int events(1),firstevent(0), runtype(6), study(0), verbose(0),batchmode(0);

  G4long myseed = 345354;
  for ( G4int i=1; i<argc; i=i+2 ) {
    if ( G4String(argv[i]) == "-seed" ) myseed    = atol(argv[i+1]);
    else if ( G4String(argv[i]) == "-o" ) outfile   = argv[i+1];
    else if ( G4String(argv[i]) == "-i" ) infile    = argv[i+1];

    else if ( G4String(argv[i]) == "-g" ) geometry  = argv[i+1];
    else if ( G4String(argv[i]) == "-l" ) lensId    = argv[i+1];
    else if ( G4String(argv[i]) == "-theta" ) geomAng   = argv[i+1];
    else if ( G4String(argv[i]) == "-phi" ) geomPhi   = argv[i+1];

    else if ( G4String(argv[i]) == "-p" ) particle  = argv[i+1];
    else if ( G4String(argv[i]) == "-m" ) momentum  = argv[i+1];
    else if ( G4String(argv[i]) == "-s" ) beamDimension  = argv[i+1];
    else if ( G4String(argv[i]) == "-e" ) events    = atoi(argv[i+1]);

    else if ( G4String(argv[i]) == "-z" ) beamZ    = argv[i+1];
    else if ( G4String(argv[i]) == "-b" ) batchmode = atoi(argv[i+1]);
    
    else if ( G4String(argv[i]) == "-u" ) lutfile   = argv[i+1];
    else if ( G4String(argv[i]) == "-h" ) radiator  = argv[i+1];
    else if ( G4String(argv[i]) == "-f" ) firstevent= atoi(argv[i+1]);
    else if ( G4String(argv[i]) == "-w" ) physlist  = argv[i+1];
    else if ( G4String(argv[i]) == "-r" ) runtype   = atoi(argv[i+1]);
    else if ( G4String(argv[i]) == "-study" ) study   = atoi(argv[i+1]);
    else if ( G4String(argv[i]) == "-c" ) mcpLayout = argv[i+1];
    else if ( G4String(argv[i]) == "-t1" ) testVal1   = argv[i+1];
    else if ( G4String(argv[i]) == "-t2" ) testVal2   = argv[i+1];
    else if ( G4String(argv[i]) == "-t3" ) testVal3   = argv[i+1];
    else if ( G4String(argv[i]) == "-gsx" ) prismStepX   = argv[i+1];
    else if ( G4String(argv[i]) == "-gsy" ) prismStepY   = argv[i+1];
    else if ( G4String(argv[i]) == "-gx" ) beamX    = argv[i+1];
    else if ( G4String(argv[i]) == "-tr" ) timeRes  = argv[i+1];
    else if ( G4String(argv[i]) == "-v" ) verbose   = atoi(argv[i+1]);
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }
  if(outfile=="") outfile = "focalplane.root";
  		   
  if(batchmode) gROOT->SetBatch(kTRUE);
  
  PrtManager::Instance(outfile,runtype);
  PrtManager::Instance()->SetEvents(events);
  PrtManager::Instance()->SetStudy(study);

  if(physlist.size()) PrtManager::Instance()->SetPhysList(atoi(physlist));
  if(geometry.size()) PrtManager::Instance()->SetGeometry(atoi(geometry));
  if(radiator.size()) PrtManager::Instance()->SetRadiator(atoi(radiator));
  if(lensId.size())   PrtManager::Instance()->SetLens(atoi(lensId));
  if(mcpLayout.size())PrtManager::Instance()->SetMcpLayout(atoi(mcpLayout));
  if(beamDimension.size())   PrtManager::Instance()->SetBeamDimension(atof(beamDimension));
  if(testVal1.size())   PrtManager::Instance()->SetShift(atof(testVal1));
  if(testVal1.size())   PrtManager::Instance()->SetTest1(atof(testVal1));
  if(testVal2.size())   PrtManager::Instance()->SetTest2(atof(testVal2));
  if(testVal3.size())   PrtManager::Instance()->SetTest3(atof(testVal3));
  if(geomAng.size())   PrtManager::Instance()->SetAngle(atof(geomAng));
  if(geomPhi.size())   PrtManager::Instance()->SetPhi(atof(geomPhi));
  if(prismStepX.size())   PrtManager::Instance()->SetPrismStepX(atof(prismStepX));
  if(prismStepY.size())   PrtManager::Instance()->SetPrismStepY(atof(prismStepY));
  if(beamX.size())   PrtManager::Instance()->SetBeamX(atof(beamX));
  if(beamZ.size())   PrtManager::Instance()->SetBeamZ(atof(beamZ));
  if(timeRes.size())   PrtManager::Instance()->SetTimeRes(atof(timeRes));

  if(runtype == 2 || runtype == 3 || runtype == 4){
    PrtLutReco * reco = new PrtLutReco(infile.c_str(),lutfile.c_str(),verbose); 
    reco->Run(firstevent, events);
    return 0;
  }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager * runManager = new G4RunManager;
#endif
  
  std::cout<<"SEED "<<myseed <<std::endl;
  G4Random::setTheSeed(myseed);

  
  // Detector construction
  runManager->SetUserInitialization(new PrtDetectorConstruction());
  // Physics list
  runManager->SetUserInitialization(new PrtPhysicsList());
  // User action initialization
  runManager->SetUserInitialization(new PrtActionInitialization(outfile));
  // Initialize G4 kernel
  runManager->Initialize();

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro.size() ) {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }else { 
    UImanager->ApplyCommand("/control/execute ../prt.mac");
  }
 
  if ( particle.size() ) {
    int pdgid = 0;
    if(particle=="mix_pip") PrtManager::Instance()->SetMix(1);
    else if(particle=="mix_pik") PrtManager::Instance()->SetMix(2);
    else{
      G4String command = "/gun/particle ";
      UImanager->ApplyCommand(command+particle);
      if(particle=="proton") pdgid = 2212;
      if(particle=="pi+") pdgid = 211;
      if(particle=="pi0") pdgid = 111;
      if(particle=="kaon+") pdgid = 321;
      if(particle=="kaon-") pdgid = -321;
      if(particle=="mu-") pdgid = 13;
      if(particle=="e-") pdgid = 11;
    }    
    PrtManager::Instance()->SetParticle(pdgid);
  }

  if(momentum.size()) UImanager->ApplyCommand( "/gun/momentumAmp "+momentum);

  if(batchmode){
    // batch mode
    UImanager->ApplyCommand("/run/beamOn 1");
  }else{
    // UI session for interactive mode
    G4UIExecutive * ui = new G4UIExecutive(argc,argv,"Qt");
    UImanager->ApplyCommand("/control/execute ../vis.mac");
    if (ui->IsGUI()) UImanager->ApplyCommand("/control/execute gui.mac");
    UImanager->ApplyCommand("/run/beamOn 1");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;
}

