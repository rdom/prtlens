#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

void fit_sim(TString infile="../build/focalplane.root", bool batch = 0){

  double p,px,py,gx,gy,z,x,theta;
  TTree *t = new TTree("data","data");
  t->Branch("px",&px,"px/D");
  t->Branch("py",&py,"py/D");
  t->Branch("gx",&gx,"gx/D");
  t->Branch("gy",&gy,"gy/D");
  t->Branch("z",&z,"z/D");
  t->Branch("x",&x,"x/D");
  t->Branch("theta",&theta,"theta/D");
  
  if(!prt_init(infile,1)) return;
  
  TH2F *h = new TH2F("h",";x [mm];y [mm]",200,-20,20,200,-20,20);
  TH1F *hx = new TH1F("hx","x [mm]",2000,-20,20);
  TH1F *hy = new TH1F("hy","y [mm]",2000,-20,20);
    
  TVector3 pos;
  
  for (int ievent=0; ievent < prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    z = prt_event->GetBeamZ();
    theta = prt_event->GetAngle();
    for(auto hit : prt_event->GetHits()){
      pos = hit.GetPosition();
      h->Fill(pos.X(),pos.Y());
      hx->Fill(pos.X());
      hy->Fill(pos.Y());
    }
  }

  double mx = hx->GetMean();
  double my = hy->GetMean();
  double sx = hx->GetStdDev();
  double sy = hy->GetStdDev();

  px = hx->GetMean();
  py = hy->GetMean();

  TFitResultPtr rx = hx->Fit("gaus", "SL", "", mx - 2 * sx, mx + 2 * sx);
  gx = rx->Parameter(2);
  TFitResultPtr ry = hy->Fit("gaus", "SL", "", my - 2 * sy, my + 2 * sy);
  gy = ry->Parameter(2);
  
  std::cout << "px " << px << " gx " << gx << " py " << py << " gy " << gy << std::endl;

  infile.ReplaceAll(".root", "_out.root");
  TFile f(infile, "recreate");
  t->Fill();
  t->Write();
  f.Close();
  
  if (!batch) {
    TString nid = Form("pos_%2.1f", z);
    prt_canvasAdd(nid, 500, 500);
    h->Draw("colz");
    prt_canvasAdd(nid + "_x", 500, 500);
    hx->Draw("colz");
    prt_canvasAdd(nid + "_y", 500, 500);
    hy->Draw("colz");
    prt_canvasSave("data/fit_sim", 0);
  }
}
