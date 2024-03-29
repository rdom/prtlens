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
  
  TH2F *h = new TH2F("h",";x [mm];y [mm]",200,-10,10,200,-5,5);
  TH1F *hx = new TH1F("hx","x [mm]",3000,-40,40);
  TH1F *hy = new TH1F("hy","y [mm]",3000,-40,40);
    
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

 
  double xmax = hx->GetXaxis()->GetBinCenter(hx->GetMaximumBin());
  int threshold = hx->GetMaximum() * 0.4;
  int firstbin = hx->FindFirstBinAbove(threshold);
  int lastbin = hx->FindLastBinAbove(threshold);
  double xf = hx->GetXaxis()->GetBinCenter(firstbin);
  double xl = hx->GetXaxis()->GetBinCenter(lastbin);
  px = xl - xf;  

  double ymax = hy->GetXaxis()->GetBinCenter(hy->GetMaximumBin());
  threshold = hy->GetMaximum() * 0.4;
  firstbin = hy->FindFirstBinAbove(threshold);
  lastbin = hy->FindLastBinAbove(threshold);
  double yf = hy->GetXaxis()->GetBinCenter(firstbin);
  double yl = hy->GetXaxis()->GetBinCenter(lastbin);
  py = yl - yf;

  TFitResultPtr rx = hx->Fit("gaus", "SL", "", xmax - px, xmax + px);
  gx = rx->Parameter(2);
  TFitResultPtr ry = hy->Fit("gaus", "SL", "", ymax - py, ymax + py);
  gy = ry->Parameter(2);

  std::cout << "px " << px << " gx " << gx << " py " << py << " gy " << gy << std::endl;

  infile.ReplaceAll(".root", "_out.root");
  TFile f(infile, "recreate");
  t->Fill();
  t->Write();
  f.Close();
  
  if (!batch) {
    TString nid = Form("pos_%2.1f", z);
    prt_canvasAdd(nid, 1600, 800);
    h->Draw("colz");
    prt_canvasAdd(nid + "_x", 1600, 800);
    hx->Draw("colz");
    prt_canvasAdd(nid + "_y", 1600, 800);
    hy->Draw("colz");
    prt_canvasSave("data/fit_sim", 0);
  }
}
