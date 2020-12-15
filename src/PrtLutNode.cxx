#include "PrtLutNode.h"

ClassImp(PrtLutNode)

// -----   Default constructor   -------------------------------------------
PrtLutNode::PrtLutNode(){ 
  fSize = 0;
}

// -----   Standard constructors   -----------------------------------------
PrtLutNode::PrtLutNode(Int_t nodeId):fDetectorId(nodeId),fSize(0){ 
}

void PrtLutNode::AddEntry(Int_t detectorId, TVector3 dir, Double_t pathid, Int_t nrefl, Double_t time, TVector3 hitpos, TVector3 digipos, Double_t weight, TVector3 vertexpos) {
  fDetectorId = detectorId; 
  fNodeArray.push_back(dir);
  fPathIdArray.push_back(pathid);
  fWeightArray.push_back(weight);
  fNRefl.push_back(nrefl);
  fTimeArray.push_back(time);
  fHitPos.push_back(hitpos);
  fVertexPos.push_back(vertexpos);
  fDigiPos = digipos;
  fSize++;
}

void PrtLutNode::AddEntry(Int_t detectorId, TVector3 dir, Double_t pathid, Int_t nrefl, Double_t time, TVector3 hitpos, TVector3 digipos, Double_t weight,
			  TVector3 d1,TVector3 d2,TVector3 d3,TVector3 d4, TVector3 d5,TVector3 d6,TVector3 d7,TVector3 d8, TVector3 vertexpos) {
  fDetectorId = detectorId; 
  fNodeArray.push_back(dir);
  fNodeArrayCs[0].push_back(dir);
  fNodeArrayCs[1].push_back(d1);
  fNodeArrayCs[2].push_back(d2);
  fNodeArrayCs[3].push_back(d3);
  fNodeArrayCs[4].push_back(d4);
  fNodeArrayCs[5].push_back(d5);
  fNodeArrayCs[6].push_back(d6);
  fNodeArrayCs[7].push_back(d7);
  fNodeArrayCs[8].push_back(d8);
		  
  fPathIdArray.push_back(pathid);
  fWeightArray.push_back(weight);
  fNRefl.push_back(nrefl);
  fTimeArray.push_back(time);
  fHitPos.push_back(hitpos);
  fVertexPos.push_back(vertexpos);
  fDigiPos = digipos;
  fSize++;
}
