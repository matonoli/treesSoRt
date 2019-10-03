#ifndef TRANSVERSESPHEROCITY__H
#define TRANSVERSESPHEROCITY__H

#include "TNamed.h"
#include "TMath.h"
#include "TH1.h"

class TransverseSpherocity : public TNamed {
 public:
  TransverseSpherocity(); //default
  ~TransverseSpherocity();
  
  void Reset() { fNtracks = 0; }
  void AddTrack(Double_t px, Double_t py) { fPx[fNtracks] = px; fPy[fNtracks] = py; fNtracks++; }
  Double_t GetTransverseSpherocity();
  Double_t GetTransverseSpherocityTracks(); //erase if not coreect
  TH1D* GetHistSpher() { return fHistSpher; } 
  Int_t GetMinimizingTrackIndex() { return fMinimizingIndex; };
  void SetMinMulti(Int_t minMulti) { fMinMulti = minMulti; }
  Int_t GetNTracks() { return fNtracks; };
 private:

  Int_t    fMinMulti;
  Int_t    fNtracks;
  Double_t *fPx; //!
  Double_t *fPy; //!
  TH1D     *fHistSpher;
  Int_t fMinimizingIndex;

  ClassDef(TransverseSpherocity, 1);
};

#endif
