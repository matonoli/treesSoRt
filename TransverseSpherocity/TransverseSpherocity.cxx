#include "TransverseSpherocity.h"

ClassImp(TransverseSpherocity)

TransverseSpherocity::TransverseSpherocity():
  fMinMulti(10),
  fNtracks(0),
  fMinimizingIndex(0)
{

  fPx = new Double_t[10000];
  fPy = new Double_t[10000];

  fHistSpher = new TH1D("histSpher","Spher; S_{0}; Counts", 50, 0, 1);
  fHistSpher->Sumw2();
  fHistSpher->SetDirectory(0);
};


//____________________________________________________________________
TransverseSpherocity::~TransverseSpherocity() 
{
  delete [] fPx;
  delete [] fPy;
  delete fHistSpher;
};

//____________________________________________________________________
Double_t TransverseSpherocity::GetTransverseSpherocity() 
{  
  if(fNtracks < fMinMulti) 
    return -1;
  
  Int_t stepSize=0.01;
  Double_t RetTransverseSpherocity = 1000;
  Double_t sumpt = 0;
  for(Int_t i = 0; i < 360/stepSize; ++i) {
    
    //Divide the whole azimuth into segments and do the projection on these segments (below)
    Double_t phiparam = ((TMath::Pi()) * i * stepSize) / 180; 
    Double_t nx = TMath::Cos(phiparam); // x component of a unitary vector n
    Double_t ny = TMath::Sin(phiparam); // y component of a unitary vector n
    
    Double_t num = 0;
    for(Int_t j = 0; j < fNtracks; j++) {
      num += TMath::Abs(ny*fPx[j] - nx*fPy[j]);

      if(i==0)
	sumpt += TMath::Sqrt(fPx[j]*fPx[j] + fPy[j]*fPy[j]);
    }

    Double_t pFull = TMath::Power((num/sumpt), 2); //Projection of sp. on the segment
    if(pFull < RetTransverseSpherocity)  //Select the lowest projection
      RetTransverseSpherocity = pFull;
  };

  RetTransverseSpherocity *= TMath::Pi()*TMath::Pi()/4.0;
  fHistSpher->Fill(RetTransverseSpherocity);
  
  return RetTransverseSpherocity;
};

//____________________________________________________________________
Double_t TransverseSpherocity::GetTransverseSpherocityTracks() 
{  
  if(fNtracks < fMinMulti) 
    return -1;
  
  //Int_t stepSize=2; //do not need this one here right?
  Double_t RetTransverseSpherocity = 1000;
  Double_t sumpt = 0;
  //const Double_t pt = 1;
  for(Int_t i = 0; i < fNtracks; i++) {
    
    Double_t pt = TMath::Sqrt(fPx[i]*fPx[i] + fPy[i]*fPy[i]);
    Double_t nx =fPx[i] / pt; // x component of a unitary vector n
    Double_t ny =fPy[i] / pt; // y component of a unitary vector n
    
    Double_t num = 0;
    for(Int_t j = 0; j < fNtracks; j++) {
      num += TMath::Abs(ny*fPx[j] - nx*fPy[j]);

      if(i==0)
	sumpt += TMath::Sqrt(fPx[j]*fPx[j] + fPy[j]*fPy[j]);
    }

    Double_t pFull = TMath::Power((num/sumpt), 2); //Projection of sp. on the segment
    if(pFull < RetTransverseSpherocity)  { //Select the lowest projection
      RetTransverseSpherocity = pFull;
      fMinimizingIndex = i;
    };
  };

  RetTransverseSpherocity *= TMath::Pi()*TMath::Pi()/4.0;
  fHistSpher->Fill(RetTransverseSpherocity);
  
  return RetTransverseSpherocity;
};

