// Pythia tree reading macro
#include <TChain.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TH1.h>
#include <TFile.h>
#include <TROOT.h>
#include <TParticlePDG.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>

using namespace std;

// GLOBALS
TChain* mChain;

bool MakeChain(const Char_t *inputFile="test.list") {

	if (!mChain) mChain = new TChain("tree");
	TString inputFileStr(inputFile);

	string const dirFile = inputFileStr.Data();
	if (dirFile.find(".lis") != string::npos)	{
		
		ifstream inputStream(dirFile.c_str());
		if (!inputStream)	{
			cout << "ERROR: Cannot open list file " << dirFile << endl;
			return false;	}

		int nFile = 0;
		string file;
		while (getline(inputStream, file))	{
	  		if (file.find(".root") != string::npos)	{
				TFile* ftmp = TFile::Open(file.c_str());
				if (ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())	{
		  			cout << " Read in tree file " << file << endl;
		  			mChain->Add(file.c_str());
		  			++nFile;	}
				if (ftmp) ftmp->Close();
	  		}
		}

	cout << " Total " << nFile << " files have been read in. " << endl;
	return true;
	}

	else if (dirFile.find(".root") != string::npos)	{
		mChain->Add(dirFile.c_str());	
		return true;	}
  	else	{
		cout << " No good input file to read ... " << endl;
		return false;	}
}

double DeltaPhi(Double_t phi1, Double_t phi2) {
	
	Double_t dphi = phi2 - phi1;
	if ( dphi > TMath::Pi() )		dphi = dphi - 2*TMath::Pi();
	if ( dphi < -1.*TMath::Pi() )	dphi = dphi + 2*TMath::Pi();

	return dphi;
}

void readTree(Int_t nEvents=100, const Char_t *inputFile="test.list", const Char_t *outputFile="test.root") {

	gROOT->ProcessLine(".x load_libraries.C");

	if (!MakeChain(inputFile)) printf("Couldn't create the chain! \n");
	else cout << "Chain created with entries: " << mChain->GetEntries() << "\n";

	// Set up output file
	TFile * fout = new TFile(outputFile, "RECREATE");

	TClonesArray* tracks = 0;
	mChain->SetBranchAddress("tracks", &tracks);


	enum { gen , genNoPt, rec, recNoPt, TSsize};
	const char* TSnames[TSsize] = { "gen" , "genNoPt", "rec", "recNoPt"};
	Float_t evSo[TSsize];
	for (int iTS = 0; iTS < TSsize; iTS++)	{
    	mChain->SetBranchAddress(Form("evSo%s",TSnames[iTS]),&evSo[iTS]);
    }

	Float_t evPtLeadgen;
    mChain->SetBranchAddress("evPtLeadgen", &evPtLeadgen);
    Float_t evPhiLeadgen;
    mChain->SetBranchAddress("evPhiLeadgen", &evPhiLeadgen);
    Float_t evEtaLeadgen;
    mChain->SetBranchAddress("evEtaLeadgen", &evEtaLeadgen);
    Float_t evPtLeadrec;
    mChain->SetBranchAddress("evPtLeadrec", &evPtLeadrec);
    Float_t evPhiLeadrec;
    mChain->SetBranchAddress("evPhiLeadrec", &evPhiLeadrec);
    Float_t evEtaLeadrec;
    mChain->SetBranchAddress("evEtaLeadrec", &evEtaLeadrec);
    Int_t evNchTrans;
    mChain->SetBranchAddress("evNchTrans", &evNchTrans);
    Int_t evNchTransRec;
    mChain->SetBranchAddress("evNchTransRec", &evNchTransRec);
    Int_t evNchCL;
    mChain->SetBranchAddress("evNchCL", &evNchCL);
    Int_t evNchCLRec;
    mChain->SetBranchAddress("evNchCLRec", &evNchCLRec);
    Int_t evNchV0M;
    mChain->SetBranchAddress("evNchV0M", &evNchV0M);

    TH1D* hEvSogen			= new TH1D("hEvSogen","",1000,0.,1.);
    TH1D* hEvSogenNoPt		= new TH1D("hEvSogenNoPt","",1000,0.,1.);
    TH1D* hEvNchV0M			= new TH1D("hEvNchV0M","",100, -1, 99);
    TH1D* hTrackPt			= new TH1D("hTrackPt","",500, -1, 25);

    int nSoCuts = 6;
    TH1D* hDPhiSogen[nSoCuts-1]; TH1D* hDPhiSogenNoPt[nSoCuts-1];
    TH1D* hDPhiSogenNeutral[nSoCuts-1]; TH1D* hDPhiSogenNoPtNeutral[nSoCuts-1];
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
    	hDPhiSogen[iH] = new TH1D(Form("hDPhiSogen_%i",iH),"",200, -3.2, 3.2);
    	hDPhiSogenNoPt[iH] = new TH1D(Form("hDPhiSogenNoPt_%i",iH),"",200, -3.2, 3.2);
    	hDPhiSogenNeutral[iH] = new TH1D(Form("hDPhiSogenNeutral_%i",iH),"",200, -3.2, 3.2);
    	hDPhiSogenNoPtNeutral[iH] = new TH1D(Form("hDPhiSogenNoPtNeutral_%i",iH),"",200, -3.2, 3.2);
    }

    // Calculate spherocity quantiles
    mChain->Draw("evSogen>>hEvSogen","evSogen>0.","goff");
    mChain->Draw("evSogenNoPt>>hEvSogenNoPt","evSogenNoPt>0.","goff");
    double quantileValues[] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
    double cutSogen[nSoCuts]; double cutSogenNoPt[nSoCuts];
    hEvSogen->GetQuantiles(nSoCuts, cutSogen, quantileValues);
    hEvSogenNoPt->GetQuantiles(nSoCuts, cutSogenNoPt, quantileValues);
    for (int i = 0; i < 6; ++i)
    {
    	printf("values for %i are %f and %f \n", i, cutSogen[i], cutSogenNoPt[i]);
    }

	nEvents = (nEvents < mChain->GetEntries() && nEvents > 0) ? nEvents : mChain->GetEntries();
	for (int iEv = 0; iEv < nEvents; ++iEv)	{

		if (iEv%10000==0) printf("Processing: %i out of total %i events...\n", iEv, nEvents);
		mChain->GetEntry(iEv);

		hEvNchV0M->Fill(evNchV0M);

		if (evSo[gen] < 0 && evSo[genNoPt] < 0) continue;

		Int_t nTracks = tracks->GetEntriesFast();
		for (int iTr = 0; iTr < nTracks; ++iTr)	{
			
			TParticle* t = (TParticle*)tracks->At(iTr);
			hTrackPt->Fill(t->Pt());

			if (TMath::Abs(t->Eta()) > 0.8) continue;

			bool isCharged1 = TMath::Abs(t->GetPDG()->Charge()) > 0.001;
			for (int iTr2 = iTr+1; iTr2 < nTracks; ++iTr2)	{
				
				TParticle* t2 = (TParticle*)tracks->At(iTr2);
				bool isCharged2 = TMath::Abs(t2->GetPDG()->Charge()) > 0.001;

				Double_t dPhi = DeltaPhi(t->Phi(), t2->Phi());
				
				for (int iSC = 0; iSC < nSoCuts-1; ++iSC)	{
					
					if (evSo[gen] > cutSogen[iSC] && evSo[gen] < cutSogen[iSC+1]
						&& isCharged1 && isCharged2)
						hDPhiSogen[iSC]->Fill(dPhi);

					if (evSo[genNoPt] > cutSogenNoPt[iSC] && evSo[genNoPt] < cutSogenNoPt[iSC+1]
						&& isCharged1 && isCharged2)
						hDPhiSogenNoPt[iSC]->Fill(dPhi);

					if (evSo[gen] > cutSogen[iSC] && evSo[gen] < cutSogen[iSC+1]
						&& !isCharged1 && !isCharged2)
						hDPhiSogenNeutral[iSC]->Fill(dPhi);

					if (evSo[genNoPt] > cutSogenNoPt[iSC] && evSo[genNoPt] < cutSogenNoPt[iSC+1]
						&& !isCharged1 && !isCharged2)
						hDPhiSogenNoPtNeutral[iSC]->Fill(dPhi);
				}
			}
		}

	}

	for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		hDPhiSogen[iH]->Scale(1./hDPhiSogen[iH]->Integral());
    	hDPhiSogenNoPt[iH]->Scale(1./hDPhiSogenNoPt[iH]->Integral());
    	hDPhiSogenNeutral[iH]->Scale(1./hDPhiSogenNeutral[iH]->Integral());
    	hDPhiSogenNoPtNeutral[iH]->Scale(1./hDPhiSogenNoPtNeutral[iH]->Integral());
    }

    TCanvas* cSo = new TCanvas("cSo","",900, 900);
    cSo->Divide(2,2,1e-04,1e-04);
    cSo->cd(1);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		hDPhiSogen[iH]->SetLineColor(2+iH);
		if (!iH) hDPhiSogen[iH]->Draw();
		else hDPhiSogen[iH]->Draw("same");
	}
	cSo->cd(2);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		hDPhiSogenNoPt[iH]->SetLineColor(2+iH);
		if (!iH) hDPhiSogenNoPt[iH]->Draw();
		else hDPhiSogenNoPt[iH]->Draw("same");
	}
	cSo->cd(3);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		hDPhiSogenNeutral[iH]->SetLineColor(2+iH);
		if (!iH) hDPhiSogenNeutral[iH]->Draw();
		else hDPhiSogenNeutral[iH]->Draw("same");
	}
	cSo->cd(4);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		hDPhiSogenNoPtNeutral[iH]->SetLineColor(2+iH);
		if (!iH) hDPhiSogenNoPtNeutral[iH]->Draw();
		else hDPhiSogenNoPtNeutral[iH]->Draw("same");
	}


	fout->Write();
}
