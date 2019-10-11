// Pythia tree reading macro
#include <TChain.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TH1.h>
#include <TFile.h>
#include <TROOT.h>
#include <TParticlePDG.h>
#include <TCanvas.h>
#include <TLegend.h>

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

void MakeNiceHistogram(TH1D* h, Int_t col) {

	h->SetLineWidth(2);
	h->SetLineColor(col);
	h->SetMarkerStyle(20);
	h->SetMarkerSize(1.3);
	h->SetMarkerColor(col);
	h->SetStats(0);

	h->GetYaxis()->SetTitleOffset(1.1);
	h->GetYaxis()->SetLabelOffset(0.0025);
	h->GetYaxis()->SetLabelSize(0.03);

	//h->SetTopMargin(0.055);
}

void MakeNiceLegend(TLegend *leg, Float_t size, Int_t columns)	{
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetMargin(0.25);
	leg->SetTextSize(size);
	leg->SetEntrySeparation(0.5);
	leg->SetNColumns(columns);

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

    TH1D* hEvSo[TSsize];
    for (int iTS = 0; iTS < TSsize; ++iTS)	{
    	hEvSo[iTS]	= new TH1D(Form("hEvSo_%s",TSnames[iTS]),"",1000,0.,1.);
    }
	
	TH1D* hEvNchV0M			= new TH1D("hEvNchV0M","",100, -1, 99);
    TH1D* hTrackPt			= new TH1D("hTrackPt","",500, -1, 25);

    int nSoCuts = 6;
    Int_t clrs[nSoCuts-1] = { kRed, kBlue, kGreen+2, kOrange+10, kMagenta};
    TH1D* hDPhiSo[TSsize][nSoCuts-1];
    TH1D* hDPhiSoNeutral[TSsize][nSoCuts-1];
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
    for (int iTS = 0; iTS < TSsize; ++iTS)	{
    	hDPhiSo[iTS][iH]		= new TH1D(Form("hDPhiSo_%s_%i",TSnames[iTS],iH),"",100, -3.2, 3.2);
    	hDPhiSoNeutral[iTS][iH]	= new TH1D(Form("hDPhiSoNeutral_%s_%i",TSnames[iTS],iH),"",100, -3.2, 3.2);
    }	}


    // Calculate spherocity quantiles
    double quantileValues[] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
    double cutSo[TSsize][nSoCuts];
    for (int iTS = 0; iTS < TSsize; ++iTS)	{
    	mChain->Draw(Form("evSo%s>>hEvSo_%s",TSnames[iTS],TSnames[iTS]),Form("evSo%s>0.",TSnames[iTS]),"goff");
    	hEvSo[iTS]->GetQuantiles(nSoCuts, cutSo[iTS], quantileValues);
    }

	nEvents = (nEvents < mChain->GetEntries() && nEvents > 0) ? nEvents : mChain->GetEntries();
	for (int iEv = 0; iEv < nEvents; ++iEv)	{

		if (iEv%10000==0) printf("Processing: %i out of total %i events...\n", iEv, nEvents);
		mChain->GetEntry(iEv);

		hEvNchV0M->Fill(evNchV0M);

		if (evNchCLRec < 27) continue;
		if (evSo[rec] < 0 && evSo[recNoPt] < 0) continue;
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
				for (int iTS = 0; iTS < TSsize; ++iTS)	{
					
					if (evSo[iTS] > cutSo[iTS][iSC] && evSo[iTS] < cutSo[iTS][iSC+1]
						&& isCharged1 && isCharged2)
						hDPhiSo[iTS][iSC]->Fill(dPhi);

					if (evSo[iTS] > cutSo[iTS][iSC] && evSo[iTS] < cutSo[iTS][iSC+1]
						&& !isCharged1 && !isCharged2 && t->GetPdgCode()!=22 && t2->GetPdgCode()!=22)
						hDPhiSoNeutral[iTS][iSC]->Fill(dPhi);
				}	}
			}
		}

	}

	for (int iH = 0; iH < nSoCuts-1; ++iH)	{
	for (int iTS = 0; iTS < TSsize; ++iTS)	{
		hDPhiSo[iTS][iH]->Scale(1./hDPhiSo[iTS][iH]->Integral());
    	hDPhiSoNeutral[iTS][iH]->Scale(1./hDPhiSoNeutral[iTS][iH]->Integral());
    }	}

    TLegend* l1 = new TLegend(0.15,0.75,0.75,0.88);
    MakeNiceLegend(l1,0.037,2);
    TCanvas* cSo_gen = new TCanvas("cSo_gen","",900, 900);
    cSo_gen->Divide(2,2,1e-04,1e-04);
    cSo_gen->cd(1);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		MakeNiceHistogram(hDPhiSo[gen][iH],clrs[iH]);
		if (!iH) {	
			hDPhiSo[gen][iH]->Draw();
			hDPhiSo[gen][iH]->GetYaxis()->SetRangeUser(0.0094,0.0110);
			hDPhiSo[gen][iH]->GetYaxis()->SetTitle("a.u.");	
		}
		else hDPhiSo[gen][iH]->Draw("same");
	}

	cSo_gen->cd(2);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		MakeNiceHistogram(hDPhiSo[genNoPt][iH],clrs[iH]);
		if (!iH) {	
			hDPhiSo[genNoPt][iH]->Draw();
			hDPhiSo[genNoPt][iH]->GetYaxis()->SetRangeUser(0.0094,0.0110);
			hDPhiSo[genNoPt][iH]->GetYaxis()->SetTitle("a.u.");	
		}
		else hDPhiSo[genNoPt][iH]->Draw("same");
	}
	cSo_gen->cd(3);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		MakeNiceHistogram(hDPhiSoNeutral[gen][iH],clrs[iH]);
		if (!iH) {	
			hDPhiSoNeutral[gen][iH]->Draw();
			hDPhiSoNeutral[gen][iH]->GetYaxis()->SetRangeUser(0.0094,0.0110);
			hDPhiSoNeutral[gen][iH]->GetYaxis()->SetTitle("a.u.");	
		}
		else hDPhiSoNeutral[gen][iH]->Draw("same");
	}
	cSo_gen->cd(4);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		MakeNiceHistogram(hDPhiSoNeutral[genNoPt][iH],clrs[iH]);
		if (!iH) {	
			hDPhiSoNeutral[genNoPt][iH]->Draw();
			hDPhiSoNeutral[genNoPt][iH]->GetYaxis()->SetRangeUser(0.0094,0.0110);
			hDPhiSoNeutral[genNoPt][iH]->GetYaxis()->SetTitle("a.u.");	
			}
		else hDPhiSoNeutral[genNoPt][iH]->Draw("same");
		//l1->AddEntry(hDPhiSoNeutral[genNoPt][iH],Form("%2.0f - %2.0f", quantileValues[iH],quantileValues[iH+1]),"l");
	}
	//l1->Draw();

	TCanvas* cSo_rec = new TCanvas("cSo_rec","",900, 900);
    cSo_rec->Divide(2,2,1e-04,1e-04);
    cSo_rec->cd(1);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		MakeNiceHistogram(hDPhiSo[rec][iH],clrs[iH]);
		if (!iH) {	
			hDPhiSo[rec][iH]->Draw();
			hDPhiSo[rec][iH]->GetYaxis()->SetRangeUser(0.0094,0.0110);
			hDPhiSo[rec][iH]->GetYaxis()->SetTitle("a.u.");	}
		else hDPhiSo[rec][iH]->Draw("same");
	}
	cSo_rec->cd(2);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		MakeNiceHistogram(hDPhiSo[recNoPt][iH],clrs[iH]);
		if (!iH) {	
			hDPhiSo[recNoPt][iH]->Draw();
			hDPhiSo[recNoPt][iH]->GetYaxis()->SetRangeUser(0.0094,0.0110);
			hDPhiSo[recNoPt][iH]->GetYaxis()->SetTitle("a.u.");	}
		else hDPhiSo[recNoPt][iH]->Draw("same");
	}
	cSo_rec->cd(3);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		MakeNiceHistogram(hDPhiSoNeutral[rec][iH],clrs[iH]);
		if (!iH) {	
			hDPhiSoNeutral[rec][iH]->Draw();
			hDPhiSoNeutral[rec][iH]->GetYaxis()->SetRangeUser(0.0094,0.0110);
			hDPhiSoNeutral[rec][iH]->GetYaxis()->SetTitle("a.u.");	}
		else hDPhiSoNeutral[rec][iH]->Draw("same");
	}
	cSo_rec->cd(4);
    for (int iH = 0; iH < nSoCuts-1; ++iH)	{
		MakeNiceHistogram(hDPhiSoNeutral[recNoPt][iH],clrs[iH]);
		if (!iH) {	
			hDPhiSoNeutral[recNoPt][iH]->Draw();
			hDPhiSoNeutral[recNoPt][iH]->GetYaxis()->SetRangeUser(0.0094,0.0110);
			hDPhiSoNeutral[recNoPt][iH]->GetYaxis()->SetTitle("a.u.");	}
		else hDPhiSoNeutral[recNoPt][iH]->Draw("same");
	}


	fout->Write();
}
