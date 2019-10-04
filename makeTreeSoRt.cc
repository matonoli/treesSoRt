// Based on main01.cc from the 
// PYTHIA event generator, but has been modified by
// Peter Christiansen (peter.christiansen@hep.lu.se at Lund University)
//
//	Modified by Oliver Matonoha 09/2019
//

#include "Pythia8/Pythia.h"

// ROOT includes
#include <TH1.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TF1.h>
#include <TRandom3.h>

#include "TransverseSpherocity/TransverseSpherocity.h"

using namespace std;
using namespace Pythia8;

const int nStrange = 7;
const int strangePDGs[nStrange] = 
	{310, 313, 333, 3122, 3312,
	3322, 3334};
bool isStrange(Int_t pdg) {
	for (int iS = 0; iS<nStrange; iS++) { 
		if ( TMath::Abs(pdg) == strangePDGs[iS] ) return true;	}
	return false;
}
double DeltaPhi(Double_t phi1, Double_t phi2) {
	
	Double_t dphi = phi2 - phi1;
	if ( dphi > TMath::Pi() )		dphi = dphi - 2*TMath::Pi();
	if ( dphi < -1.*TMath::Pi() )	dphi = dphi + 2*TMath::Pi();

	return dphi;
}

bool IsTrans(Double_t phi1, Double_t phiTrig) {

	Double_t dphi = DeltaPhi(phi1,phiTrig);
	if (TMath::Abs(dphi) < TMath::Pi()/3.) 		return false;
	if (TMath::Abs(dphi) > 2.*TMath::Pi()/3.) 	return false;
	
	return true;
}

int WhatRegion(Double_t phi1, Double_t phiTrig) {

	Double_t dphi = DeltaPhi(phi1,phiTrig);
	
	if (TMath::Abs(dphi) < TMath::Pi()/3.)			return 1;
	else if (TMath::Abs(dphi) > 2.*TMath::Pi()/3.) 	return 2;
	else return 0;
}

int main(int argc, const char **argv) {

	const char *defaults[6] = {"","2000","pytest.root","1"};
	if ( argc < 3 ) {
		argv = defaults;
		argc = 3;
		cout << "Using default arguments..." << endl;
	}

	Int_t nEvents = stoi(argv[1]);// InFileName = argv[0];
	cout << "Number of events: " << nEvents << endl;
	
	TString OutFileName = argv[2];
	cout << "Output file is "<< OutFileName.Data() << endl;

	Bool_t showInfo = (Bool_t)stoi(argv[3]);
	cout << "showInfo is " << showInfo << endl;

	// Set up output file
	TFile * fout = new TFile(OutFileName.Data(), "RECREATE");

	// Analysis parameters
	const Bool_t enhanceRt = false;
	const Int_t minTracks = 10;
	const Float_t cutEta = 0.8;
	const Float_t ptLeadCut = 5.0;
	enum { pi, k, p, k0s, l, phi, xi, partSize};
	int PDGs[partSize] = { 211, 321, 2212, 310, 3122, 333, 3312 };
	TF1* pEffi[partSize];
	pEffi[pi]	= new TF1("pi_pos_Eff", "(x>=0&&x<0.2)*(0.65)+(x>=0.2&&x<[0])*([1]+x*[2])+(x>=[0]&&x<[3])*([4]+[5]*x+[6]*x*x)+(x>=[3])*([7])", 0.0, 20.0);
	pEffi[pi]->SetParameters(0.4,0.559616,0.634754,1.5,0.710963,0.312747,-0.163094,0.813976);
	pEffi[k]	= new TF1("k_pos_Eff", "(x>=0&&x<0.25)*(0.06)+(x>=0.25&&x<[0])*([1]+[2]*x+[3]*x*x+[4]*x*x*x+[5]*x*x*x*x)+(x>=[0]&&x<[6])*([7]+[8]*x)+(x>=[6])*([9])", 0.0, 20.0);
	pEffi[k]->SetParameters(1.0,-0.437094,4.03935,-6.10465,4.41681,-1.21593,3.0,0.694729,0.0238144,0.784834);        
	pEffi[p]	= new TF1("p_pos_Eff", "(x>=0&&x<0.3)*(0.05)+(x>=0.3&&x<[0])*([1]+[2]*x)+(x>=[0]&&x<[3])*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)+(x>=[3])*([9])", 0.0, 20.0);
	pEffi[p]->SetParameters(0.35,-1.49693,5.55626,2.2,0.477826,1.10708,-1.08169,0.419493,-0.0565612,0.802333);
	pEffi[k0s] 	= new TF1("pEffi_k0s","([0]*x^[1]*exp(-x) - [2]*x+[3])*(x>[4])", 0, 20);
	pEffi[k0s]->SetParameters(-7.17967e-01,1.27198e-01,1.51269e-02,5.19441e-01,0.1);
	pEffi[l] 	= new TF1("pEffi_l","([0]*x^[1]*exp(-x) - [2]*x+[3])*(x>[4])", 0, 20);
	pEffi[l]->SetParameters(-5.46924e-01,8.50191e-03,1.49558e-02,3.96383e-01,0.35);
	pEffi[phi]	= new TF1("pEffi_phi","([0]*x^[1]*exp(-x) - [2]*x +[3])*(x>[4])", 0, 20);
	pEffi[phi]->SetParameters(-5.25411e-01,5.42076e-02,-2.99137e-03,3.03124e-01,2.77814e-02); 
	pEffi[xi]	= new TF1("xiEff", "0*(x<[0]) + ([1]*(x-[0])+[2]*(x-[0])*(x-[0]))*([0]<=x&&x<[3]) +[4]*(1-[5]/x)*([3]<=x)", 0.0, 20.0);
	pEffi[xi]->SetParameters(0.643, 0.114, -0.00594, 3.06, 0.283, 0.489);
	for (int iF = 0; iF < partSize; iF++) pEffi[iF]->Write();


  	// Initialize PYTHIA minbias Generator.
	Pythia8::Pythia pythia;
	pythia.readString("Beams:eCM = 13000."); // 7 TeV pp
  	if (enhanceRt) pythia.readString("PhaseSpace:pTHatMin = 4.5");
  	/*
	SoftQCD:all = on                   ! Allow total sigma = elastic/SD/DD/ND
	Optionally only study one or a few processes at a time.
	SoftQCD:elastic = on               ! Elastic
	SoftQCD:singleDiffractive = on     ! Single diffractive
	SoftQCD:doubleDiffractive = on     ! Double diffractive
	SoftQCD:centralDiffractive = on    ! Central diffractive
	SoftQCD:nonDiffractive = on        ! Nondiffractive (inelastic)
	SoftQCD:inelastic = on             ! All inelastic
   	*/
	pythia.readString("SoftQCD:nonDiffractive = on");    
	pythia.readString("SoftQCD:doubleDiffractive = on");   
	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = 0");

	pythia.init();

	// Create histograms and other analysis objects
	enum { gen , genNoPt, rec, recNoPt, TSsize};
	const char* TSnames[TSsize] = { "gen" , "genNoPt", "rec", "recNoPt"};
	TransverseSpherocity* TS[TSsize];
	for (int iTS = 0; iTS < TSsize; iTS++)	{
		TS[iTS] = new TransverseSpherocity();
		TS[iTS]->SetMinMulti(minTracks);
	}
	TRandom3* random = new TRandom3();

	// Create tree and branches
	TTree* tree = new TTree("tree", "PYTHIA Track Tree");
	const Int_t maxSize = 10e3;					// max size of particle arrays / event
	TClonesArray trackArray("TParticle", maxSize);
	tree->Branch("tracks", &trackArray);		// why bronch?
    Float_t evSo[TSsize];
    for (int iTS = 0; iTS < TSsize; iTS++)	{
    	tree->Branch(Form("evSo%s",TSnames[iTS]),&evSo[iTS],
    		Form("evSo%s/f",TSnames[iTS]));
    }
    Float_t evPtLeadgen;
    tree->Branch("evPtLeadgen", &evPtLeadgen, "evPtLeadgen/F");
    Float_t evPhiLeadgen;
    tree->Branch("evPhiLeadgen", &evPhiLeadgen, "evPhiLeadgen/F");
    Float_t evEtaLeadgen;
    tree->Branch("evEtaLeadgen", &evEtaLeadgen, "evEtaLeadgen/F");
    Float_t evPtLeadrec;
    tree->Branch("evPtLeadrec", &evPtLeadrec, "evPtLeadrec/F");
    Float_t evPhiLeadrec;
    tree->Branch("evPhiLeadrec", &evPhiLeadrec, "evPhiLeadrec/F");
    Float_t evEtaLeadrec;
    tree->Branch("evEtaLeadrec", &evEtaLeadrec, "evEtaLeadrec/F");
    Int_t evNchTrans;
    tree->Branch("evNchTrans", &evNchTrans, "evNchTrans/I");
    Int_t evNchTransRec;
    tree->Branch("evNchTransRec", &evNchTransRec, "evNchTransRec/I");

	// Event loop
	int   nRealEvents = 0;
	for (int iEvent = 0; iEvent < nEvents; ++iEvent)	{
	
		int nTr  = 0;
		if (!pythia.next()) continue;
		nRealEvents++;

		trackArray.Clear();
		for (int iTS = 0; iTS < TSsize; iTS++) {
			TS[iTS]->Reset();
			evSo[iTS] = -1;
		}
		evPtLeadgen = -1.; evPhiLeadgen = 0; evEtaLeadgen = 0;
		evPtLeadrec = -1.; evPhiLeadrec = 0; evEtaLeadrec = 0;
		evNchTrans = -1; evNchTransRec = -1;

		Int_t nChargedFinal = 0;
		Int_t nChargedFinalRec = 0;
		Int_t nTransCh = 0;
		Int_t nTransChRec = 0;
		std::vector<Double_t> phis;
		std::vector<Double_t> phisRec;

		// Particle loop
		for (int iP = 0; iP < pythia.event.size(); ++iP)	{
	  
	  		Particle& p = pythia.event[iP];
	  	
	  		Bool_t saveTrack = false;
	  		if (p.isFinal())	saveTrack = true;				//save final-state charged: pi+-, K+-, p, e
																//save final-state neutrals: gamma, K0L, n
	  		if ( isStrange(p.id()) ) 	saveTrack = true;		//save final-state strangeness: phi, K0s, L, Xi, Omega
	  		if (TMath::Abs(p.id()) == 111)	saveTrack = true;	// save also pi0
	  		if (!saveTrack) continue;	

			TParticle* track =	new( trackArray[nTr++] ) TParticle();
			track->SetPdgCode(p.id());
			track->SetFirstMother(p.mother1());
	  		track->SetLastMother(pythia.event[p.mother1()].id());	// storing mother PDG as last mother :-X
	  		track->SetFirstDaughter(p.daughter1());
	  		track->SetLastDaughter(pythia.event[p.daughter1()].id());	// storing daughter PDG as last daughter :-X
			track->SetMomentum(p.px(), p.py(), p.pz(), p.e());
			track->SetProductionVertex(p.xProd(), p.yProd(), p.zProd(), p.tProd());

			// generated
			// calculate spherocities, multiplicities
			Bool_t chargedFinal = p.isFinal() && p.isHadron() && p.isCharged() && (TMath::Abs(p.eta())<cutEta);
			if (chargedFinal) {
				nChargedFinal++;
				TS[gen]->AddTrack(p.px(), p.py());
				if (p.pT()>0) TS[genNoPt]->AddTrack(p.px()/p.pT(), p.py()/p.pT());
			}

			// apply efficiency
			Double_t mcRec = random->Uniform(0.,1.);
			Bool_t isReco = false;
			for (int iPdg = 0; iPdg < 3; iPdg++) {	// for pi,k,p only, non-finals need a different approach
				if (TMath::Abs(PDGs[iPdg])==p.id()) isReco = ( mcRec < pEffi[iPdg]->Eval(p.pT()) );
			}

			// reconstructed
			Bool_t chargedFinalRec = chargedFinal && isReco;
			// calculate spherocities
			if (chargedFinalRec) {
				nChargedFinalRec++;
				TS[rec]->AddTrack(p.px(), p.py());
				if (p.pT()>0) TS[recNoPt]->AddTrack(p.px()/p.pT(), p.py()/p.pT());
			}

			// calculate multiplicities

			// calculate rt generated
			if (chargedFinal) {
				phis.push_back(p.phi());
				if (p.pT() > evPtLeadgen) {
					evPtLeadgen = p.pT();
					evPhiLeadgen = p.phi();
					evEtaLeadgen = p.eta();
				}
			}
			// calculate rt reconstructed
			if (chargedFinalRec) {
				phisRec.push_back(p.phi());
				if (p.pT() > evPtLeadrec) {
					evPtLeadrec = p.pT();
					evPhiLeadrec = p.phi();
					evEtaLeadrec = p.eta();
				}
			}

		}
	
		// Fill other event info
		if (nChargedFinal > minTracks) {
			evSo[gen] = TS[gen]->GetTransverseSpherocityTracks();
			evSo[genNoPt] = TS[genNoPt]->GetTransverseSpherocityTracks();
		}
		if (nChargedFinalRec > minTracks) {
			evSo[rec] = TS[rec]->GetTransverseSpherocityTracks();
			evSo[recNoPt] = TS[recNoPt]->GetTransverseSpherocityTracks();
		}
		if (evPtLeadgen > ptLeadCut) {
			for (auto iPhi : phis) if (IsTrans(iPhi,evPhiLeadgen)) nTransCh++;

				//cout << "ntr " << nTransCh << endl;
				//cout << "ev " << evNchTrans << endl;
			
			evNchTrans = nTransCh;
		}
		if (evPtLeadrec > ptLeadCut) {
			for (auto iPhi : phisRec) if (IsTrans(iPhi,evPhiLeadrec)) nTransChRec++;
			
			evNchTransRec = nTransChRec;
		}
		
		
		tree->Fill();	// Update tree for this event
	
	} // End of event loop.
  
	if (showInfo) pythia.stat();	// write PYTHIA summary to screen

	// Write and close output file
	fout->Write();
	fout->Close();
  
	// Check to see that we got most of the events we wanted
	cout << "Real events/simulated: " << nRealEvents << "/ " << nEvents << endl;
  
	return 0;
}
