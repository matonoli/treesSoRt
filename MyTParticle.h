// expanded class of the root TParticle
// OliverM 2019 Lund

#ifndef MYTPARTICLEH
#define MYTPARTICLEH

#include "TParticle.h"
#include "TObject.h"

class MyTParticle: public TParticle {

	public:
		
		MyTParticle() : fIsReco(0), fRegion(-1), fRegionRec(-1) { }	
		~MyTParticle() { }

		Bool_t IsReco() 			{	return fIsReco;};
		//Int_t GetRegion()			{	return fRegion;};
		//Int_t GetRegionRec()		{	return fRegionRec;};

		void SetIsReco(Bool_t b)	{	fIsReco = b;};
		//void SetRegion(Int_t r)		{	fRegion = r;};
		//void SetRegionRec(Int_t r)	{	fRegionRec = r;};
		
		ClassDef(MyTParticle,1);

	protected:

		Bool_t fIsReco;
		//Int_t fRegion;
		//Int_t fRegionRec;

		
};
#endif