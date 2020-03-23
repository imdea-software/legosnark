#ifndef CPHAD_H
#define CPHAD_H

#include "snark.h"
#include "sumcheck.h"
#include "poly.h"

using HadRand = vector<CommRand>;
using HadField = CommRand;


struct HadIn
{
  vector<HadField> a, b, c;
};


struct HadRel
{
  size_t n;

  HadRel(size_t _n) : n(_n) { } ;

  using InT = HadIn;

};



struct HadKey {
	long n;
	long d; // logarithm of n
	CPPoly *cppoly;
	SumcheckKey *scCrs;
	HadKey(long _n, long _d, CPPoly *_cppoly, SumcheckKey *_scKey) :
	 	n(_n), d(_d), cppoly(_cppoly), scCrs(_scKey) {}
};

struct HadPf {
	// randomness r
	HadRand r;

	PolyPf polyProof;
	Comm uProdEvalCm;

	SumcheckPf *sumcheckPf;


	HadPf(long d)
	{
			r.resize(d);
	}

};




// TODO: make this a subclass of CompositeCPSnark
class CPHad : public CPSnark<HadRel, HadKey, HadPf, CPPIn, CPVIn>
{
public:

  using RelT = HadRel;
  using InT = RelT::InT; // XXX: These lines should be moved earlier in the class hierarchy

  CPHad(CommScheme *_commScm, CPPoly *_cppoly) :
    CPSnark(_commScm), cpsumcheck(_commScm, _cppoly), cppoly(_cppoly)
  {
		addBenchmarkSlave(&cpsumcheck, "Sumcheck in Hadamard");
  }
  ~CPHad() {}

  virtual HadKey* keygen(const HadRel *rel) override;
  virtual HadPf* prove(const HadKey *crs, const CPPIn &in) override;
  virtual bool verify(const HadKey *crs, const CPVIn &, const HadPf *pf) override;

	CPSumcheck cpsumcheck;
	CPPoly *cppoly;
protected:


};



#endif
