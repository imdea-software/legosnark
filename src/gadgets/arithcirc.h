#ifndef ARITH_CIRC_H
#define ARITH_CIRC_H

#include "snark.h"
#include "lipmaa.h"
#include "matrix.h"
#include "subspace.h"

using RCPairs = vector<pair<size_t, size_t>>;

using InT = IScalar;


// generic input structure of arithmetic circuit

struct ACPIn {
	CommOut aCOut, bCOut, cCOut, uCOut;
};

struct ACVIn {
	Comm aC, bC, cC, uC;
};




class ACRel
{
public:
	vector<InT> Wtrs;
	Ins tgtV;
	long n_mul_gates, r, c;
    size_t input_size;
    bool has_tgtV;
	
	ACRel(size_t n, size_t _r, size_t _c, size_t u_sz, const vector<InT> &_Wtrs, const Ins &_tgtV) : 
		n_mul_gates(n), r(_r), c(_c), input_size(u_sz), Wtrs(_Wtrs), tgtV(_tgtV), has_tgtV(_tgtV.size() != 0) {
			
	}

};

struct ACKey {
	HadLKey *hadkey;
	SubspaceKey *sskey;
	vector<LG1> g1TgtV;
	Fqk<def_ec> tgtV_C_precomp;
 	bool has_tgtV;
    size_t input_size;

};

struct ACPf {
	HadLPf hadpf;
	SubspacePf *sspf;
};

class CPAC : public CPSnark<ACRel, ACKey, ACPf, ACPIn, ACVIn>
{
	Interpolator *interp = nullptr;
public:
  CPHadL cphadl;
  SubspaceSnark ss;

  
  void prepareInputs(const Ins &a, const Ins &b, const Ins &c, ACPIn &pIn, ACVIn &vIn)
  {
	  auto ics = getCommScheme();
	  
	  
	  pIn.aCOut = ics->commit(a);
	  applyBenchmarkFrom(*ics, "commit", "commit_a");  // HACK: We multiply by 3 later
	  pIn.bCOut = ics->commit(b);
	  pIn.cCOut = ics->commit(c);
	  
	  vIn.aC = pIn.aCOut.c;
	  vIn.bC = pIn.bCOut.c;
	  vIn.cC = pIn.cCOut.c;
	  
	  
  }
  
  
  void prepareInputs(const Ins &a, const Ins &b, const Ins &c, const Ins &u, ACPIn &pIn, ACVIn &vIn)
  {
	  auto ics = getCommScheme();
	  
	  
	pIn.aCOut = ics->commit(a);
	applyBenchmarkFrom(*ics, "commit", "commit_a"); // HACK: We multiply by 3 later
	pIn.bCOut = ics->commit(b);
	pIn.cCOut = ics->commit(c);
    pIn.uCOut = ics->commit(u);
    applyBenchmarkFrom(*ics, "commit", "commit_u");

	  
	vIn.aC = pIn.aCOut.c;
	vIn.bC = pIn.bCOut.c;
	vIn.cC = pIn.cCOut.c;
    vIn.uC = pIn.uCOut.c;

	  
	  
  }
	
  CPAC(CommScheme *_commScm) :
    CPSnark( _commScm)
  {
	addBenchmarkSlave(&cphadl, "CPAC's Hadamard SNARK");
	addBenchmarkSlave(getCommScheme(), "CPAC's Commitment Scheme");

  }
  
  virtual ACKey* keygen(const ACRel *rel) override;
  virtual ACPf* prove(const ACKey *crs, const ACPIn &in) override;
  virtual bool verify(const ACKey *crs, const ACVIn &, const ACPf *pf) override;
};


#endif
