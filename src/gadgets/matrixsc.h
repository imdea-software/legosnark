#ifndef CPMATSC_H
#define CPMATSC_H

#include "snark.h"
#include "sumcheck.h"
#include "poly.h"

using MatRel = long;
using MatRand = vector<CommRand>;
using MatField = CommRand;

struct MatKey {
  long n;
  long d; // logarithm of n
  CPPoly *cppoly;
  SumcheckKey *scCrs;
  MatKey(long _n, long _d, CPPoly *_cppoly, SumcheckKey *_scKey) :
          n(_n), d(_d), cppoly(_cppoly), scCrs(_scKey) {}
};

struct MatPf {
  // randomness r
  MatRand r, s;

  PolyPf polyProof;
  Comm uProdEvalCm;

  SumcheckPf *sumcheckPf;


  MatPf(long d)
  {
    r.resize(d >> 1);
    s.resize(d >> 1);
  }

  MatRand getRndConcat() const
  {
    auto out = r;
    out.insert(out.end(), s.begin(), s.end());
    return out;
  }

  // size in terms of group elements
  size_t getSize() const
  {
    auto cmSize = 2;
    return sumcheckPf->getSize() + polyProof.getSize() + cmSize;
  }

};

// TODO: make this a subclass of CompositeCPSnark
class CPMat : public CPSnark<MatRel, MatKey, MatPf, CPPIn, CPVIn>
{
public:
  CPMat(CommScheme *_commScm, CPPoly *_cppoly) :
          CPSnark(_commScm), cpsumcheck(_commScm, _cppoly), cppoly(_cppoly)
  {
    addBenchmarkSlave(&cpsumcheck, "Sumcheck in Matrix");
  }

  virtual MatKey* keygen(const MatRel *rel) override;
  virtual MatPf* prove(const MatKey *crs, const CPPIn &in) override;
  virtual bool verify(const MatKey *crs, const CPVIn &, const MatPf *pf) override;

  MatPf* proveOutputMatrixInClear(const MatKey *crs, const CPPIn &in) ;

  bool verifyOutputMatrixInClear(const MatKey *crs, const CPVIn &in, const Ins &outMtx, const MatPf *pf);

  CPSumcheckMatrix cpsumcheck;
  CPPoly *cppoly;
protected:


};

#endif
