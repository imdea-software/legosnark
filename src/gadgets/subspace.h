#ifndef GADGETS_SUBSPACE_H
#define GADGETS_SUBSPACE_H

#include "matrix.h"
#include "snark.h"
#include "commit.h"
#include "interp.h"

struct SubspaceRel {
  vector<ColG1> M;
  vector<ColFr> sM;

  int l; // # of rows
  int t; // # of cols
  
  size_t C_precomp_sz = 0;
  Interpolator *interp;

  bool scalarsAvailable = true; // for efficient keygen when possible

  SubspaceRel() {}

  SubspaceRel &withNRows(int _l) {
    l = _l;
    return *this;
  }
  SubspaceRel &withNCols(int _t) {
    t = _t;
    return *this;
  }

  SubspaceRel &withMatrix(const vector<ColG1> &_M) {
    M = _M;
    return *this;
  }

  SubspaceRel &withoutScalars() {
    scalarsAvailable = false;
    return *this;
  }

  SubspaceRel(int _l, int _t, vector<ColFr> _sM, const vector<ColG1> &_M, size_t _C_precomp_sz, Interpolator *_interp) :
	l(_l), t(_t), sM(_sM), M(_M), interp(_interp), C_precomp_sz(_C_precomp_sz)
  {

  }

};

struct SubspaceKey {

  const SubspaceRel *rel;

  vector<LG1> P;
  
  vector<LG2> C;
  vector<G2_precomp<def_ec>> C_precomp;
  
  LG2 a;
  G2_precomp<def_ec> a_precomp;
  
  void print_size() {
    fmt::print("Size of subspace's PK: {} G1 \n", P.size() );
    fmt::print("Size of subspace's VK: {} G2 \n", 4+1 ); // XXX

  }

};

typedef LG1 SubspacePf;

// This class proves "x = Mw" for witness w
class SubspaceSnark :
  public Snark<SubspaceRel, SubspaceKey, SubspacePf, vector<LFr>, vector<LG1>>
{
public:
  SubspaceSnark() :
    Snark<SubspaceRel, SubspaceKey, SubspacePf, vector<LFr>, vector<LG1>>()
  {

  }

  virtual SubspaceKey* keygen(const SubspaceRel *rel) override;
  // NB: In method prove the ScalarVec is a col vector
  virtual SubspacePf* prove(const SubspaceKey *crs, const vector<LFr> &) override;
  virtual bool verify(const SubspaceKey *crs, const vector<LG1> &x, const SubspacePf *pf) override;
  bool verifyLin3or4(const SubspaceKey *crs, const vector<LG1> &cs, const Fqk<def_ec> *aux_precomp, const SubspacePf *pf);
};


#endif
