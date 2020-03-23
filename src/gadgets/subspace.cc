#include "subspace.h"
#include "dbgutil.h"
#include "fmt/format.h"
#include "util.h"
#include "sparsemexp.h"

#include <libff/common/default_types/ec_pp.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
using namespace libff;

#include <set>
using namespace std;
using namespace bn;



// NB: m is col-major
void mtxmultiexp(vector<LG1> &out, const vector<LFr> &exps, const vector<ColG1> &m)
{
	out.resize(m.size());
	int i = 0;
	for (const ColG1 &c : m) {
		out[i++] = simplesparsemexp(c, exps);
	}
}

void mtxmultiexp(vector<LG1> &out, const vector<LFr> &exps, const vector<ColFr> &m, Interpolator *interp)
{
	vector<LFr> tmp(m.size());
	int i = 0;
	for (const ColFr &c : m) {
		tmp[i++] = sparseinnerproduct(exps, c);
	}
	out = interp->mkG1Exp(tmp);
}

SubspaceKey* SubspaceSnark::keygen(const SubspaceRel *rel)
{

  SubspaceKey *key = new SubspaceKey;
  // sample vector k of size of rows
  vector<IScalar> k; 
  sampleRandomVector(k, rel->l); 
  
  cpdbg::print(k, "k", "keygen");
 
  key->a = LFr::random_element()*LG2::one();
  key->a_precomp = def_ec::precompute_G2(key->a);
  cpdbg::print(key->a, "a", "keygen");

  if (rel->scalarsAvailable) {
    // this is more efficient but requires to know exponents, so not always possible
    mtxmultiexp(key->P, k, rel->sM, rel->interp);
  } else {
    mtxmultiexp(key->P, k, rel->M);
  }
  
  assert(key->P.size() == rel->t);
  cpdbg::print(key->P, "P", "keygen");

  // C = a * k (C is a vector, k is a vector, a is Ec2 point)
  key->C = cputil::simpleBatchExp<LG2, LFr>(key->a, k);
  cpdbg::print(key->C, "C", "keygen");
  
  // for CPLin3or4
  if (rel->C_precomp_sz != 0) {
    key->C_precomp.resize(rel->C_precomp_sz);
    for (auto i = 0; i < rel->C_precomp_sz; i++) {
      key->C_precomp[i] = def_ec::precompute_G2(key->C[i]);
    }
  }

  key->rel = rel;

  return key;
}

SubspacePf* SubspaceSnark::prove(const SubspaceKey *crs, const vector<LFr> &w)
{
  SubspacePf *pf = new SubspacePf;
  startBenchmark("prove");
  *pf = multiExpMA<LG1>(crs->P, w);
  stopBenchmark("prove");
  return pf;
}


inline LGT e(const LG1 &x, const LG2 &y)
{
  return def_ec::reduced_pairing(x, y); // XXX: should be optimized with miller_loop
}

LGT eVecFromSet(const vector<LG1> &xs, const vector<LG2> &ys, const set<int> &s)
{
  LGT res = LGT::one();
  for (auto it = s.begin(); it != s.end(); it++) {
    auto x = xs[*it];
    auto y = ys[*it];
    res = res * e(x, y);
  }
  return res;
}



bool SubspaceSnark::verify(const SubspaceKey *crs, const vector<LG1> &xStdVec, const SubspacePf *pf)
{

  cpdbg::print(xStdVec, "xStdVec", "subspace_verify");

  // optimize by ignoring non zero elements
  set<int> nonZeroSet;
  for (auto i = 0; i < xStdVec.size(); i++) {
    if (!xStdVec[i].is_zero()) {
      nonZeroSet.insert(i);
    }
  }


  startBenchmark("verify");
  //EcTgt lhs = AlgPP::eVec(xStdVec, crs->C);
  // skip non zero elements for pairing
  auto lhs = eVecFromSet(xStdVec, crs->C, nonZeroSet);
  auto rhs = e(*pf, crs->a);
  stopBenchmark("verify");

  cpdbg::print(lhs, "lhs", "subspace_verify");
  cpdbg::print(*pf, "*pf", "subspace_verify");
  cpdbg::print(crs->a, "a", "subspace_verify");
  cpdbg::print(rhs, "rhs", "subspace_verify");

  return (lhs == rhs);
}

bool SubspaceSnark::verifyLin3or4(
	const SubspaceKey *crs, 
	const vector<LG1> &cs,
	const Fqk<def_ec> *aux_precomp,
	const SubspacePf *pf)
{
  const unsigned nComms = cs.size();
  vector<G1_precomp<def_ec>> com_precomp(nComms);
  for (auto i = 0; i < nComms; i++) {
	  com_precomp[i] = def_ec::precompute_G1(cs[i]);
  }


  Fqk<def_ec> lhs;
  
  // pairings on commitment and C[0..2]
  if (nComms==3) {
    auto eComm0 = def_ec::miller_loop(com_precomp[0], crs->C_precomp[0]);
    auto eComm12 = def_ec::double_miller_loop(com_precomp[1], crs->C_precomp[1], com_precomp[2], crs->C_precomp[2]);
    lhs = eComm0*eComm12;
  }
  if (nComms==4) { // pairings on commitment and C[0..3] 
    auto eComm01 = def_ec::double_miller_loop(com_precomp[0], crs->C_precomp[0], com_precomp[1], crs->C_precomp[1]);
    auto eComm23 = def_ec::double_miller_loop(com_precomp[2], crs->C_precomp[2], com_precomp[3], crs->C_precomp[3]);
    lhs = eComm01*eComm23;
  }
  if (aux_precomp != nullptr){
    lhs = lhs* *aux_precomp;
  }
  
  auto pf_precomp = def_ec::precompute_G1(*pf);
  auto rhs = def_ec::miller_loop(pf_precomp, crs->a_precomp);
  
  auto shouldBeOne = def_ec::final_exponentiation(lhs*rhs.unitary_inverse());
  bool isGood = (shouldBeOne == LGT::one());
  return isGood;
}
