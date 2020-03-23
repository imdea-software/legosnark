#include <cstdio>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <utility>
#include <type_traits>
using namespace std;

#include "subspace.h"
#include "util.h"



// Function that interfaces relations for "CPLink" to those of CPSubspace
// XXX: Don't do this at home if you care about efficiency: it is returning/copying something potentially large
SubspaceRel makeLinkingRel(
        const CommScheme &cmScm,
        const vector<LG1> &F)
{
  const size_t nRows = 2;
  const size_t nCols = 2*F.size();

  SubspaceRel ssRel;
  ssRel
    .withNRows(nRows)
    .withNCols(nCols)
    .withoutScalars(); // this line states we do not know the coefficients for the grp elems in M

  vector<ColG1> M(nCols);

  // This corresponds to the h-s in paper
  insertAsColMajor(0, 0, cmScm.getBlindingH(), M);
  insertRowAsColMajor(0, 2, cmScm.getBases1(), M);

  // This corresponds to the f-s in paper
  insertRowAsColMajor(1, 1, F, M);

  ssRel.withMatrix(M);
  return ssRel;
}

// utility functions for random vectors
void init_as_random(vector <LFr> &dst)
{
  for (auto i = 0; i < dst.size(); i++) {
    dst[i] = LFr::random_element();
  }
}

void init_as_random(vector <LG1> &dst)
{
  vector<LFr> tmp(dst.size());
  init_as_random(tmp);
  for (auto i = 0; i < dst.size(); i++) {
    dst[i] = tmp[i]*LG1::one();
  }
}

// wrapper functions to simplify interface to commitment in this ctx
void justCommit(CommScheme &cmScm, const Ins &u, LG1 &cm, LFr &opn)
{
  auto cmOut = cmScm.commit(u);
  cm = cmOut.c.c;
  opn = cmOut.r;
}

void justCommit(const vector<LG1> &bases, const Ins &u, LG1 &cm, LFr &opn)
{
  opn = LFr::random_element();

  auto hiding_base = bases[0];
  vector<LG1> bases_no_hiding(bases.begin()+1, bases.end());
  cm = multiExpMA<LG1>(bases_no_hiding, u) + opn*hiding_base;
}



int main(int argc, char **argv)
{
  default_ec_pp::init_public_params();

  // example input size
  const long N = 1 << 10;

  // commitment scheme (corresponding to the h-s in paper)
  CommScheme cmScm;
  cmScm.keygen(N);

  // random generators
  vector<LG1> F(N+1);
  init_as_random(F);

  SubspaceRel ssRel = makeLinkingRel(cmScm, F);

  // Let's prove cH and cF commit to the same input u
  vector<LFr> u(N);
  init_as_random(u);
  LG1 cH, cF;
  LFr rH, rF; // respective openings

  // Prepare inputs
  justCommit(cmScm, u, cH, rH);
  justCommit(F, u, cF, rF);

  // w, i.e. witness
  vector<LFr> w({rH, rF});
  w.insert(w.end(), u.begin(), u.end());

  // Keygen/Prove/Verify link
  SubspaceSnark ss;
  auto crs = ss.keygen(&ssRel);
  auto pf = ss.prove(crs, w);
  MYREQUIRE(ss.verify(crs, {cH,cF}, pf));

  return 0;
}