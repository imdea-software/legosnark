#include "sumcheck.h"
#include <utility>
using std::swap;

SumcheckKey* CPSumcheck::keygen(const SumcheckRel *rel)
{
  SumcheckKey *crs = new SumcheckKey(*rel);
  return crs;
}

SumcheckPf *CPSumcheck::prove(const SumcheckKey *crs, const CPPIn &in)
{
  /*
   * Input format:
   * publicSlot: rho (fixed randomness in beta)
   * comSlot:
   *  y: alleged result
   *  a: first vector
   *  b: second vector
   *
   */

  auto c = CommOut::toComms(in.commSlot);
  auto u = CommOut::toCommittedVals(in.commSlot);

  SumcheckField y = u[0][0]; // alleged result
  auto &yComOut = in.commSlot[0];

  const Ins &a = u[1];
  const Ins &b = u[2];

  const Ins &rho = in.publicSlot;
  assert(rho.size() == crs->d);

  startBenchmark("prove");

  // initialize dynamic programming engine for poly optimizations
  shared_ptr<DPBeta> beta = init_beta(crs->d, rho);
  // for i in 0 to n_cm_polys-1...
  vector<shared_ptr<DPMle>> mles;
  init_mles(mles, crs->d, rho, a, b);

  // generate randomnesses
  SumcheckRand r(crs->d);

  for(auto i=0; i < crs->d; i++) {
      r[i] = CommRand::random_element();
	}

  vector<PolyT> h(crs->d); // vector of univariate polynomial

  Scalars z(crs->d+1); // vector of results
  z[0] = y;
  for (auto i = 0; i < crs->d; i++) {
    // h_i(x_i) := "sum of g(r0, ..., r(i-1), x_i, ..., x(d-1))"
    h[i] = make_new_h_poly(crs->d, i, beta, mles);

    // next point
    z[i+1] = h[i].eval(r[i]);

    // update engine (at all iterations except last)
    if (i < crs->d - 1) {
      beta->pushRandomness(r[i], i);
      for (shared_ptr<DPMle> mle : mles) {
        mle->pushRandomness(r[i], i);
      }
    }
  }

  /* Below we make hCom and zk-equality proofs */

  // commit to h[i]
  auto comScm = getCommScheme();
  // make CommOuts
  auto mkComOut = [&comScm](const PolyT &p) { return p.commit(comScm); };
  auto hComOut = cputil::map<PolyT, CommOuts>(h, mkComOut);


  SumcheckPf::EqProofs eqPfs(crs->d);

  CommOuts zComOut(crs->d+1);
  zComOut[0] = yComOut;

  for (auto i = 0; i < crs->d; i++) {
    auto hComEvalOn0 = PolyT::evalAsPolyOn(hComOut[i], In::zero());
    auto hComEvalOn1 = PolyT::evalAsPolyOn(hComOut[i], In::one());
    auto vComOut = hComEvalOn0 + hComEvalOn1;

    eqPfs[i] = make_shared<ZKEqProof>(comScm, vComOut, zComOut[i]);

    zComOut[i+1] = PolyT::evalAsPolyOn(hComOut[i], r[i]);
  }

  /* Poly Proofs */
  CommOuts cmout_g(n_cm_polys), cmout_eval(n_cm_polys);
  vector<PolyPf> polypf(n_cm_polys);
  for (auto i = 0; i < n_cm_polys; i++) {
    auto &mle_poly = mles[i]->getV();
    cppoly->computeAnswer(cmout_eval[i], r, mle_poly);
    cmout_g[i] = cppoly->commitPoly(mle_poly);
    cppoly->prove(mle_poly, cmout_eval[i], r, polypf[i]);
  }

  /* Product Proof */
  auto betaEval = beta->evalOnPoint(rho, r);
  auto lhsProd = cmout_eval[0]*betaEval;
  ZKPrdProof prdPf(comScm, lhsProd, cmout_eval[1], zComOut[crs->d]);

  stopBenchmark("prove");

  /* Final wrap up */

  // just extract commitments frp, hComOut
  auto mkCom = [](const CommOuts &cOuts) { return CommOut::toComms(cOuts); };
  auto hCom = cputil::map<CommOuts, Comms>(hComOut, mkCom);

  return new SumcheckPf(
          r, hCom, eqPfs, CommOut::toComms(cmout_g),
          CommOut::toComms({lhsProd, cmout_eval[1]}), polypf, prdPf);

}


bool scCheckCommit(CPPoly *cppoly, const Comm &cm)
{
  return cppoly->checkCommit(cm);
}


bool checkEqPf(Comm c0, Comm c1, shared_ptr<ZKEqProof> eqPf)
{
  return (c0.c == eqPf->c0) && (c1.c == eqPf->c1) && eqPf->verify(); // possibly move c0.c and c1.c as parameters?
}

bool CPSumcheck::verify(const SumcheckKey *crs, const CPVIn &in, const SumcheckPf *pf)
{
  // verifier gear
  vector <bool> checks;
  auto addCheck = [&checks](bool b) { checks.push_back(b); };
  auto idFn = [](bool b) { return b; };
  auto checkAll = [checks, idFn]() { return all_of(begin(checks), end(checks), idFn); };

  // setup inputs
  // inputs are: commitment to y and commitment to g
  auto &c = in.commIn;
  auto yCom = c[0];

  auto &hCom = pf->hCom;

  startBenchmark("verify");
  // Verification loop
  Comms zCom(crs->d+1);
  zCom[0] = yCom;

  for (auto i = 0; i < crs->d; i++)
  {
    auto hComEvalOn0 = PolyT::evalAsPolyOn(hCom[i], In::zero());
    auto hComEvalOn1 = PolyT::evalAsPolyOn(hCom[i], In::one());
    auto vCom = hComEvalOn0 + hComEvalOn1;
    addCheck(checkEqPf(vCom, zCom[i], pf->eqPfs[i]));

    zCom[i+1] = PolyT::evalAsPolyOn(hCom[i], pf->r[i]);
  }


  // check commitments and proofs
  for (auto i = 0; i < n_cm_polys; i++) {
    //addCheck(checkCommit(cppoly, c[i]));
    addCheck(scCheckCommit(cppoly, pf->polycm_evalg[i]));
    addCheck(cppoly->verify(pf->polycm_g[i], pf->polycm_evalg[i], pf->r, pf->polypf_g[i]));
  }

  addCheck(pf->finalPrdPf.verify());

  stopBenchmark("verify");


  return checkAll();
}
