#include "hadamardsc.h"
#include "util.h"


CPPIn makeSCPIn(const vector<CommRand> &rho, const CommOut &aCmOut,
				const CommOut &bCmOut, const CommOut &yCmOut)
{
	CPPIn scInput;
	scInput.publicSlot = rho;
	scInput.commSlot.push_back(yCmOut); // NB: This is not initialized yet
  scInput.commSlot.push_back(aCmOut);
  scInput.commSlot.push_back(bCmOut);

  return scInput;
}

CPVIn makeSCVIn(const vector<CommRand> &rho,
				const Comm &aCm, const Comm &bCm, const Comm &yCm)
{
	CPVIn scIn;
	scIn.publicIn = rho;
	scIn.commIn = {yCm, aCm, bCm };
	return scIn;
}


void cppolyProve(CPPoly *cppoly, const Ins &v, const HadRand &point, CommOut &cOut, PolyPf &polyPf)
{
	cppoly->computeAnswer(cOut, point, v);
	cppoly->prove(v, cOut, point, polyPf);

}


// verifyPolyPf(cppoly, pf->r, c[0], pf->polyAnsComms[0], pf->polyProofs[0]);
bool verifyPolyPf(
	CPPoly *cppoly, const HadRand &r,
	const Comm &cmPoly, const Comm &cmVal,
	const PolyPf &polyPf)
{
	return cppoly->verify(cmPoly, cmVal, r, polyPf);
}

// rel here is just long*
HadKey* CPHad::keygen(const HadRel *rel)
{
	size_t logRel = cputil::log2ceiled(rel->n);

	auto scCrs = cpsumcheck.keygen(&logRel);

	return new HadKey(rel->n, logRel, cppoly, scCrs);
}

HadPf* CPHad::prove(const HadKey *crs, const CPPIn &in)
{
	HadPf *hadPf = new HadPf(crs->d);

	/*
	 * The algorithm:
	 * - Prove with cp-poly the result of \tilde{u}_0(r) (call this result t)
	 * - Prove with cp-sc that t is also the result of sumcheck on the three poly-s:
	 *   1) eq_r(X) // equality predicate with r fixed
	 *   2) \tilde{u}_1
	 *   3) \tilde{u}_2
	 *  -
	 */

	auto c = CommOut::toComms(in.commSlot);
	auto u = CommOut::toCommittedVals(in.commSlot);

  auto &aCmOut = in.commSlot[1];
  auto &bCmOut = in.commSlot[2];

  CommOut uProdEvalCmOut;

	// sample r (this will be rho for cpsumcheck)
	for(auto i=0;i<hadPf->r.size();i++) {
		hadPf->r[i] = CommRand::random_element();
	}
  CPPIn scInput = makeSCPIn(hadPf->r, aCmOut, bCmOut, uProdEvalCmOut);

	// NB: the first (0-th input is the result of product)

	// poly proof for result
	startBenchmark("prove_cppoly");
	cppolyProve(cppoly, u[0], hadPf->r, uProdEvalCmOut, hadPf->polyProof);
	stopBenchmark("prove_cppoly");
	// sumcheck proof
	scInput.commSlot[0] = uProdEvalCmOut; // Because it's just been initialized
  hadPf->sumcheckPf = cpsumcheck.prove(crs->scCrs, scInput);
	applyBenchmarkFrom(cpsumcheck, "prove", "prove_sc");

	// "Flush" commitments into proof we return
  hadPf->uProdEvalCm = uProdEvalCmOut.c;
  //stopBenchmark("prove");

	return hadPf;
}


bool CPHad::verify(const HadKey *crs, const CPVIn &in, const HadPf *pf)
{
	bool isGdPf;
	vector <bool> checks;
	auto addCheck = [&checks](bool b) { checks.push_back(b); };
	auto idFn = [](bool b) { return b; };
	auto checkAll = [checks, idFn]() { return all_of(begin(checks), end(checks), idFn); };

	// NB: the first (0-th input is the result of product)

	CPPoly *cppoly = crs->cppoly;

	auto &c = in.commIn;

	CPVIn scIn = makeSCVIn( pf->r, c[1], c[2], pf->uProdEvalCm );


	// Verify Sumcheck Proof
	bool isSCPfGd = cpsumcheck.verify(crs->scCrs, scIn, pf->sumcheckPf);
	addCheck(isSCPfGd);
	applyBenchmarkFrom(cpsumcheck, "verify", "verify_sc");

	// Verify  Poly Proof
	startBenchmark("verify_cppoly");
	bool isPf0Gd = verifyPolyPf(cppoly, pf->r, c[0], pf->uProdEvalCm, pf->polyProof);
	addCheck(isPf0Gd);

	stopBenchmark("verify_cppoly");

	return checkAll();
}


