#include "arithcirc.h"
#include "matrix.h"
#include "util.h"

#include <utility>
#include <algorithm>
using namespace std;

#include <libff/common/default_types/ec_pp.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
using namespace libff;



template<typename T>
void expandTripletsFromPairs(
	size_t r_offset, size_t c_offset,
	const RCPairs &rcs, const vector<T> &vals,
	vector<vector<CoeffPos<T>>> &M)
{
	MYREQUIRE(rcs.size() == vals.size());
	for (auto i = 0; i < rcs.size(); i++) {
		T v = vals[i];
		size_t r = rcs[i].first;
		size_t c = rcs[i].second;
		insertAsColMajor(r_offset+r, c_offset+c, v, M);
	}
}


// M is the output and its size should be initialized as (F.rows()+nComms, F.cols()+nComms)
// NB: this version assumes that the blinding factor h0 and the base is the same for all commitments

template<typename T>
void mkSubspaceMatrixVEq(
	size_t nComms,
	const RCPairs &posW, const vector<T> &valsW,
	const T &h0, const vector<vector<T>> &bases,
	vector<vector<CoeffPos<T>>> &M)
{
  
  /* Goal matrix M should look like the following: (this pic is for nComms==2)
  -----------------------------
  | h0 0 |  Base        0      |
  | 0 h0 |   0         Base    |
  | - - - - - - - - - - - - - -|
  |      |                     |
  |   0  |          F          |
  |      |                     |
   -----------------------------
  */
  
  
  auto tot_base_sz = 0;
  for (auto i = 0; i < nComms; i++) {
    tot_base_sz += bases[i].size();
  }
  // set right # of cols
  M.resize(nComms + tot_base_sz);
  

  // top-left
  for (auto i = 0; i < nComms; i++) {
	insertAsColMajor(i, i, h0, M);
  }

  // top-right
  size_t colOffset = nComms;
  for (auto i = 0; i < nComms; i++) {
    insertRowAsColMajor(i, colOffset, bases[i], M);
    colOffset += bases[i].size();
  }

  // bottom-right
  expandTripletsFromPairs(nComms, nComms, posW, valsW, M); 


}

void mkAuxPrecompTgtV(Fqk<def_ec> &aux_precomp, const vector<LG1> &V, const vector<G2_precomp<def_ec>> &C)
{
	MYREQUIRE(V.size() == C.size()-3);
	size_t t = V.size();
	size_t thalf = t/2;
	aux_precomp = Fqk<def_ec>::one();
	// We go two-by-two as far as we can
	for (auto i =0; i < 2*thalf; i +=2) {
		auto iC = i+3; // shift index for C
		auto Vi_precomp = def_ec::precompute_G1(V[i]);
		auto Vip1_precomp = def_ec::precompute_G1(V[i+1]);
		auto tmp = 
			def_ec::double_miller_loop(
				Vi_precomp, C[iC], Vip1_precomp, C[iC+1]);
		aux_precomp = aux_precomp * tmp;
	}
	// if t is odd then we need one more pairing
	if (t % 2 != 0) {
		auto V_last_precomp = def_ec::precompute_G1(V[t-1]);
		auto tmp = def_ec::miller_loop(V_last_precomp, C[t-1+3]);
		aux_precomp = aux_precomp * tmp;
	}
	
}


ACKey* CPAC::keygen(const ACRel *rel) 
{
	fmt::print("n={}; r={}; c={}\n", rel->n_mul_gates, rel->r, rel->c);
	auto n = rel->n_mul_gates;
	ACKey *crs = new ACKey;
	
	startBenchmark("keygen_comm");
	auto ics = dynamic_cast<InterpCommScheme *>(getCommScheme());
	if (!ics) {

	}
	
	auto nEntries = rel->Wtrs.size();
	fmt::print("nEntries={}\n", nEntries);
	
	// begin joint kg for commitment/hadamard
	long unsigned N = max( (long unsigned)n, nEntries);
	N = max(N, rel->tgtV.size()); 
	interp = new Interpolator(n, N); 
	// Trapdoors
	IScalar chi = IScalar::random_element();
	IScalar gamma = IScalar::random_element();

	ics->keygen(n, *interp, chi, gamma);
	stopBenchmark("keygen_comm");
	
	cphadl.keygen(n, *interp, chi, gamma);
	crs->hadkey = &cphadl.key;
	// end joint kg for commitment/hadamard
	
	applyBenchmarkFrom(cphadl, "keygen", "keygen_had");
	
	cout << "Finished KG Had\n";
	
	// let us compute W in G1

	RCPairs rcs(nEntries);
	Ins scalarsW(nEntries);
	for (auto i = 0; i < nEntries; i++) {
		rcs[i] = make_pair(rel->Wtrs[i].row(), rel->Wtrs[i].col());
		scalarsW[i] = rel->Wtrs[i].value();
	}
	
	
	startBenchmark("keygen_prep_veq");
  if (rel->has_tgtV) {
    crs->g1TgtV = interp->mkG1Exp(rel->tgtV);
  }
  vector<LG1>  groupvalsW = interp->mkG1Exp(scalarsW); // XXX: Check if you can use smaller window.
	stopBenchmark("keygen_prep_veq");


	cout << "Done multiexp for W's values\n";
  auto nComms = (rel->input_size > 0) ? 4 : 3;
	auto ssR = rel->r+nComms;
	auto ssC = rel->c+nComms;
	
  size_t C_precomp_sz;
  if (rel->has_tgtV) {
    C_precomp_sz = ssR;
  } else {
    C_precomp_sz = nComms;
  }
  crs->has_tgtV = rel->has_tgtV;
  crs->input_size = rel->input_size;

  
  // simple book-keeping of data structures here
	SubspaceRel *ssrel = 
    new SubspaceRel(ssR, ssC, vector<ColFr>(ssC), vector<ColG1>(ssC), C_precomp_sz, interp); 
  vector<vector<LG1>> g1_bases;
  vector<vector<LFr>> sc_bases;
  for (auto i = 0; i < nComms; i++) {
    if (i == 3) { // nComms == 4 
      // put special base
      vector<LG1> g1_ubase(ics->key.lg1.begin(), ics->key.lg1.begin()+rel->input_size);
      g1_bases.push_back(g1_ubase);
      vector<LFr> sc_ubase(ics->key.l.begin(), ics->key.l.begin()+rel->input_size);
      sc_bases.push_back(sc_ubase);
    } else {
      g1_bases.push_back(ics->key.lg1);
      sc_bases.push_back(ics->key.l);
    }
  }
  
	mkSubspaceMatrixVEq(nComms, rcs, groupvalsW, ics->key.zg1, g1_bases, ssrel->M);
	// do same for scalar
	mkSubspaceMatrixVEq(nComms, rcs, scalarsW, ics->key.z, sc_bases, ssrel->sM);

	
	startBenchmark("keygen_ss");
	crs->sskey = ss.keygen(ssrel);
	stopBenchmark("keygen_ss");

	startBenchmark("keygen_aux");
  if (rel->has_tgtV) {
    mkAuxPrecompTgtV(crs->tgtV_C_precomp, crs->g1TgtV, crs->sskey->C_precomp);
  }
	stopBenchmark("keygen_aux");

	cout << "Done KG VEq" << endl;
  
  ics->print_key_size();
  crs->sskey->print_size();
  cphadl.print_key_size();

	return crs;
}
ACPf* CPAC::prove(const ACKey *crs, const ACPIn &in) 
{
	ACPf *pf = new ACPf;

	// Hadamard
	pf->hadpf = cphadl.prove(in.aCOut, in.bCOut, in.cCOut);
	applyBenchmarkFrom(cphadl, "prove", "prove_had");

	// VEq
	// We first make a subspace witness prepending openings
	vector<LFr> w{in.aCOut.r, in.bCOut.r, in.cCOut.r};
  if (crs->input_size != 0) { // if we have u
    w.push_back(in.uCOut.r);
  }
	vector<LFr> abc =
		cputil::concat3(in.aCOut.xs, in.bCOut.xs, in.cCOut.xs);
  if (crs->input_size != 0) { // if we have u
    abc.insert(abc.end(), in.uCOut.xs.begin(), in.uCOut.xs.end());
  }
	w.insert(w.end(), abc.begin(), abc.end()); 
	auto aSize = in.aCOut.xs.size();

	
	startBenchmark("prove_veq");
	pf->sspf = ss.prove(crs->sskey, w);
	stopBenchmark("prove_veq");
	
	return pf;
}

bool CPAC::verify(const ACKey *crs, const ACVIn &in, const ACPf *pf)
{
  auto ics = getCommScheme();
  startBenchmark("verify_commit");
  bool b;
  b = b && ics->verify(in.aC);
  b = b && ics->verify(in.bC);
  b = b && ics->verify(in.cC);
  if (crs->input_size != 0) { // if we have u
    b = b && ics->verify(in.uC);  }
  stopBenchmark("verify_commit");
  
	bool hadVfy = cphadl.verify(pf->hadpf, in.aC, in.bC, in.cC);
	applyBenchmarkFrom(cphadl, "verify", "verify_had");
	
	vector<LG1> x {in.aC.c, in.bC.c, in.cC.c};
  if (crs->input_size != 0) { // if we have u
    x.push_back(in.uC.c);
  }
	startBenchmark("verify_veq");
  bool ssVfy;
  if(crs->has_tgtV) 
    ssVfy = ss.verifyLin3or4(crs->sskey, x, &crs->tgtV_C_precomp, pf->sspf);
  else 
    ssVfy = ss.verifyLin3or4(crs->sskey, x, nullptr, pf->sspf);
	stopBenchmark("verify_veq");

	return hadVfy && ssVfy;
	
}
