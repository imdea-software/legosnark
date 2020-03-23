#ifndef GLOBL_H

#define GLOBL_H

#include <libff/common/default_types/ec_pp.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libfqfft/polynomial_arithmetic/basic_operations.hpp>


using namespace libfqfft;
using namespace libff;

#include <vector>
#include <cstdlib>
#include <algorithm>
using std::runtime_error;
using std::vector;
using std::begin;
using std::end;
using std::min;

using uint64 = unsigned long long;

//default_ec_pp::init_public_params();
using def_ec = default_ec_pp;

// XXX: Might be changed later
using LG1 = G1<def_ec>;
using LG2 = G2<def_ec>;
using LFr = Fr<def_ec>;
using LGT = GT<def_ec>;



inline void MYREQUIRE(bool b) {
  if (!b) {
    throw runtime_error("Error!");
  }
}

template<typename T>
inline T mulEcByScalar(const T &pt, const LFr &x) {
	return x*pt;
}


template<typename G>
G multiExp(const vector<G> &gs, const  vector<LFr> &xs)
{
  size_t n = min(gs.size(), xs.size());
	#ifdef MULTICORE
    const size_t chunks = omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call omp_set_num_threads()
	#else
    const size_t chunks = 1;
	#endif

	return libff::multi_exp<G, LFr, libff::multi_exp_method_bos_coster>(
	  gs.begin(), gs.begin()+n,
	  xs.begin(), xs.begin()+n,
	  chunks);
}

template<typename G>
G multiExpMA(const vector<G> &gs, const vector<LFr> &xs)
{
  size_t n = min(gs.size(), xs.size());
	#ifdef MULTICORE
    const size_t chunks = omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call omp_set_num_threads()
	#else
    const size_t chunks = 1;
	#endif
      printf("NCHUNKS : %d\n", chunks);

	return libff::multi_exp_with_mixed_addition<G, LFr, libff::multi_exp_method_BDLO12>(
	  gs.begin(), gs.begin()+n,
	  xs.begin(), xs.begin()+n,
	  chunks);
}


 inline vector<LFr> hadamard(const vector<LFr> &a, const vector<LFr> &b)
{
  MYREQUIRE(a.size() == b.size());
  auto n = a.size();
  vector<LFr> ret(n);
  for (auto i = 0; i < n; i++) {
	ret[i] = a[i]*b[i];
  }
  return ret;
}


// returns if e(a1, a2) == e(b1, b2)
inline bool simple_pairing_check(LG1 a1, LG2 a2, LG1 b1, LG2 b2)
{
   auto a1_pre = def_ec::precompute_G1(a1);
   auto b1_pre = def_ec::precompute_G1(b1);
   auto a2_pre = def_ec::precompute_G2(a2);
   auto b2_pre = def_ec::precompute_G2(b2);

   auto lhs = def_ec::miller_loop(a1_pre, a2_pre);
   auto rhs = def_ec::miller_loop(b1_pre, b2_pre);
   auto out = def_ec::final_exponentiation(lhs*rhs.unitary_inverse());
   return out == LGT::one();
}



#endif
