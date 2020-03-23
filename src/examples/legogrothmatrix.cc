#include <cstdio>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <utility>
#include <fstream>
#include <type_traits>
using namespace std;

#include "globl.h"
#include "benchmark.h"

#include <libff/common/default_types/ec_pp.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>

#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

using namespace libsnark;
using namespace libff;

const int nreps = 1;

template<typename ppT, typename FieldT>
void run_it(protoboard<FieldT> &pb, const int n_sqrd);

unsigned rand32b()
{
	return rand() % 0xFFFFFFFF;
}

unsigned rand16b()
{
	return rand() % 0xFFFF;
}


unsigned rand2b()
{
	return rand() % 2;
}


vector<LG1> random_base(size_t m)
{
	vector<LG1> ret(m);
	for (auto i = 0; i < m; i++) {
		ret[i] = LFr::random_element()*LG1::one();
	}
	return ret;
}

vector<LFr> random_scalars(size_t m)
{
	vector<LFr> ret(m);
	for (auto i=0; i < m; i++) {
		ret[i] = LFr::one()*rand32b();
	}
	return ret;
}

template<typename ppT>
bool run_sha_snark(const int n)
{
  using FieldT = Fr<ppT>;
  using VArray = pb_variable_array<FieldT>;
  
  auto n_sqrd = n*n;

  protoboard<FieldT> pb;
  
  vector<VArray> M_rows(n), N_cols(n);
  vector<vector<pb_variable<FieldT>>> U; // result matrix
  vector<vector<inner_product_gadget<FieldT>>> inner_products;
  
  
  for (auto i = 0; i < n; i++) {
    M_rows[i].allocate(pb, n, "M_row");
    N_cols[i].allocate(pb, n, "N_col");
    U.push_back(vector<pb_variable<FieldT>>(n));
    inner_products.push_back(vector<inner_product_gadget<FieldT>>());
  }
   
  for (auto r = 0; r < n; r++) {
    for (auto c = 0; c < n; c++) {
      U[r][c].allocate(pb, "U_elt");
      inner_products[r].push_back(
        inner_product_gadget<FieldT>(
          pb, M_rows[r], N_cols[c], U[r][c], "inner_product"));
    }
  }

  // set U as the primary input
  pb.set_input_sizes(0);

  // gen constraints 
  for (auto r = 0; r < n; r++) {
    for (auto c = 0; c < n; c++) {
      inner_products[r][c].generate_r1cs_constraints();
    }
  }
  
  // gen witness
  for (auto i = 0; i < n; i++) {
    for (auto j = 0; j < n; j++) {
      pb.val(M_rows[i][j]) = FieldT::one()*rand32b(); //FieldT::random_element();
      pb.val(N_cols[i][j]) = FieldT::one()*rand32b(); //FieldT::random_element();
    }
  }
  
  for (auto r = 0; r < n; r++) {
    for (auto c = 0; c < n; c++) {
      inner_products[r][c].generate_r1cs_witness();
    }
  }


  assert(pb.is_satisfied());
  
  
  run_it<ppT,FieldT>(pb, n_sqrd);
}



template<typename ppT, typename FieldT>
void run_it(protoboard<FieldT> &pb, const int n_sqrd)
{
  vector<LG1> u1 = random_base(3*n_sqrd+1);
  vector<LG1> u2 = random_base(3*n_sqrd+2);
  vector<LFr> e1 = random_scalars(3*n_sqrd+1);
  vector<LFr> e2 = random_scalars(3*n_sqrd+2);

  r1cs_gg_ppzksnark_keypair<ppT> keypair = r1cs_gg_ppzksnark_generator<ppT>(pb.get_constraint_system());
  r1cs_gg_ppzksnark_processed_verification_key<ppT> pvk;
  r1cs_gg_ppzksnark_proof<ppT> proof;
  bool ans;


  auto kgFn = [&]() {
    r1cs_gg_ppzksnark_keypair<ppT> keypair1 = r1cs_gg_ppzksnark_generator<ppT>(pb.get_constraint_system());
    pvk = r1cs_gg_ppzksnark_verifier_process_vk<ppT>(keypair.vk);
	};
  
  LG1 res1, res2;
	auto prvFn = [&]() {
    proof =	r1cs_gg_ppzksnark_prover<ppT>(
      keypair.pk, pb.primary_input(), pb.auxiliary_input());
    res1 = multiExpMA<LG1>(u1, e1);
    res2 = multiExpMA<LG1>(u2, e2);
	};
  
  auto nComms = 2;
  vector<LG1> v{LFr::random_element()*LG1::one(), LFr::random_element()*LG1::one()};
  vector<LG2> v2{LFr::random_element()*LG2::one(), LFr::random_element()*LG2::one()};
  vector<G2_precomp<def_ec>> v2p(2);
  for (auto i = 0; i < nComms; i++) {
      v2p[i] = def_ec::precompute_G2(v2[i]);
  }
  auto pf = LFr::random_element()*LG1::one();
	auto verFn = [&]() {
    ans = r1cs_gg_ppzksnark_online_verifier_strong_IC<ppT>(pvk, pb.primary_input(), proof);
    vector<G1_precomp<def_ec>> com_precomp(nComms);
    for (auto i = 0; i < nComms; i++) {
      com_precomp[i] = def_ec::precompute_G1(v[i]);
    }
    Fqk<def_ec> lhs;
      
    
    lhs= def_ec::double_miller_loop(com_precomp[0], v2p[0], com_precomp[1], v2p[1]);
    lhs= lhs * def_ec::miller_loop(com_precomp[0], v2p[0]);
    auto pf_precomp = def_ec::precompute_G1(pf);
    auto rhs = def_ec::miller_loop(pf_precomp, v2p[0]);
    auto shouldBeOne = def_ec::final_exponentiation(lhs*rhs.unitary_inverse());
	};
 
  libff::print_header("R1CS GG-ppzkSNARK Generator");
  printf("KG Time: (%f s)\n", TimeDelta::runAndAverage(kgFn, nreps)/1000000);

  libff::print_header("R1CS GG-ppzkSNARK Prover");
  printf("Prv Time: (%f s)\n", TimeDelta::runAndAverage(prvFn, nreps)/1000000);


  libff::print_header("R1CS GG-ppzkSNARK Online Verifier");
  printf("Ver Time: (%f s)\n", TimeDelta::runAndAverage(verFn, nreps)/1000000);


}

int main () {
  default_r1cs_gg_ppzksnark_pp::init_public_params();
  
  auto lmin = 2;
  auto lmax = 7;
  for (auto l = lmin; l <= lmax; l++) {
    auto n = 1 << l;
    cout << "Measuring Time for n = " << n << endl;
    run_sha_snark<def_ec>(n);
  }

  return 0;
}
