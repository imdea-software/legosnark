#include "commit.h"
#include "dbgutil.h"

#include "hadamardsc.h"
#include "lipmaa.h"
#include "benchmark.h"

#include <chrono>
#include <thread>
#include <cstdlib>
using namespace std;

#include "fmt/format.h"


void mle_test()
{
  // v = [1, 0 ,0 ,2, 0, 0, 1, 0]
  In f0 = In::zero();
  In f1 = In::one();

  Ins v(8, f0);
  v[0] = f1; // eq_{0,0,0} = [-1, 1, 1, -1, 1, -1, -1, 1]
  v[4] = f1*2; // eq_{1,0,0} = [1, -1, -1, 1, 0,...,0]
  v[6] = f1; // eq_{1,1,0} = [-1,1, 0, ..., 0]

}

void sc_test(CommScheme *cmScm)
{

  // v = [1, 0 ,0 ,2, 0, 0, 1, 0]
  In f0 = In::zero();
  In f1 = In::one();

  const size_t v_sz = 8;
  const size_t d = 3;
  Ins v(v_sz, f0);
  v[0] = f1; // eq_{0,0,0} = [-1, 1, 1, -1, 1, -1, -1, 1]
  v[4] = f1*2; // eq_{1,0,0} = [1, -1, -1, 1, 0,...,0]
  v[6] = f1; // eq_{1,1,0} = [-1,1, 0, ..., 0]


}

void hadlipmaa(const Ins &a, const Ins &b, const Ins &c)
{
  const uint64 n = a.size();

  InterpCommScheme ics;
  CPHadL cphadl;
  auto pBm = make_shared<Benchmark>();
  cphadl.setBenchmark(pBm, "CPHadLipmaa");

  auto interp = new Interpolator(n);
  // Trapdoors
  IScalar chi = IScalar::random_element();
  IScalar gamma = IScalar::random_element();

  ics.keygen(n, *interp, chi, gamma);

  cphadl.keygen(n, *interp, chi, gamma);

  auto cmOuta = ics.commit(a);
  auto cmOutb = ics.commit(b);
  auto cmOutc = ics.commit(c);

  auto pf = cphadl.prove(cmOuta, cmOutb, cmOutc);

  cout << "## ---" << endl;
  print_bm("##had_lipmaa Prove", "prove", cphadl);

  bool isGd = cphadl.verify(pf, cmOuta.c, cmOutb.c, cmOutc.c);
  print_bm("##had_lipmaa Verify", "verify", cphadl);

}

void hadsc(const Ins &a, const Ins &b, const Ins &c)
{
  const long n = a.size();
  CommScheme *commScm = new CommScheme;
  commScm->keygen(n);

  CPPIn proverInput;
  CPVIn verifInput;
  CPInputFmt::init_no_pub(proverInput, verifInput, commScm, {c, a, b});


  // actual keygen, proving, verifying
  CPHad had(commScm, new CPPoly(commScm)); // XXX: Check how we are passing CPPoly
  auto pBm = make_shared<Benchmark>();
  had.setBenchmark(pBm, "CPHadSumcheck");
  auto crs = had.keygen(new HadRel(n));

  auto pf = had.prove(crs, proverInput);
  bool isGdPf = had.verify(crs, verifInput, pf);

  print_bm("##had_sc (Sumcheck) Prove", "prove_sc", had);
  print_bm("##had_sc (CPPoly) Prove", "prove_cppoly", had);
  print_sum_bm("##had_sc TOTAL Prove", "prove_sc", "prove_cppoly", had);
  cout << "##" << endl;
  print_bm("##had_sc (Sumcheck) Verify", "verify_sc", had);
  print_bm("##had_sc (CPPoly) Verify", "verify_cppoly", had);
  print_sum_bm("##had_sc TOTAL Verify", "verify_sc", "verify_cppoly", had);
  cout << "## ## ##" << endl;
}

int main(int argc, char **argv){

  default_ec_pp::init_public_params();

  size_t MIN_D, MAX_D;
  MIN_D = MAX_D = 3;


  if (argc == 2) {
    MIN_D = MAX_D = stoi(argv[1]);
  } else if (argc == 3) {
    MIN_D = stoi(argv[1]);
    MAX_D = stoi(argv[2]);
  }

  for (size_t d = MIN_D; d <= MAX_D; d++) {
    const uint64 n = 1 << d; // 2**d

    cout << "## Vector size: " << n << endl;


    // prepare input
    Ins u(n);
    Ins uSqrd(n);
    for (auto i = 0; i < n; i++) {
      u[i] = LFr::one() * i;
      uSqrd[i] = LFr::one() * i * i;
    }

    hadsc(u,u, uSqrd);
    hadlipmaa(u,u,uSqrd);
  }

  return 0;
}
