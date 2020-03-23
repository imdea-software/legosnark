#define OUTPUT_MATRIX_IN_CLEAR


#include "commit.h"
#include "dbgutil.h"

#include "matrixsc.h"
#include "lipmaa.h"
#include "benchmark.h"

#include <chrono>
#include <thread>
#include <cstdlib>
using namespace std;

#include "fmt/format.h"

void matsc(const Ins &a, const Ins &b, const Ins &c)
{
  const long n = a.size();
  CommScheme *commScm = new CommScheme;
  commScm->keygen(n);

  CPPIn proverInput;
  CPVIn verifInput;
  CPInputFmt::init_no_pub(proverInput, verifInput, commScm, {c, a, b});


  // actual keygen, proving, verifying
  CPMat mat(commScm, new CPPoly(commScm)); // XXX: Check how we are passing CPPoly
  auto pBm = make_shared<Benchmark>();
  mat.setBenchmark(pBm, "CPMatSumcheck");
  auto crs = mat.keygen(&n);

  auto pf = mat.proveOutputMatrixInClear(crs, proverInput);
  bool isGdPf = mat.verifyOutputMatrixInClear(crs, verifInput, c, pf);

  print_bm("##mat_sc (Sumcheck) Prove:", "prove_sc", mat);
  print_bm("##mat_sc (eval) Prove:", "prove_cppoly", mat);
  print_sum_bm("##mat_sc TOTAL Prove:", "prove_sc", "prove_cppoly", mat);
  cout << "##" << endl;
  print_bm("##mat_sc (Sumcheck) Verify:", "verify_sc", mat);
  print_bm("##mat_sc (eval) Verify:", "verify_eval", mat);
  print_sum_bm("##mat_sc TOTAL Verify:", "verify_sc", "verify_eval", mat);
  cout << "##" << endl;
  cout << "##mat_sc Proof Size: " << pf->getSize() << endl;
  cout << "## ## ##" << endl;
}

unsigned rand32b()
{
  return rand() % 0xFFFFFFFF;
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
    auto N = n*n;

    cout << "## Matrix size: " << N << " (n=" << n << ")\n";


    // prepare input
    Ins A(N), B(N), C(N);
    for (auto i = 0; i < N; i++) {
      A[i] = LFr::one() * rand32b();
      B[i] = LFr::one() * rand32b();
    }
    for (auto i = 0; i < n; i++) {
      for (auto j = 0; j < n; j++) {
        auto idx_C = i*n + j;
        C[idx_C] = LFr::zero();
        for (auto k = 0; k < n; k++) {
          C[idx_C] = C[idx_C] + A[i*n+k]*B[k*n+j];
        }
      }
    }

    matsc(A,B, C);
  }

  return 0;
}
