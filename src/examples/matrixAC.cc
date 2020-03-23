#include <cstdio>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <utility>
#include <type_traits>
using namespace std;

#include "globl.h"
#include "lipmaa.h"
#include "bp_circuits.h"
#include "matrix.h"
#include "dbgutil.h"
#include "subspace.h"
#include "arithcirc.h"
#include "util.h"



const int NREPS = 1;


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

void addI(size_t offset_r, size_t offset_c, vector<InT> &trs, size_t m)
{
  for (auto i = 0; i < m; i++) {
    addToTriplets(offset_r+i, offset_c+i, LFr::one(), trs);
  }
}

void addMI(size_t offset_r, size_t offset_c, vector<InT> &trs, size_t m)
{
  for (auto i = 0; i < m; i++) {
    addToTriplets(offset_r+i, offset_c+i, LFr::one()*(-1), trs);
  }
}

void mkMatrixMultF(ScalarMat &F,  vector<InT> &trs, size_t n)
{
  // Make WL by introducing a large I in the middle
  addI(n*n, 0, trs, n*n*n);
  // Make WR by introducing a large I in the end
  addI(n*n+n*n*n, n*n*n, trs, n*n*n);
  
  vector<LFr> n_ones(n, LFr::one());
  vector<LFr> n_neg_ones(n, LFr::one()*(-1));

  auto c_WO_offset = 2*n*n*n;
  
  // mk WO with n ones at each row, positioned "diagonally:
  MYREQUIRE(F.rows() == n*n+2*n*n*n);
  for (auto r = 0; r < n*n; r++) {
    addRowToTriplets(r, r*n + c_WO_offset, n_ones, trs);
  }
  
  // mk WU which is basically -I, where I is the identity matrix
  auto c_WU_offset = 3*n*n*n;

  for (auto r = 0; r < n*n; r++) {
    addToTriplets(r, r + c_WU_offset, LFr::one()*(-1), trs);
  }
  
  // M
  auto r_u_offset = n*n;
  auto c_u_offset = n*n + c_WU_offset;
  for (auto c = 0; c < n*n; c++) {
    addColToTriplets(r_u_offset+c*n, c_u_offset + c, n_neg_ones, trs);
  }
  
  auto r_M_offset = n*n*n + r_u_offset;
  auto c_M_offset = n*n + c_u_offset;
  
  // N
  // do n times a descending scale of -I of size n
  for (auto i = 0; i < n; i++)  {
    for (auto j = 0; j < n; j++) {
        addMI(r_M_offset+i*n*n + j*n, c_M_offset+j*n, trs, n);
    }
  }
  
      
  // finally mk the matrix
  initFromTriplets(F, trs);
}

Ins get_random_vec(long n)
{
  Ins ret(n);
  for (auto i = 0; i < n; i++) {
    ret[i] = LFr::one()*rand32b();
  }
  return ret;
}


Ins findTarget(const ScalarMat &F, const Ins &a, const Ins &b, const Ins &c)
{
  // concatenate a, b and c
  Ins abc = cputil::concat3(a, b, c);
  // convert the result to a col vector
  auto abcVec = mkColVector(abc.size(), abc);

  // do the multiplication
  ScalarVec tgt = F*abcVec;

  //bool isGd = satisfyLinTransform(abcVec, F, tgt);
  //MYREQUIRE(isGd);

  // reconvert to std vector
  Ins ret;
  colToStdVec(ret, tgt);
  return ret;
}


void print_time(string msgtag, double t)
{

	fmt::print("{}: {} micros ({} s)\n", msgtag, t, t/1000000);
}

void print_time_string_onerep(string msgtag, string bmlbl, const Benchmarkable &obj)
{
	double t;
  if (bmlbl == "commit") {	// NB: Special handling

    t =  obj.getTimingInMicrosFor(bmlbl + "_a");
		t *= 3;
    print_time(msgtag + " abc", t);
    t =  obj.getTimingInMicrosFor(bmlbl + "_u");
    print_time(msgtag + " u", t);
    return;

  } else {
    t =  obj.getTimingInMicrosFor(bmlbl);
  }
	
	print_time(msgtag, t);
}

// our matrix is n by n
void benchmark_bp_circuit(size_t n)
{
  // generate F 
  cout << "Building up Ws..." << endl;
  vector<InT> trsF;
  auto nMulGates = n*n*n;
  auto nConstraints = n*n+2*n*n*n;
  auto nInSize = 3*n*n ; // the product matrix
  auto nrowsF = nConstraints;
  auto ncolsF = 3*nMulGates + nInSize;
  ScalarMat F(nrowsF, ncolsF);
  mkMatrixMultF(F, trsF, n);
  
  //cout << F << endl;

  // generate fake inputs
  cout << "Generating Inputs..." << endl;
  vector<vector<LFr>> M, N;
  M.resize(n);
  N.resize(n);
  for (auto i = 0; i < n; i++) {
    M[i] = get_random_vec(n); // rows of M
    N[i] = get_random_vec(n); // cols of N
  }
  
  // make a = (M_1, M_1,..., M_1, M_2,..., M_2, ..., M_n, ..., M_n)
  vector<LFr> a;
  // for each row in M
  for (auto r = 0; r < n; r++) {
    // insert it n times
    for (auto i = 0; i < n; i++) {
      a.insert(a.end(), M[r].begin(), M[r].end());
    }
  }
  
  // make b = (col(N,1), col(N,2),..., col(N,n)) x n times
  vector<LFr> b;
  // insert n times all cols
  for (auto i = 0; i < n; i++) {
    for (auto c = 0; c < n; c++) {
      b.insert(b.end(), N[c].begin(), N[c].end());
    }
  }
    
  MYREQUIRE(a.size() == n*n*n);
  MYREQUIRE(b.size() == n*n*n);
  
  // c is (M_{i,k} * N_{k,j}) for all possible (i,j,k)
  auto c = hadamard(a,b);
  
  // mk u from chunks of summands in c
  vector<LFr> u(n*n);
  for (auto i = 0; i < n*n; i++) {
    auto chunk_beg = c.begin()+i*n;
    u[i] = accumulate(chunk_beg, chunk_beg+n, LFr::zero());
  }
  // append M and N
  for (auto i = 0; i < n; i++) {
    u.insert(u.end(), M[i].begin(), M[i].end());
  }
  // insert col-by-col
  for (auto i = 0; i < n; i++) {
    u.insert(u.end(), N[i].begin(), N[i].end());
  }
  
  // check code
  /*
  auto check_inp_stdv = cputil::concat3(a,b,c);
  check_inp_stdv.insert(check_inp_stdv.end(), u.begin(), u.end());
  ScalarVec check_inp = mkColVector(check_inp_stdv.size(), check_inp_stdv);
  ScalarMat res = F*check_inp;
  cout << res << endl;
   */
  
  vector<LFr> tgtV; // NB: Our tgtV is all of zeros

  // keygen, proof and verification below
  
  ACRel pRel(nMulGates, F.rows(), F.cols(), u.size(), trsF, tgtV);
  
  ACPIn proverInput;
  ACVIn verifInput;

  auto pBm = make_shared<Benchmark>();
  InterpCommScheme ics;
  CPAC snarkAC(&ics);
 
  
  snarkAC.setBenchmark(pBm, "cpAC");
  double hadProveT, veqProveT;
  double hadVerT, veqVerT, comVerT;

  // Keygen
  cout << "Calling Keygen..." << endl;
  auto crs = snarkAC.keygen(&pRel); 
  print_time_string_onerep("Commitment KG Time", "keygen_comm", snarkAC);
  print_time_string_onerep("Hadamard KG Time", "keygen_had", snarkAC);
  print_time_string_onerep("Prep VEq KG Time", "keygen_prep_veq", snarkAC);
  print_time_string_onerep("Subspace KG Time", "keygen_ss", snarkAC);
  print_time_string_onerep("AuxGen KG Time", "keygen_aux", snarkAC);



  // Pack inputs
  snarkAC.prepareInputs(a, b, c, u, proverInput, verifInput);
  print_time_string_onerep("Commitment Prove Time", "commit", snarkAC);

  // NB: Prev operation must be done after commitment is ready

  // Proving
  ACPf *pf;
  auto bmProveFn = [&] {
    cout << "Proving..." << endl;
    pf = snarkAC.prove(crs, proverInput);
  };
  snarkAC.runAndAverage2(
	bmProveFn, 
	"prove_had", hadProveT,
	"prove_veq", veqProveT,
	NREPS);
  print_time("Hadamard Prove Time", hadProveT);
  print_time("VEq Prove Time", veqProveT);
  // ss provetime

  // Verifying
  auto bmVerifFn = [&] {
    cout <<  "Verifying..." << endl;
    bool isGdPf = snarkAC.verify(crs, verifInput, pf);
  };
  snarkAC.runAndAverage3(
	bmVerifFn,
  "verify_commit", comVerT,
	"verify_had", hadVerT,
	"verify_veq", veqVerT,
	NREPS);
  print_time("Commitment Verification Time", comVerT);

  print_time("Hadamard Verification Time", hadVerT);
  print_time("VEq Verification Time", veqVerT);

 }


int main()
{
  default_ec_pp::init_public_params();
  
  // n must be a power of 2
  auto lmin = 2;
  auto lmax = 8;
  for (auto l = lmin; l <= lmax; l++) {
    auto n = 1 << l;
    cout << "Measuring Time for n = " << n << endl;
    benchmark_bp_circuit(n);
  }
  
}
