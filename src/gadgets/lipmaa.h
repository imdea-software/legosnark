#ifndef LIPMAA_H
#define LIPMAA_H

#include "interp.h"
#include "commit.h"



template <typename T>
void print_poly(const vector<T> &c) {
  for (size_t i = 0; i < c.size(); i++)
    {
      unsigned long coefficient = c[i].as_ulong();

      if (i == 0) std::cout << coefficient << " + ";
      else if (i < c.size()-1) std::cout << coefficient << "x^" << i << " + ";
      else std::cout << coefficient << "x^" << i << std::endl;
    }
}


class CPHadL;

struct InterpCommKey {
  LG1 zg1;
  vector<LG1> lg1;
  
  LFr z;
  vector<LFr> l;

  LG2 gammazg2;
  vector<LG2> gammalg2;

};

class InterpCommScheme : public CommScheme {
public:
  void keygen(long _n, Interpolator &interp, IScalar chi, IScalar gamma)
  {
		n = _n;
		key.l = interp.getAllLagrangianPolys(chi);

		key.lg1 = interp.mkG1Exp(key.l);
		
		key.z = interp.mkZ(chi);
		key.zg1 = key.z*LG1::one();

		key.gammazg2 = gamma*interp.mkZ(chi)*LG2::one();

		auto mulByGamma = [gamma](const IScalar &x) { return x*gamma; };
		auto lEvalOnChiGamma = cputil::map<IScalar,IScalar>(key.l, mulByGamma);
		key.gammalg2 = interp.mkG2Exp(lEvalOnChiGamma);
  }

  friend void LGlobalKeygen(long n, InterpCommScheme &ics, CPHadL &cphadl);
  InterpCommKey key;


  InterpCommScheme() : CommScheme()
  {
       // Nothing to be done here
  }


  CommOut commit(const IScalars &v) override;
  
  bool verify(const Comm &c);
  
  void print_key_size() {
    fmt::print("Commitment key Size: {} G1 + {} G2\n", key.lg1.size()+1, key.gammalg2.size()+1);
  }

};

using HadLPf = LG1;

struct HadLKey {
  LG2 gammazg2;
  vector<LG1> chipowsg1;
  libff::G1_precomp<def_ec> g1_precomp;
  libff::G2_precomp<def_ec> gammazg2_precomp;
};

class CPHadL : public Benchmarkable
{
public:
 void keygen(long _n, Interpolator &interp, IScalar chi, IScalar gamma)
  {
		n = _n;
	
		key.gammazg2 = gamma*interp.mkZ(chi)*LG2::one();
		
		startBenchmark("keygen");

		key.g1_precomp = def_ec::precompute_G1(LG1::one());
		key.gammazg2_precomp = def_ec::precompute_G2(key.gammazg2);
		
		vector<IScalar> chiPows(n);
		chiPows[0] = chi; // this will get shifted later at pos. 1
		for (auto i = 1; i < n; i++) {
			chiPows[i] = chi*chiPows[i-1];
		}
			
		key.chipowsg1 = interp.mkG1Exp(chiPows);
		stopBenchmark("keygen");

		key.chipowsg1.insert(key.chipowsg1.begin(), LG1::one()); // NB: Would be O(1) with simple optimization
  }

	friend void LGlobalKeygen(long n, InterpCommScheme &ics, CPHadL &cphadl);

	HadLKey key;
	long n;

	CPHadL() {}
  
  void print_key_size() {
    fmt::print("Had PK Size: {} G1 \n", key.chipowsg1.size());
    fmt::print("Had VK Size: {} G2\n", 1);
  }
	
	HadLPf prove(const CommOut &aCOut, const CommOut &bCOut, const CommOut &cCOut);
	bool verify(HadLPf, const Comm&, const Comm&, const Comm&);
};


#endif
