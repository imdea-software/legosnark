#include "lipmaa.h"

#include <algorithm>
using std::multiplies;
using std::transform;

#include <iostream>
using namespace std;

void LGlobalKeygen(long n, InterpCommScheme &ics, CPHadL &cphadl)
{
	Interpolator interp(n);
	// Trapdoors
   IScalar chi = IScalar::random_element();
   IScalar gamma = IScalar::random_element();

	ics.keygen(n, interp, chi, gamma);
	
	cphadl.keygen(n, interp, chi, gamma);	
}

CommOut InterpCommScheme::commit(const IScalars &v)
{
	auto r = IScalar::random_element();
	startBenchmark("commit");
	LG1 c = r*key.zg1 + multiExpMA<LG1>(key.lg1, v); // commitment: multiexp of lg1[i]-s on v[i]-s (+ randomness)
	KCT kc = r*key.gammazg2 + multiExpMA<LG2>(key.gammalg2, v); // XXX: can be optimized with kc_method
	stopBenchmark("commit");
	return CommOut(Comm(c, kc), r, v);
}


bool InterpCommScheme::verify(const Comm &c)
{
	GT<def_ec> lhs =
		def_ec::reduced_pairing(c.c, key.gammazg2);	
	GT<def_ec> rhs =	
		def_ec::reduced_pairing(key.zg1, c.kc);

	return lhs == rhs;
}

void check_correctness(const IScalars &coeffs, const IScalars &pts)
{
	// Compute poly on \omega^1
	auto omega = get_root_of_unity<IScalar>(pts.size());
	auto evalOnOmega = IScalar::zero();

	auto tmp = IScalar::one();
	for (auto i = 0; i < coeffs.size(); i++) {
		evalOnOmega += tmp*coeffs[i];
		tmp *= omega;

	}

	MYREQUIRE(pts[1] == evalOnOmega);

}

void mkZ(long n, IScalars &p, const IScalar r)
{
	p.resize(n+1);
	fill(p.begin(), p.end(), IScalar::zero());
	p[n] = r;
	p[0] = -r;
}

void divisionFast(long n, domain_ptr domain, IScalars &piCoeffs, const IScalars &QCoeffs, const IScalars &ZCoeffs)
{
	// identity poly
	IScalars oneCoeffs(n, IScalar::zero());
	oneCoeffs[0] = IScalar::one();

	IScalars onePts = oneCoeffs; // XXX: to optimize
	domain->cosetFFT(onePts, IScalar::multiplicative_generator);

	IScalars ZinvPts = onePts; // XXX: to optimize 
	domain->divide_by_Z_on_coset(ZinvPts);

	IScalars ZinvCoeffs = ZinvPts;
	domain->icosetFFT(ZinvCoeffs, IScalar::multiplicative_generator);

	// Just a check
	IScalars tstPoly;
	_polynomial_multiplication(tstPoly, ZinvCoeffs, ZCoeffs);
	cout << "TstPoly : ";
	print_poly(tstPoly);
	cout << endl;

	_polynomial_multiplication(piCoeffs, QCoeffs, ZinvCoeffs);
}


HadLPf CPHadL::prove(const CommOut &aCOut, const CommOut &bCOut, const CommOut &cCOut)
{
	using FieldT = IScalar;
	
	startBenchmark("prove");

	

	auto domain = get_evaluation_domain<IScalar>(n);
	std::vector<FieldT> coefficients_for_H(domain->m+1, FieldT::zero());
      {

	// Point-form of a,b and c
	auto &aPts = aCOut.xs;
	auto &bPts = bCOut.xs;
	auto &cPts = cCOut.xs;

	// randomness
	auto d1 = aCOut.r;
	auto d2 = bCOut.r;
	auto d3 = cCOut.r;
	
    std::vector<FieldT> aA(domain->m, FieldT::zero()), aB(domain->m, FieldT::zero());

    for (size_t i = 0; i < n; ++i)
    {
        aA[i] = aPts[i];
        aB[i] = aPts[i];
    }
    domain->iFFT(aA);
    domain->iFFT(aB);
#ifdef MULTICORE
#pragma omp parallel for
#endif
    /* add coefficients of the polynomial (d2*A + d1*B - d3) + d1*d2*Z */
    for (size_t i = 0; i < domain->m; ++i)
    {
        coefficients_for_H[i] = d2*aA[i] + d1*aB[i];
    }
    coefficients_for_H[0] -= d3;
    domain->add_poly_Z(d1*d2, coefficients_for_H);

    domain->cosetFFT(aA, FieldT::multiplicative_generator);

    domain->cosetFFT(aB, FieldT::multiplicative_generator);
    std::vector<FieldT> &H_tmp = aA; // can overwrite aA because it is not used later
#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        H_tmp[i] = aA[i]*aB[i];
    }
    std::vector<FieldT>().swap(aB); // destroy aB

    std::vector<FieldT> aC(domain->m, FieldT::zero());
    for (size_t i = 0; i < n; ++i)
    {
        aC[i] += cPts[i];
    }

    domain->iFFT(aC);
    domain->cosetFFT(aC, FieldT::multiplicative_generator);

#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        H_tmp[i] = (H_tmp[i]-aC[i]);
    }

    domain->divide_by_Z_on_coset(H_tmp);

    domain->icosetFFT(H_tmp, FieldT::multiplicative_generator);

#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < domain->m; ++i)
    {
        coefficients_for_H[i] += H_tmp[i];
    }
  }
    
	fmt::print("SZ of multiexp in Hadamard: {}\n", key.chipowsg1.size());
	auto ret = multiExpMA<LG1>(key.chipowsg1, coefficients_for_H); 

	stopBenchmark("prove");
	
	return ret; 
}

bool CPHadL::verify(HadLPf pf, const Comm &ca, const Comm &cb, const Comm &cc)
{
	startBenchmark("verify");
	
	auto cck_precomp = def_ec::precompute_G2(cc.kc);
	auto pf_precomp = def_ec::precompute_G1(pf);
	
	auto rhs = def_ec::double_miller_loop(
		key.g1_precomp, cck_precomp, 
		pf_precomp, key.gammazg2_precomp);
	
	auto ca_precomp = def_ec::precompute_G1(ca.c);
	auto cbk_precomp = def_ec::precompute_G2(cb.kc);
	auto lhs = def_ec::miller_loop(ca_precomp, cbk_precomp);

	
	auto shouldBeOne = def_ec::final_exponentiation(lhs*rhs.unitary_inverse());
	bool isGood =  (shouldBeOne == LGT::one());
	stopBenchmark("verify");
	return isGood;
}
