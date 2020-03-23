#ifndef INTERP_H
#define INTERP_H

#include "globl.h"

#include <vector>
#include <memory>
using std::vector;
using std::shared_ptr;

using IScalar = LFr;
using IScalars = vector<IScalar>;
using IPoly = vector<IScalar>;

#include <libfqfft/evaluation_domain/get_evaluation_domain.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/common/utils.hpp>

using domain_ptr = std::shared_ptr<evaluation_domain<IScalar> >;

class Interpolator {
	public:
	domain_ptr domain;
	long n;
	const size_t fldBitSz;

	
	size_t g1_exp_count;
	size_t g1_window;
	window_table<LG1> g1_table;

	size_t g2_exp_count;
	size_t g2_window;
	window_table<LG2> g2_table;
	
	template<typename T>
	void setupExp(
		const long n,
		size_t &g_exp_count,
		size_t &g_window,
		window_table<T> &g_table)
	{
		// setup multiexp
		g_exp_count = n;
		g_window = get_exp_window_size<T>(g_exp_count);
		g_table = get_window_table(fldBitSz, g_window, T::one());
	}
	
	vector<LG1> mkG1Exp(const IScalars &xs)
	{	
		assert(xs.size() <= n);
		return batch_exp(fldBitSz, g1_window, g1_table, xs);
	}

	vector<LG2> mkG2Exp(const IScalars &xs)
	{	
		assert(xs.size() <= n);
		return batch_exp(fldBitSz, g2_window, g2_table, xs);
	}
	
	Interpolator(long _n, long N) : n(_n), fldBitSz(LFr::size_in_bits()) {
		domain = get_evaluation_domain<IScalar>(n);
		setupExp<LG1>(N, g1_exp_count, g1_window, g1_table);
		setupExp<LG2>(N, g2_exp_count, g2_window, g2_table);
	}
	Interpolator(long _n) : Interpolator(_n, _n) {}
	
IScalars getAllLagrangianPolys(const IScalar tgtPt) {
		auto lagrangianPolysOnTgt = domain->evaluate_all_lagrange_polynomials(tgtPt);
		return lagrangianPolysOnTgt;
	}
	
	
	IScalar mkZ(const IScalar t)
	{
		auto Zt = domain->compute_vanishing_polynomial(t);
		return Zt;
	}

};


#endif
