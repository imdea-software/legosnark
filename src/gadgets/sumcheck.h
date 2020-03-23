#ifndef CP_SUMCHECK_H
#define CP_SUMCHECK_H

#include "snark.h"
#include "sigma.h"
#include "poly.h"
#include "polytools.h"
#include "mle.h"
#include <iostream>
using namespace std;


using SumcheckRand = vector<CommRand>;
using SumcheckField = PolyTField;

struct SumcheckPf {
	const SumcheckRand r;

	using EqProofs = vector<shared_ptr<ZKEqProof>>;
	const EqProofs eqPfs;
 	const vector<Comms> hCom;
	const Comms polycm_g, polycm_evalg;
	const vector<PolyPf> polypf_g;

	ZKPrdProof finalPrdPf;

	SumcheckPf(
					SumcheckRand &_r, vector<Comms> &_hCom, const EqProofs &_eqPfs, const Comms &cg,
					const Comms &cevg, const vector<PolyPf> &pf_g, const ZKPrdProof &prdPf) :
						r(_r), hCom(_hCom), eqPfs(_eqPfs), polycm_g(cg),
						polycm_evalg(cevg), polypf_g(pf_g), finalPrdPf(prdPf)
	{
	}

	size_t getSize() const {
		size_t polypf_sz = 0;
		for (auto &ppf : polypf_g) {
			polypf_sz += ppf.getSize();
		}
		return eqPfs.size()*eqPfs[0]->getSize() + hCom.size()*hCom[0].size() +
			2*polycm_g.size() + 2*polycm_evalg.size() + polypf_sz + finalPrdPf.getSize();
	}

};


using SumcheckRel = size_t;

struct SumcheckKey {
  size_t d;
	// NB: d is the size of the sumcheck (# of variables) but not necessarily the log of the vectors on which we work
	// (in the matrix case d is not the log of the size of the vectors but its half)

  SumcheckKey(size_t _d) : d(_d) { }
};


class CPSumcheck : public CPSnark<SumcheckRel, SumcheckKey, SumcheckPf, CPPIn, CPVIn>
{
public:

  static const size_t n_cm_polys = 2; // number of committed polys

	CPSumcheck(CommScheme *_commScm, CPPoly *_cppoly) :
    CPSnark(_commScm), cppoly(_cppoly)
  {

}

  virtual SumcheckKey* keygen(const SumcheckRel *rel) override;
  virtual SumcheckPf* prove(const SumcheckKey *crs, const CPPIn &in) override;
  virtual bool verify(const SumcheckKey *crs, const CPVIn &, const SumcheckPf *pf) override;

	virtual shared_ptr<DPBeta> init_beta(size_t d, const Ins &rho)
	{
		return make_shared<DPBeta>(d, rho);
	}

  virtual void init_mles(vector<shared_ptr<DPMle>> &mles, size_t d, const Ins &rho, const Ins &a, const Ins &b)
	{
		mles.push_back(make_shared<DPMle>(d, 1 << d, a));
		mles.push_back(make_shared<DPMle>(d, 1 << d, b));
	}

	PolyT make_new_h_poly(size_t d, size_t j, const shared_ptr<DPBeta> beta, const vector<shared_ptr<DPMle>> &mles)
	{
		// D = 1 + # of mle poly-s // If you have beta, a and b, then it's 3
		// out_poly = empty_poly of degree D
		size_t D = CPSumcheck::n_cm_polys;
		PolyT out_poly = PolyT::zero(D);
		auto bound_p = 1 << (d-j-1);
		for (auto p = 0; p < bound_p; p++) {
			PolyT beta_poly = beta->getBetaPoly(j, p);

			// poly_p: polynomial increment depending on p
			auto poly_p = beta_poly;

			for (const shared_ptr<DPMle> mle : mles) {
				PolyT mle_poly = mle->getMLEPoly(j, p);
				poly_p = poly_p.mul(mle_poly);
			}
			// update out_poly
			out_poly = out_poly.add(poly_p);
		}
		return out_poly;
	}

  CPPoly *cppoly;
};

// each vector represents a matrix
class CPSumcheckMatrix : public CPSumcheck
{
public:
	CPSumcheckMatrix(CommScheme *_commScm, CPPoly *_cppoly)
	 : CPSumcheck(_commScm, _cppoly) {}


	virtual shared_ptr<DPBeta> init_beta(size_t d, const Ins &rho) override
	{
		return make_shared<DPBetaDummy>();
	}

	virtual void init_mles(vector<shared_ptr<DPMle>> &mles, size_t d, const Ins &rho, const Ins &a, const Ins &b) override
	{
		// a and b represent matrices but we want sumcheck only on sqr(a.size()) bits
		mles.push_back(make_shared<DPMatrixMle>(d, 1 << d, a, rho));
		mles.push_back(make_shared<DPMatrixMle>(d, 1 << d, b, rho));
	}

};



#endif
