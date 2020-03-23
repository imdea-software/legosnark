#ifndef POLYTOOLS_H
#define POLYTOOLS_H

#include "globl.h"
#include "commit.h"
#include <iostream>
using namespace std;


using PolyTField = In;
using Scalars = vector<PolyTField>;


class PolyT
{
protected:
  PolyT(size_t D) {
    vRepr.resize(D+1);
    fill(vRepr.begin(), vRepr.end(), PolyTField::zero());
  }


public:
    Scalars vRepr;

    PolyT(const Scalars &_v) : vRepr(_v) {

    }

    PolyT() : PolyT(0) { }

    // return zero poly of certain deg
    static PolyT zero(size_t D)
    {
      return PolyT(D);
    }

    // simply the poly "x"
    static PolyT x()
    {
      return PolyT({PolyTField::zero(), PolyTField::one()});
    }

    static PolyT one()
    {
      return PolyT({PolyTField::one()});
    }

  static PolyT one_minus_x()
  {
    return PolyT({PolyTField::one(), -PolyTField::one()});
  }

    PolyT mul(const PolyT &p) const {
      auto this_d = static_deg();
      auto p_d = p.static_deg();
      PolyT res = PolyT::zero(this_d+p_d);
      for (auto i = 0; i <= this_d; i++) {
        for (auto j = 0; j <= p_d; j++) {
          res.vRepr[i+j] = res.vRepr[i+j] + vRepr[i]*p.vRepr[j];
        }
      }
      return res;
    }

    PolyT mul(PolyTField s) const {
      auto scaleFn = [s](PolyTField x) {return x*s; };
      Scalars resRepr = cputil::map<PolyTField, PolyTField>
        (vRepr, scaleFn);
      return PolyT(resRepr);
    }

    PolyT add(const PolyT &p) {
      auto this_d = static_deg();
      auto p_d = p.static_deg();
      auto D = max(this_d, p_d);
      PolyT res = PolyT::zero(D);
      for (auto i = 0; i <= D; i++) {
        if (i <= this_d)
          res.vRepr[i] = res.vRepr[i] + vRepr[i];
        if (i <= p_d)
          res.vRepr[i] = res.vRepr[i] + p.vRepr[i];
      }
      return res;
    }

    // return size of vRepr-1
    size_t static_deg() const {
      return vRepr.size()-1;
    }

    PolyTField eval(const PolyTField pt) const {
      PolyTField res = PolyTField::zero();
      PolyTField ptPow = PolyTField::one();

      for (auto coeff : vRepr ) {
        res = res + ptPow*coeff;
        ptPow = ptPow * pt;
      }
      return res;
    }


    static Comm evalAsPolyOn(const Comms &comms, const PolyTField &pt)
    {
        auto rslt = comms[0];
        auto scaling = pt;
        for (auto i = 1; i < comms.size(); i++) {
            rslt = rslt + comms[i]*scaling;
            scaling = scaling*pt;
        }
        return rslt;
    }

    static CommOut evalAsPolyOn(const CommOuts &cmouts, const PolyTField &pt)
    {
        auto rslt = cmouts[0];
        auto scaling = pt;
        for (auto i = 1; i < cmouts.size(); i++) {
            rslt = rslt + cmouts[i]*scaling;
            scaling = scaling*pt;
        }
        return rslt;
    }

    CommOuts commit(CommScheme *comScm) const
    {
        auto mkComLambda =  [&comScm](auto coeff) {
            return comScm->commit(coeff);
        };

        return cputil::map<PolyTField,CommOut>(vRepr, mkComLambda);
    }

};


// NB: We assume this is a multilinear polynomial

class MultiVPolyT {
public:
  Scalars vRepr;
  size_t nVars;

  MultiVPolyT(const Scalars &_vRepr, size_t _nVars) : vRepr(_vRepr), nVars(_nVars)
  {
    assert((1 << nVars) == vRepr.size());
  }

  MultiVPolyT() {}


  inline static bool isVarInSet(unsigned i, uint64 S, int d)
  {
    auto e = d-i-1;
    return  S & (1 << e);
  }

  // set S of variables represented as bit-string
  inline PolyTField getCoeff(uint64 S) const
  {
    return vRepr[S];
  }

  static const Scalars mkBeta(uint64 S, int d)
  {
    auto eq_bit = [](bool b) {
      PolyTField bfld = b ? PolyTField::one() : PolyTField::zero();
      return Scalars {(bfld+bfld) - PolyTField::one(), PolyTField::one() - bfld} ;
    };

    vector<Scalars> eqs(d);
    for (auto k = 0; k < d; k++) {
      eqs[k] = eq_bit(isVarInSet(k, S, d));
    }

    Scalars beta(1 << d);
    Scalars hlpBuff(1 << d); // helper buffer
    // init
    beta[0] = eqs[d-1][0];
    beta[1] = eqs[d-1][1];
    for (auto k = d-2; k >= 0; k--) {
      auto cur_beta_sz = 1 << (d-k-1);
      for (auto j = 0; j < cur_beta_sz; j++) {
        // first half
        hlpBuff[j] = beta[j]*eqs[k][0];
        // second half
        hlpBuff[cur_beta_sz+j] = beta[j]*eqs[k][1];
      }
      swap(beta, hlpBuff);
    }
    return beta;
  }

  static const vector<Scalars> mkBetas(int d)
  {
    const uint64 N  = 1 << d;
    vector<Scalars> betas(N);
    for (uint64 S = 0; S < N; S++) {
      betas[S] = mkBeta(S, d);
    }

    return betas;
  }

  // evaluates MLE of v on r
  static PolyTField evalMLE(const Scalars &v, const Scalars &r)
  {
    auto d = r.size();
    uint64 N = v.size();
    assert(N == 1 << d);

    // compute all betas monomials of mle
    Scalars products(N);
    products[0] = PolyTField::one();

    auto f1 = PolyTField::one();

    size_t idx = 1;
    for (auto i = 0; i < d; i++) {
      uint64 pBound = 1 << i;
      for (uint64 p = 0; p < pBound; p++) {
        products[p+idx] = products[p]*r[i];
        products[p] = products[p]*(f1-r[i]);
      }
      idx += 1 << i;
    }

    auto out = PolyTField::zero();
    for (uint64 p = 0; p < v.size(); p++) {
      out = out + v[p]*products[p];
    }
    return out;
  }

  Scalars getVRepr() const
  {
    return vRepr;
  }

  PolyTField sumOverAllBinValues() const
  {
    Scalars A(nVars); // assignment
    return sumOverAllBinValues_backtrack(A, 0);
  }

  PolyTField sumOverAllBinValues_backtrack(Scalars &A, int curLen) const
  {
    if (curLen < nVars) {
      // increase assignment depth first with 0, then with 1
      A[curLen] = LFr::zero();
      auto out0 = sumOverAllBinValues_backtrack(A, curLen+1);
      A[curLen] = LFr::zero();
      auto out1 = sumOverAllBinValues_backtrack(A, curLen+1);

      return out0+out1;
    }

    // curLen == nVars: let's return evaluation for current assignment
    return evalOn(A);
  }

  enum class WhichHalf {NoFirstVar, YesFirstVar};
  void initPolyFromHalfCoeffs(MultiVPolyT &p, WhichHalf wh) const
  {
    p.nVars = nVars-1;

    Scalars::const_iterator fst, lst;
    if (wh == WhichHalf::NoFirstVar) {
      // We take first half of vRepr
      fst = begin(vRepr);
      lst = begin(vRepr) + (1 << (nVars-1));
    } else { // YesFirstVar
      // We take second half of vRepr
      fst = begin(vRepr) + (1 << (nVars-1));
      lst = end(vRepr);
    }

    p.vRepr = Scalars(fst, lst);
  }

  MultiVPolyT getAllPossibleAssigmentsPolyExceptFirstVar() const
  {
    MultiVPolyT pY, pN;
    initPolyFromHalfCoeffs(pN, WhichHalf::NoFirstVar);
    initPolyFromHalfCoeffs(pY, WhichHalf::YesFirstVar);

    auto vN = pN.sumOverAllBinValues();
    auto vY = pY.sumOverAllBinValues();

    MultiVPolyT outPoly({vN, vY}, 1);
    return outPoly;
  }

  // return p(pt, x2...xd) (se p supporta d variabili)
  void initPartialEvaluationPoly(MultiVPolyT &pOut, PolyTField pt) const {
    pOut.nVars = nVars-1;
    pOut.vRepr.resize(1 << (nVars-1));

    Scalars::const_iterator fst, mid, lst;
    fst = begin(vRepr);
    mid = begin(vRepr) + (1 << (nVars-1));
    lst = end(vRepr);

    // right part refers to that has "x0 on"
    auto i = 0;
    for (
            auto itLeft = fst, itRight = mid;
            itLeft != mid;
            itLeft++, itRight++, i++
            ) {
      pOut.vRepr[i] = (*itLeft) + (*itRight) * pt;
    }

  }


   // [l,...r)
   static PolyTField eval_impl(size_t j, uint64 l, uint64 r, const Scalars &pts, const Scalars &coeffs, size_t d)
  {
    if (j == d) {
      assert(l == r-1);
      return coeffs[l];
    }
    auto mid = l+(r-l)/2;

    PolyTField evalWithout = eval_impl(j+1, l, mid, pts, coeffs, d);
    PolyTField evalWith = eval_impl(j+1, mid, r, pts, coeffs, d);

    return pts[j]*evalWith + evalWithout;
  }

  PolyTField evalOn(const Scalars &pts) const {
    return eval_impl(0, 0, 1 << nVars, pts, vRepr, nVars);
  }

  PolyTField evalMonomial(uint64 S, const Scalars &pts) const
  {
    return getCoeff(S)*evalLiteralMonomial(S, pts);
  }

  // evaluation of monomial discarding coefficient
  PolyTField evalLiteralMonomial(uint64 S, const Scalars &pts) const
  {
    PolyTField out = LFr::one();
    for (auto i = 0; i < nVars; i++) {
      if (isVarInSet(i, S, nVars)) {
        out *= pts[i];
      }
    }
    return out;
  }

};

#endif