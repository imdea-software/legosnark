#ifndef MLE_H
#define MLE_H

#include "globl.h"
#include "polytools.h"


LFr eqbit(LFr x, LFr r);
LFr eqbit(bool b, LFr r);
PolyT eqbit_poly(bool b);
PolyT eqbit_poly(LFr r);
In evalBetaOnPoint(const Ins &rho, const Ins &s);


class DPBeta {
public:
  size_t d;
  Ins rho;
  Ins rhoInvs;

  DPBeta(size_t _d, const Ins &_rho) : d(_d), rho(_rho)
  {
    assert(d == rho.size());
    precomputeAll();
  }


  In getBetaPre(size_t j) const // j in {-1,...,d-2}
  {
    if(j == -1) {
      return LFr::one();
    } else { // j >= 0
      return beta_pre_rho[j];
    }
    // should never get here
  }

  virtual void pushRandomness(In r, size_t j) {
    // Expects 0 <= j < d-1

    // update beta_pre_rho
    beta_pre_rho[j] = getBetaPre(j-1) * eqbit(r, rho[j]);

    // update beta_suff_rho
    if (j >= d-2) { // no need to update at last step
      return;
    }
    swap(beta_suff_rho_cur, beta_suff_rho_old);
    uint64 pBound = 1 << (d-j-2);
    for (uint64 p = 0; p < pBound; p++) {
      auto old_p = pBound + p; // 1 concat p
      beta_suff_rho_cur[p] = beta_suff_rho_old[old_p]*rhoInvs[j+1];
    }
    cur_suff_j++;

  }

  In getBetaSuff(size_t j, uint64 p) const
  {
    // Expects 1 <= j <= d
    if (j > d-1) {
      return In::one();
    }

    // check we are not asking before it's time
    if (j > cur_suff_j) {
      throw std::runtime_error("Requested beta_suffix is not available yet");
    }

    // j \in {1,...,d-1}
    return beta_suff_rho_cur[p];


  }

  virtual PolyT getBetaPoly(size_t j, uint64 p) const
  {
    // Expects 0 <= j < d-1
    auto pre = getBetaPre(j-1);
    auto suff = getBetaSuff(j+1, p);
    return eqbit_poly(rho[j]).mul(pre*suff);
  }

  // beta_pre_rho[j] = eqbit(r[0], rho[0]) \cdot ... \cdot eqbit(r[j], rho[j])
  Ins beta_pre_rho;

  // beta_suff_rho[j][p_j,...,p_{d-1}] = eqbit(p_j, rho[j]) \cdot ... \cdot eqbit(p_{d-1}, rho[d-1])
  Ins beta_suff_rho_cur, beta_suff_rho_old;

  size_t cur_suff_j; // debug variable

  // at the end dst[p] = eq(p, r) for all p. tmp is used as auxiliary
  static void compute_eq_tbl(size_t d, Ins &dst, Ins &tmp, const Ins &r)
  {
    dst[0] = eqbit(false, r[0]);
    dst[1] = eqbit(true, r[0]);

    for (auto j = 1; j < d; j++) {
      for (uint64 p = 0; p < 1 << (j+1); p++) {
        bool msb = (p >= (1 << j)); // most significant bit of p
        tmp[p] = eqbit(msb, r[j])*dst[p >> 1];
      }
      swap(tmp, dst);
    }
  }

  virtual In evalOnPoint(const Ins &rho, const Ins &s)
  {
    return evalBetaOnPoint(rho, s);
  }

protected:

  DPBeta() {}

  void cacheInvsOfRho()
  {
    auto invFn = [](In x) { return x.inverse(); };
    rhoInvs = cputil::map<In,In>(rho, invFn);
  }

  void precomputeAll()
  {
    cacheInvsOfRho();

    // init beta_pre
    beta_pre_rho.resize(d+1); // XXX: check this size is good // NB:actually we just need last value
    //

    // init_beta_suff
    beta_suff_rho_cur.resize(1<< d);
    beta_suff_rho_old.resize(1<< d);

    compute_eq_tbl(d, beta_suff_rho_old, beta_suff_rho_cur, rho);

    // here beta_suff_rho_old has everything. But we remove rho[0] from _cur to be ready to be used
    uint64 pBound = 1 << (d-1);
    for (uint64 p = 0; p < pBound; p++) {
      auto old_p = p + pBound; // 1 concat p
      beta_suff_rho_cur[p] = beta_suff_rho_old[old_p]*rhoInvs[0]; // XXX: are we sure we are doing old_p the right way? (maybe it's computed vice versa with the lower bit)
    }

    cur_suff_j = 1;

  }
};

class DPBetaDummy : public DPBeta
{
public:
  DPBetaDummy() {}

  virtual PolyT getBetaPoly(size_t j, uint64 p) const override {
    return PolyT::one(); // return identity polynomial
  }


  virtual void pushRandomness(In r, size_t j) override {
    // do nothing
  }

  virtual In evalOnPoint(const Ins &rho, const Ins &s) override
  {
    return In::one();
  }

  };


class DPMle {
protected:
  //Ins v;
  size_t d;
  uint64 n;

  Ins curVTable, oldVTable;
  Ins v; // v contains the "original" vector of points we are doing mle on
  // NB: v is not the whole matrix in the matrix version of DPMle but a processed variant

public:
  DPMle(size_t _d, uint64 _n) : d(_d), n(_n)
  {
    oldVTable.resize(1 << d);
    curVTable.resize(1 << d);
    v.resize(1 << d);
    fill(curVTable.begin(), curVTable.end(), In::zero());
    fill(v.begin(), v.end(), In::zero());
  }

  DPMle(size_t _d, uint64 _n, const Ins &_v) : d(_d), n(_n), curVTable(_v), v(_v)
  {
    oldVTable.resize(1 << d);
  }

  const Ins getV() const {
    return v;
  }

  void pushRandomness(In r, size_t j) {
    // We have 0 <= j <= d-1
    swap(curVTable, oldVTable);
    // int64 pStart = 1 << (j+1);
    uint64 pBound = 1 << (d-j-1);
    for (uint64 p = 0; p < pBound; p++ ) {
      auto p0 = p; // 0 concat p
      auto p1 = p0 + pBound; // 1 concat p
      curVTable[p] = oldVTable[p0]*eqbit(false, r) +
                     oldVTable[p1]*eqbit(true, r);
    }
  }

  In getVTable(size_t j, uint64 p) const {
    // we ignore j
    return curVTable[p];
  }

  PolyT getMLEPoly(size_t j, uint64 p) const
  {
    // We have 0 <= j <= d-1
    auto p0 = p; // 0 concat p
    auto v0 = eqbit_poly(0).mul(getVTable(j, p0));
    auto p1 = p0 + (1 << (d-j-1)); // 1 concat p
    auto v1 = eqbit_poly(1).mul(getVTable(j, p1));
    return v0.add(v1);
  }
};


/*
 * // XXX: Below should be done with the first d bits or not according to the situation
 * Like DPMle, but allows for preprocessing of vectors for CPMatMult:
 * \tilde{A}_\rho(i) =
 * = \sum_{p \in {0,1}^{2*_d} } eq(p1...p_d , rho) eq(p_d+1...p_2_d, i) A_p
 * = \sum_r eq(r, i) (\sum_l A_{l||r} eq(l, rho)
 * = \sum_r eq(r, i) v_r
 * */
class DPMatrixMle: public DPMle
{
public:
  // A is a vectorized matrix of size _n x _n, with _n = 2^_d
  DPMatrixMle(size_t _d, uint64 _n, const Ins &_A, const Ins &rho) :
    DPMle(_d, _n)
  {
    const uint64 N = _n*_n;
    //const size_t D = 2*_d;

    Ins eqTbl(_n), tmp(_n);
    DPBeta::compute_eq_tbl(_d, eqTbl, tmp, rho);

    // we preprocess the vector scaling every element _A[p]

    for (uint64 r = 0; r < _n; r++) {
      for (uint64 l = 0; l < _n; l++) {
        auto p = (l << _d) + r;
        auto inc = _A[p] * eqTbl[l];
        v[r] = curVTable[r] = curVTable[r] + inc;
      }
    }

  }

};


#endif