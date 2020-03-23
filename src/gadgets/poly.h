#ifndef CP_POLY
#define CP_POLY

#include "globl.h"
#include "polytools.h"


struct PolyPf
{
    vector<LG1> witness;
    vector<LG1> witnessa;

    // size in group elements
    size_t getSize() const
    {
      return witness.size() + witnessa.size();
    }
};


class CPPoly : public Benchmarkable {

  CommScheme *cmScm;
public:
    CPPoly(CommScheme *_cmScm) : cmScm(_cmScm) { }

    void keygen(int d) { } // XXX: For now separated from other commitment's kg

    // NB: For now we basically commit to the vector
    CommOut commitPoly(const Scalars &v)  const {
      return cmScm->commit(v);
    }

    bool checkCommit(const Comm &cm) const {
      return simple_pairing_check(cm.c, LG2::one(), LG1::one(), cm.kc);
    }


    void computeAnswer(CommOut &cmoutAns, const Ins &input, const Ins &v) const
    {
        auto ans = MultiVPolyT::evalMLE(v, input);
        cmoutAns = cmScm->commit(ans);
    }

    void prove(const Ins &v, const CommOut &cmoutAns, const Ins &input, PolyPf &pf) const
    {
        auto &r = input;
        size_t d = r.size();
        uint64 N = v.size();

        // witness coefficients
        Scalars w_coeffs(1 << d);
        Scalars tmp_v = v;

        size_t start = 0;
        for (auto i = 0; i < d; i++) {
            uint64 pBound = 1 << (d-i-1);
            // p is sort of an offset
            for (uint64 p = 0; p < pBound; p++) {
                uint64 p0 = p << 1;
                uint64 p1 = (p << 1) + 1;
                w_coeffs[start+p] = -tmp_v[p0] + tmp_v[p1];
                tmp_v[p] = -tmp_v[p0]*(r[i]-1) + tmp_v[p1]*r[i];
            }
            tmp_v.resize(pBound);
            start += pBound;
        }

        pf.witness.resize(d);
        pf.witnessa.resize(d);

        // NB: Putting g1 as bases; benchmark purposes only.
        auto g1s = cmScm->getBases1();

        // make multiexps
        start = 0;
        for (auto i = 0; i < d; i++) {
            uint64 pBound = 1 << (d-i-1);
            Scalars tmp_e(pBound);
            for (uint64 p = 0; p < pBound; p++) {
                tmp_e[p] = w_coeffs[start+p];
            }
            pf.witness[i] = multiExpMA<LG1>(g1s, tmp_e); // NB: bases for benchmarking purposes only
            if (i != 0) {
                pf.witnessa[i] = multiExpMA<LG1>(g1s, tmp_e); // NB: bases for benchmarking purposes only
            }
            start += pBound;
        }


    }

    bool verify(const  Comm &cmPoly, const Comm& cmAns, const Ins &pts, const PolyPf &pf) const
    {
      auto d = pts.size();

      auto g2_precomp = def_ec::precompute_G2(LG2::one());
      auto g2a_precomp = def_ec::precompute_G2(LG2::one()); //NB: Should use g^a. Benchmarking purposes only

      bool isGd = true;
      using preG1T = decltype(def_ec::precompute_G1(LG1::one()));
      vector<preG1T> wl(d);

      /* This is check commit basically */
      for (auto i = 1; i < d; i++) {
          wl[i] = def_ec::precompute_G1(pf.witness[i]);
          auto wr = def_ec::precompute_G1(pf.witnessa[i]);
          auto lhs = def_ec::miller_loop(wl[i], g2_precomp);
          auto rhs = def_ec::miller_loop(wr, g2a_precomp);
          auto out = def_ec::final_exponentiation(lhs*rhs.unitary_inverse());
          isGd = isGd && (out == LGT::one());
       }


      auto acc = libff::Fqk<def_ec>::one();
      for (auto i = 0; i < pf.witness.size(); i++) {
        auto base2_precomp = def_ec::precompute_G2(pts[i]*LG2::one());
        acc = acc*def_ec::miller_loop(wl[i], base2_precomp);
      }
      auto cmPolyAns_pre = def_ec::precompute_G1(cmPoly.c-cmAns.c);
      acc = acc.unitary_inverse() * def_ec::miller_loop(cmPolyAns_pre, g2_precomp);
      auto out = def_ec::final_exponentiation(acc);
      isGd = isGd && (out == LGT::one());

      return isGd;
    }

};

#endif