#ifndef COMMIT_H
#define COMMIT_H

#include "globl.h"
#include "util.h"
#include "benchmark.h"

#include <vector>
using std::vector;

using KCT = LG2; // Knowledge Curve

/* == Basic types for inputs and randomness == */

using CommRand = LFr;
class Comm;
using Comms = vector<Comm>;
using In = LFr;
using Ins = vector<LFr>;


/* == Commitment-related classes == */

// TODO: This should be an abstract class
class Comm
{
public:
  // Actual commitment and its knowledge counterpart
  LG1 c;
  KCT kc;

  Comm() : Comm(LG1::zero(), KCT::zero()) { }
  Comm(LG1 _c, KCT _kc) : c(_c), kc(_kc) { }

  friend Comm operator+(const Comm& a, const Comm& b)  {
    return Comm(a.c+b.c, a.kc+b.kc);
  }
  friend Comm operator-(const Comm& a, const Comm& b) {
    return Comm(a.c-b.c, a.kc-b.kc);
  }

   Comm operator*(const CommRand b) const {
    auto cOut = b*c;
    auto kcOut = b*kc;
    return Comm(cOut, kcOut);
  }

  void set(LG1 _c, KCT _kc) {
    c = _c;
    kc = _kc;
  }

};

struct CommOut;
using CommOuts = vector<CommOut>;


struct CommOut {
  Comm c;
  CommRand r;

  Ins xs;
  int lenXs;

  CommOut() {}

  CommOut(const Comm _c, const CommRand _r, const Ins _xs) :
    c(_c), r(_r), xs(_xs), lenXs(_xs.size())
  {

  }
  CommOut(const Comm _c, const CommRand _r, const In _x)
  : c(_c), r(_r), xs(vector<In> {_x}), lenXs(1)
  {
  }

  In val() const {
    if (lenXs != 1) {
      throw runtime_error("val() is only for commitments to single elements.");
    }
    return xs[0];
  }

  CommOut operator*(const CommRand b) const
  {
      if (lenXs != 1) {
          throw runtime_error("Operator * on CommOut is only for commitments to single elements.");
      }
      return CommOut(c*b, r*b, xs[0]*b);
  }

  friend CommOut operator+(const CommOut& a, const CommOut& b);
  friend CommOut operator-(const CommOut& a, const CommOut& b);

  static Comms toComms(const CommOuts &outs)
  {
    auto commFn = [](auto comOut) { return comOut.c;};

    return cputil::map<CommOut, Comm>(outs, commFn);
  }

  static vector<CommRand> toOpenings(const CommOuts &outs)
  {
    auto openFn = [](auto comOut) { return comOut.r;};

    return cputil::map<CommOut, CommRand>(outs, openFn);
  }

  static vector<Ins> toCommittedVals(const CommOuts &outs)
  {
    auto commValFn = [](auto comOut) { return comOut.xs;};

    return cputil::map<CommOut, Ins>(outs, commValFn);
  }
};



class CommScheme : public Benchmarkable {
public:
    long n;

    CommScheme()
    {
        // Nothing to be done here
    }

    virtual void keygen(long _n)
    {
        n = _n;

        // NB: or now no ZK
        g1s.resize(n);
        fill(g1s.begin(), g1s.end(), LG1::one());

        g2s.resize(n);
        fill(g2s.begin(), g2s.end(), LG2::one());
    }

    LG1 getBlindingH() const {
        return LG1::one(); // XXX: Should actually be computed at kg time
    }

    vector<LG1> getBases1() const {
      return g1s;
    }

    virtual CommOut commit(const Ins &v) {
        assert(g1s.size() != 0);
        //auto r = CommRand::random_element(); // XXX: Ignored
        CommRand r = CommRand::zero();

        LG1 c = multiExpMA<LG1>(g1s, v) + r*getBlindingH();
        LG2 kc = multiExpMA<LG2>(g2s, v);

        return CommOut(Comm(c, kc), r, v);
    }

    virtual CommOut commit(const In &v) {
        auto r = In::random_element(); // XXX: Ignored
        LG1 c = v*g1s[0];
        LG2 kc = v*g2s[0];

        return CommOut(Comm(c, kc), r, v);
    }


protected:
    vector<LG1> g1s;
    vector<LG2> g2s;


};



#endif
