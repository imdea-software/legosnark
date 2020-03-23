#ifndef GADGETS_SNARK_H
#define GADGETS_SNARK_H

#include "commit.h"
#include "benchmark.h"
#include "lipmaa.h"


using Rand = LFr;

// pure virtual class
template
<typename RelTyp, typename KTyp, typename PfTyp, typename PInTyp, typename VInTyp>
class Snark : public Benchmarkable
{
public:
  Snark()
  {

  }

  virtual KTyp* keygen(const RelTyp *rel) = 0;
  virtual PfTyp* prove(const KTyp *crs, const PInTyp &u) = 0;
  virtual bool verify(const KTyp *crs, const VInTyp &x, const PfTyp *pf) = 0;

protected:

  void sampleRandomElement(Rand &r)
  {
    r = LFr::random_element();
  }

  void sampleRandomVector(vector<Rand> &rs, int l)
  {
    rs.reserve(l);
    for (auto i = 0; i < l; i++) {
      Rand tmp;
      sampleRandomElement(tmp);
      rs.push_back(tmp);
    }
  }
};



// helper for CPSnark prover
struct CPPIn
{
  Ins publicSlot; // public (uncommitted input) // this is x in the paper
  CommOuts commSlot; // committed input // this is (u,\omega) in the paper

  bool hasPublicInput() { return publicSlot.size() != 0; }
};

// helper for CPSnark verifier
struct CPVIn
{
  Ins publicIn; // public (uncommitted input) // this is x in the paper
  Comms commIn; // committed input // this is (u,\omega) in the paper

  bool hasPublicInput() { return publicIn.size() != 0; }
};

// helper class
class CPInputFmt
{
public:
    static void init_no_pub(CPPIn &prvIn, CPVIn &vrfIn, CommScheme *commScm, const vector<Ins> &vecIns)
    {
      auto vecComFn = [&commScm](const Ins &v) { return commScm->commit(v); };
      CommOuts committedIns = cputil::map<Ins, CommOut>(vecIns, vecComFn);
      prvIn.commSlot = committedIns;
      vrfIn.commIn = CommOut::toComms(prvIn.commSlot);
    }

    static void init(CPPIn &prvIn, CPVIn &vrfIn, CommScheme *commScm, const vector<Ins> &vecIns, const Ins &pub)
    {
      init_no_pub(prvIn, vrfIn, commScm, vecIns);
      prvIn.publicSlot = vrfIn.publicIn = pub;
    }

};

template <typename RelTyp, typename KTyp, typename PfTyp, typename PIn, typename VIn>
class CPSnark : public Snark<RelTyp, KTyp, PfTyp, PIn, VIn> {
public:
  CommScheme *commScm;
  CommScheme *getCommScheme() const { return commScm; }

  CPSnark(CommScheme *_commScm) :
    Snark<RelTyp, KTyp, PfTyp, PIn, VIn>(),
    commScm(_commScm)
  {

  }

};



#endif
