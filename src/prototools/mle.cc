#include "mle.h"


LFr eqbit(LFr x, LFr r)
{
  auto f1 = LFr::one();
  auto f2 = LFr::one()+LFr::one();

  return x*(f2*r-f1) + f1 - r;
}

LFr eqbit(bool b, LFr r)
{
  return b ? r : LFr::one()-r;
}

PolyT eqbit_poly(bool b)
{
  return b ? PolyT::x() : PolyT::one_minus_x();
}


PolyT eqbit_poly(LFr r)
{
  auto f1 = LFr::one();
  auto f2 = LFr::one()+LFr::one();

  return PolyT({f1-r, f2*r-f1});
}


In evalBetaOnPoint(const Ins &rho, const Ins &s)
{
  In out = In::one();
  for (auto j = 0; j < rho.size(); j++) {
    out = out*eqbit(rho[j], s[j]);
  }
  return out;
}