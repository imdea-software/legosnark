#include "sigma.h"

ZKEqProof::ZKEqProof(CommScheme *_comScm, CommOut cOut0, CommOut cOut1) :
  comScm(_comScm)
{
  c0 = cOut0.c.c; // ignoring knowledge part
  c1 = cOut1.c.c;

  h = comScm->getBlindingH();

  CommRand r = CommRand::random_element(); // this is Prover's randomness; no need for RO
  a = r*h;

  c = CommRand::random_element();

  z = c*(cOut0.r-cOut1.r) + r;
}

bool ZKEqProof::verify()
{
  auto lhs = z*h;
  auto rhs = c*(c0-c1) + a;
  return lhs == rhs;
}


ZKPrdProof::ZKPrdProof(CommScheme *_comScm, CommOut cOut0, CommOut cOut1, CommOut cOutPrd) :
        comScm(_comScm)
{
  c0 = cOut0.c.c; // ignoring knowledge part
  c1 = cOut1.c.c;
  cPrd = cOutPrd.c.c;

  h = comScm->getBlindingH();

  auto x = cOut0.val();
  auto y = cOut1.val();
  auto rx = cOut0.r;
  auto ry = cOut1.r;
  auto rz = cOutPrd.r;

  vector<CommRand> bs(6); // we discard the first
  for (auto i = 1; i < bs.size(); i++) {
    bs[i] = CommRand::random_element();
  }

  c  = CommRand::random_element();

  z1 = bs[1] + c*x;
  z2 = bs[2] + c*rx;
  z3 = bs[3] + c*y;
  z4 = bs[4] + c*ry;
  z5 = bs[5] + c*(rz-rx*y);

}

bool ZKPrdProof::verify() const
{
  vector<LG1> lhs { alpha + c*c0, beta + c*c1, delta + c*cPrd };
  vector<LG1> rhs { ghPow(z1, z2), ghPow(z3, z4), ghPow(z3, z5)};

  for (auto i = 0; i < lhs.size(); i++) {
    if (lhs[i] != rhs[i]) {
      return false;
    }
  }
  return true;

}
