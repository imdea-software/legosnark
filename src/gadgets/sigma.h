#ifndef SIGMA_H
#define SIGMA_H

#include "commit.h"

// TODO: Make it inherit from Snark

class ZKEqProof {
public:
  CommScheme *comScm;
  CommRand c, z;
  LG1 a;
  LG1 c0,c1;

  LG1 h; // blinding

  ZKEqProof(CommScheme *comScm, CommOut cOut0, CommOut cOut1);
  bool verify();

  size_t getSize() const {
    return 3; // does not consider field elements
  }
};

class ZKPrdProof {
public:
  CommScheme *comScm;
  CommRand c, z1, z2, z3, z4, z5;
  LG1 alpha, beta, delta;
  LG1 c0, c1, cPrd;

  LG1 h; // blinding

  // returns g^x \cdot h^r
  LG1 ghPow(In x, CommRand r) const {
    return x*LG1::one() + r*h;
  }

  ZKPrdProof(CommScheme *comScm, CommOut cOut0, CommOut cOut1, CommOut cOutPrd);
  bool verify() const;

  size_t getSize() const {
    return 6; // does not consider field elements
  }
};


#endif
