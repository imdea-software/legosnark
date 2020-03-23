<h1 align="center">LegoSNARK</h1>

<div align="center">
  <strong>Composable Commit-and-Prove zkSNARKs</strong>
</div>

<div align="center">
   :envelope:  + &#x1F9F1; + :wrench:  &#8594; <code>zk</code>&#129416;
</div>

<div align="center" >
   <sub>commitments, gadgets and a framework for commit-and-prove zkSNARKs</sub>
</div>



<br/>

This codebase is part of the [LegoSNARK paper](https://eprint.iacr.org/2019/142.pdf).

**What this codebase includes:** example and benchmark implementations in C++17 for some of the schemes in the LegoSNARK paper (plus others, e.g. multivariate polynomial commitments, algorithms for multilinear extensions, a product scheme from [eprint:2014/396](https://eprint.iacr.org/2014/396.pdf)).

**What this codebase is not:** it is not for production use; it is not extensively tested; it is not a full-fledged API or EDSL* for commit-and-prove SNARKs.

<sup>(*We are still considering an EDSL for commit-and-prove but moved  our focus from C++  to Rust as an implementation language as we found the latter to be a superior match)</sup>


## Overview

This repo includes commit-and-prove gadgets for the following relations:
- **matrix multiplication** (CPmmp in paper): <code>[src/examples/matrixsc.cc](src/examples/matrixsc.cc)</code>
- **generalized sumcheck** (CPsc in paper): <code>[src/gadgets/sumcheck.h](src/gadgets/sumcheck.h)</code>
- **"Linking" Pedersen commitments to vectors** in different bases, i.e. showing that they have the same opening (CPlink in paper): <code>[src/examples/cplink.cc](src/examples/cplink.cc)</code>
- **Hadamard product** (CPhad in paper): <code>[src/gadgets/hadamardsc.h](src/gadgets/hadamardsc.h)</code>
<!-- - **general arithmetic circuits** -->

It also includes code for:
- **multivariate polynomial commitments** (CPpoly in paper, partly based on an implementation of the scheme in [vSQL](https://web.eecs.umich.edu/~genkin/papers/vsql.pdf)  by Yupeng Zhang): <code>[src/gadgets/poly.h](src/gadgets/poly.h)</code>
- an **additional Hadamard product** based on the scheme in [Lipmaa's Commit-and-Prove paper](https://eprint.iacr.org/2014/396.pdf): <code>[src/gadgets/lipmaa.h](src/gadgets/lipmaa.h)</code>
- an **R1CS for matrix multiplication**: <code>[src/examples/legogrothmatrix.cc](src/examples/legogrothmatrix.cc)</code>



## Setup and Building Instructions

First, install the libraries and utilities required by libsnark and legosnark (see [here](https://github.com/scipr-lab/libsnark) for more detailed requirements). On several Ubuntu systems this can be done directly through the following command:
~~~~~~
sudo apt-get install build-essential cmake git libgmp3-dev libprocps-dev python-markdown libboost-all-dev libssl-dev
~~~~~~

 Clone the repo and set up submodules:
~~~~~~
 git clone https://github.com/imdea-software/legosnark.git
 cd legosnark
 git submodule update --init --recursive
~~~~~~

Build all dependencies:
~~~~~
mkdir -p build
cd build
cmake ..
cd depends
make -j8
sudo make -C libsnark install
~~~~~

To build library and executables:
~~~~~
cd ../src # Assuming you were in build/depends from the steps above
make -j8
~~~~~

To try an example, run e.g.:
~~~~~
examples/cplink
~~~~~

<!-- ### Using it as a library -->

## License

This code is licensed under either of the following licenses, at your discretion.

 * [Apache License Version 2.0](LICENSE-APACHE)
 * [MIT License](LICENSE-MIT)

Unless you explicitly state otherwise, any contribution that you submit to this library shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

[legosnark]: https://eprint.iacr.org/2019/142.pdf

## Reference paper

[LegoSNARK: Modular Design and Composition of Succinct Zero-Knowledge Proofs][legosnark]     
[Matteo Campanelli](https://www.github.com/matteocam), [Dario Fiore](https://github.com/dariofiore), [Ana&#239;s Querol](https://github.com/querolita)

CCS 2019

## Acknowledgements

This work has been supported by the Spanish Government under projects Datamantium (ref. RTC-2016-4930-7), SCUM (ref. RTI2018-102043-B-I00), and CRYPTOEPIC (refs. ERC2018-092822, EUR2019-103816), by the Madrid Regional Government under project BLOQUES (ref. S2018/TCS-4339) and by Protocol Labs. The project that gave rise to these results received the support of a fellowship from “la Caixa” Foundation (ID 100010434). The fellowship code is LCF/BQ/ES18/11670018.