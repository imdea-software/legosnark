#ifndef SPARSE_MEXP
#define SPARSE_MEXP

/* Definitions for sparse matrix exponentiations */

#include "globl.h"
#include "matrix.h"

#include <libff/algebra/scalar_multiplication/multiexp.hpp>
using namespace libff;

template<typename T, typename FieldT, multi_exp_method Method>
T sparsemexpS(const vector<T> &bases, const vector<CoeffPos<FieldT>> &cps, const size_t chunks)
{

    const FieldT zero = FieldT::zero();
    const FieldT one = FieldT::one();
    
    std::vector<FieldT> p;
    std::vector<T> g;

    T acc = T::zero();

    size_t num_skip = 0;
    size_t num_add = 0;
    size_t num_other = 0;

    for (auto cp : cps)
    {
				auto curBase = bases[cp.pos];
        if (cp.val == zero)
        {
            // do nothing
            ++num_skip;
        }
        else if (cp.val == one)
        {
          
#ifdef USE_MIXED_ADDITION
            acc = acc.mixed_add(curBase);
#else

            acc = acc + (curBase);
#endif
            ++num_add;
        }
        else
        {
            p.push_back(cp.val);
            g.push_back(curBase);
            ++num_other;
        }
    }
    print_indent(); printf("* Elements of w skipped: %zu (%0.2f%%)\n", num_skip, 100.*num_skip/(num_skip+num_add+num_other));
    print_indent(); printf("* Elements of w processed with special addition: %zu (%0.2f%%)\n", num_add, 100.*num_add/(num_skip+num_add+num_other));
    print_indent(); printf("* Elements of w remaining: %zu (%0.2f%%)\n", num_other, 100.*num_other/(num_skip+num_add+num_other));

    return acc + multi_exp<T, FieldT, Method>(begin(g), end(g), begin(p), end(p), chunks);
}


template<typename T, typename FieldT, multi_exp_method Method>
T sparsemexpG(const vector<CoeffPos<T>> &cps, const vector<FieldT> &exps, const size_t chunks)
{

    const T zero = T::zero();
    const T one = T::one();
    
    std::vector<FieldT> p;
    std::vector<T> g;

    FieldT acc = FieldT::zero();

    for (auto cp : cps)
    {
      auto curExp = exps[cp.pos];
    
      if (cp.val == zero) {
        // do nothing
      } else if (cp.val == one) {
        // increase acc
        acc += curExp;
      } else {
        g.push_back(cp.val);
        p.push_back(curExp);
      }
    }
  
    return acc*one + multi_exp<T, FieldT, Method>(begin(g), end(g), begin(p), end(p), chunks);
}

LG1 simplesparsemexp(const vector<LG1> &bases, const vector<CoeffPos<LFr>> &cps);
LG1 simplesparsemexp(const vector<CoeffPos<LG1>> &cps, const vector<LFr> &exps);

LFr sparseinnerproduct(const vector<LFr> &x, const vector<CoeffPos<LFr>> &cps);

#endif
