#ifndef GADGETS_MATRIX_H
#define GADGETS_MATRIX_H

/* Sparse Matrix types and functions */


#include <time.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>
#include <string>
#if !__has_include("optional")
	#include <experimental/optional>
	#include <experimental/functional>
	using namespace std::experimental;
#else
	#include <optional>
	#include <functional>
#endif

#include "globl.h"

#include "fmt/format.h"



using namespace std;
using namespace bn;



template<typename FieldT>
struct CoeffPos
{
	FieldT val;
	size_t pos;
	
	CoeffPos(const FieldT v, const size_t p) : val(v), pos(p) { }
};

using ColG1 = vector<CoeffPos<LG1>>;
using ColFr = vector<CoeffPos<LFr>>;


template<typename T>
void insertAsColMajor(const size_t r, const size_t c, const T &v, vector<vector<CoeffPos<T>>> &m)
{
	// NB: we assume that position (r,c) is not taken yet
	m[c].push_back(CoeffPos<T>(v, r));
}

template<typename T>
void insertRowAsColMajor(const size_t r, const size_t offset_c, const vector<T> &v, vector<vector<CoeffPos<T>>> &m)
{
	for (auto i = 0; i < v.size(); i++) {
		insertAsColMajor(r, offset_c+i, v[i], m);
	}
}

#endif
