#ifndef CP_UTIL_H
#define CP_UTIL_H

#include <vector>
#include <algorithm>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <stdexcept>

#include "globl.h"


/* General utility classes/functions */

namespace cputil {
  using namespace std;

  template<typename T, typename U>
  vector<U> map(const vector <T> &src, function<U(T)> fn)
  {
    vector<U> res;
    transform(src.begin(), src.end(), back_inserter(res), fn);
    return res;
  }

  // XXX: For now only three vectors
  template<typename T>
  vector<T> flatten(const vector<T> &a, const vector<T> &b, const vector<T> &c)
  {
    vector<T> dst;

    auto concatFn = [&dst](const vector<T> &v) {
        dst.insert( dst.end(), v.begin(), v.end() ); };

    concatFn(a);
    concatFn(b);
    concatFn(c);

    return dst;
  }



  template<typename T>
  vector<T> flatten(const vector<T> &a, const vector<T> &b)
  {
    vector<T> emptyVec;
    return flatten(a, b, emptyVec);
  }

  template<typename T>
  void dumpIntoFile(const std::vector<T> &v, ofstream &outFile)
  {
    outFile << v.size() << "\n";
    auto outFn = [&outFile](const T& x) { outFile << x << "\n"; };
    std::for_each(std::begin(v), std::end(v), outFn);
  }

  template<typename T>
  void dumpIntoFile(const std::vector<T> &v, std::string fn)
  {
    std::ofstream outFile(fn);
    if (outFile.fail()) {
      // We failed to open the file: throw an exception.
      throw runtime_error("Failed opening " + fn);
    }
    dumpIntoFile(v, outFile);
  }

  template<typename T>
  void loadFromFile(std::vector<T> &v, ifstream &inFile)
  {
    std::size_t sz;
    inFile >> sz;
    v.resize(sz);

    for (auto i = 0; i < sz; i++) {
      inFile >> v[i];
    }
  }

  template<typename T>
  void loadFromFile(std::vector<T> &v, std::string fn)
  {
    std::ifstream inFile(fn);
    if (inFile.fail()) {
      // We failed to open the file: throw an exception.
      throw runtime_error("Failed opening " + fn);
    }
    loadFromFile(v, inFile);
  }

  inline long log2ceiled(long x)
  {
    auto l = log2(x);
    return ceil(l);
  }
  
  template<typename T>
  vector<T> concat3(const vector<T> &a, const vector<T> &b, const vector<T> &c)
  {
	vector<T> abc;
	auto append2abc = [&abc](auto &v) { 
		abc.insert(abc.end(), v.begin(), v.end());
	};
	append2abc(a);
	append2abc(b);
	append2abc(c);
	return abc;
  }



template<typename T, typename FieldT>
vector<T> simpleBatchExp(const T &base, const vector<FieldT> &exps)
{
	long n;
	size_t g_exp_count;
	size_t g_window;
	window_table<T> g_table;
	// setup multiexp
	const size_t fldBitSz = FieldT::size_in_bits();
	n = exps.size();
	g_exp_count = n;
	g_window = get_exp_window_size<T>(g_exp_count);
	g_table = get_window_table(fldBitSz, g_window, base);
	
	return batch_exp(fldBitSz, g_window, g_table, exps);
}



inline void populate_from_file_dist(vector<LFr> &x, size_t n, string fn)
{
	vector<LFr> fileDump;
	//LFr tmp;
	
	ifstream f(fn);

	string line;
	getline(f, line);
	while(!f.eof()) {
	  LFr tmp(line.c_str());

	  fileDump.push_back(tmp);
	  getline(f, line);

	}
	
	
	
	// actually populate here
	int idx = 0;
	for (auto j = 0; j < n; j++) {
	  x.push_back(fileDump[idx]);
	  idx = (idx + 1) % fileDump.size();
	}
	
}



} // end namespace cputil


#endif
