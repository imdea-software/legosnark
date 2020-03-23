#ifndef DBG_UTIL
#define DBG_UTIL

#include <vector>
#include <iostream>
#include <string>
#include "globl.h"
#include "matrix.h"

#define NO_CPDBG

namespace cpdbg {

  #ifdef NO_CPDBG
    template <typename T>
    void print(T, const std::string &, const std::string &) {

    }
  #endif

  #ifndef NO_CPDBG
    void printheader(const std::string &name, const std::string &ctx);

    template <typename T>
    void print(const std::vector<T> &v, const std::string &name, const std::string &ctx = "")
    {
      printheader(name, ctx);

      int cnt = 0;
      std::cout << "\t[ ";
      for (auto x : v) {
        std::cout << x;
        std::cout << (cnt < v.size()-1 ? "," : "");
        cnt++;
      }
      std::cout << "] " << std::endl;

    }

    //template<typename T>
    void print(const Ec1 x, const std::string &name, const std::string &ctx = "")
    {
      printheader(name, ctx);
      std::cout << "\t" << x << std::endl;
    }

    template<typename T>
    void print(const MyMat<T> &x, const std::string &name, const std::string &ctx = "")
    {
      printheader(name, ctx);
      #ifdef USE_SPARSE
        auto xDense = sparse2dense(x); // we convert to ensure order
        auto cols = xDense.cols();
        forEachItem<T>(xDense, [cols](int r, int c, T x) {cout << x << (c < cols-1 ? " " :"\n");} );
        //std::cout << "\t" << xDense << std::endl;

      #endif

      #ifdef USE_DENSE
        std::cout << "\t" << x << std::endl;
      #endif
    }


    void print(const Scalar x, const std::string &name, const std::string &ctx = "");

    template<typename T>
    void print(const T x, const std::string &name, const std::string &ctx  )
    {
      printheader(name, ctx);
      std::cout << "\t" << x << std::endl;
    }
  #endif

} // end of namespace cpdbg




#endif
