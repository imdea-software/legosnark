#include "dbgutil.h"

namespace cpdbg {
  void printheader(const std::string &name, const std::string &ctx)
  {
    std::cout << std::endl;
    if (ctx != "") {
      std::cout << "-- " << ctx << " --" << std::endl;
    }
    std::cout << "<Value of " << name << ">" << std::endl;
  }


#ifndef BN_SUPPORT_SNARK
  void print(const Scalar x, const std::string &name, const std::string &ctx  )
  {
    printheader(name, ctx);
    std::cout << "\t" << x << std::endl;
  }
 #endif

}
