#include "commit.h"

using namespace std;

CommOut operator+(const CommOut& a, const CommOut& b)  {
    if (a.lenXs != 1 || b.lenXs != 1) {
    	throw runtime_error("Operations CommOut are only for commitments of single elements.");
    }
    return CommOut(a.c+b.c, a.r+b.r, a.xs[0]+b.xs[0]);
}

CommOut operator-(const CommOut& a, const CommOut& b)  {
	if (a.lenXs != 1 || b.lenXs != 1) {
    	throw runtime_error("Operations CommOut are only for commitments of single elements.");
   }
   return CommOut(a.c-b.c, a.r-b.r, a.xs[0]-b.xs[0]);
}

