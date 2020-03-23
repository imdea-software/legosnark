#include "benchmark.h"

void fmt_time(string msgtag, double t)
{
  fmt::print("{}: {} micros ({} s)\n", msgtag, t, t/1000000);
}

void print_bm(string msgtag, string bmlbl, const Benchmarkable &obj)
{
  double t;

  t =  obj.getTimingInMicrosFor(bmlbl);

  fmt_time(msgtag, t);
}

void print_sum_bm(string msgtag, string bmlbl1, string bmlbl2, const Benchmarkable &obj)
{
  double t1, t2;

  t1 =  obj.getTimingInMicrosFor(bmlbl1);
  t2 =  obj.getTimingInMicrosFor(bmlbl2);

  fmt_time(msgtag, t1+t2);
}