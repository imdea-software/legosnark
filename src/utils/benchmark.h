#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <memory>
#include <string>
#include <unordered_map>
#include <set>
#include <chrono>
#include <utility>
#include <iostream>
#include <string>
#include <vector>
#include <functional>

#include "util.h"
#include "fmt/format.h"

using namespace std;

class Benchmarkable;

void print_bm(string msgtag, string bmlbl, const Benchmarkable &obj);
void print_sum_bm(string msgtag, string bmlbl1, string bmlbl2, const Benchmarkable &obj);
void fmt_time(string msgtag, double t);



struct TimeDelta {
   chrono::high_resolution_clock::time_point begin, end;
   bool started = false;
   bool over = false;

   static unique_ptr<TimeDelta> fromNow()
   {
     auto td = make_unique<TimeDelta>();
     td->begin = chrono::high_resolution_clock::now();
  
     td->started = true;
     return move(td);
    }
    
    static double timeFunction(function<void()> fn)
    {
		auto tdPtr = fromNow();
		fn();
		tdPtr->stop();
		return tdPtr->getDelta();
	}
	
	static double runAndAverage(function<void()> fn, unsigned T)
	{
		double res = 0;
		for (auto i = 1; i <= T; i++) {
			res += timeFunction(fn);
		}
		double out = res/T;
		return out;
	}

    void stop()
    {
      end = chrono::high_resolution_clock::now ();
 
      over = true;
    }

    long getDelta() const
    {
      if (!startedAndOver()) {
        throw runtime_error("Not started or not over!");
      }
      return chrono::duration_cast < chrono::microseconds >
        (end - begin).count ();
    }

    bool startedAndOver() const
    {
      return started && over;
    }

};




class Benchmark {
protected:
  unordered_map<string, unique_ptr<TimeDelta>> idmap;

    void _registerAndStartTiming(string obj, string session)
    {
      auto obj_session = obj+session;
      idmap[obj_session] = TimeDelta::fromNow();
    }

    void _stopTiming(string obj, string session)
    {
      auto obj_session = obj+session;
      idmap[obj_session]->stop();
    }

    long _getTiming(string obj, string session)
    {
      auto obj_session = obj+session;
      // Throw if we cannot find obj_session
      if (idmap.find(obj_session) == idmap.end()) {
		  throw runtime_error(
			fmt::format("No session {} in object {}", session, obj));
	  }
      return idmap[obj_session]->getDelta();
    }

    void _copyTimeDelta(
      string obj1, string session1,
      string obj2, string session2, bool okToReplace)
    {
      auto obj_session1 = obj1+session1;
      auto obj_session2 = obj2+session2;

      assert(alreadyExists(obj1, session1));
      assert(okToReplace || !alreadyExists(obj2, session2));

      TimeDelta td = *idmap[obj_session1];
      idmap[obj_session2] = make_unique<TimeDelta>(td);
    }

public:

  void startTiming(string obj, string session, bool okToReplace)
  {
    // we forbid timing something we already timed
    if (!okToReplace && alreadyExists(obj, session)) {
      throw runtime_error(
        fmt::format("({},{}) already exists as benchmark.", obj, session));
    }
    _registerAndStartTiming(obj, session);
  }

  void stopTiming(string obj, string session)
  {
    // we forbid timing something we already timed
    if (!alreadyExists(obj, session)) {
      auto msg =
        fmt::format(
          "({},{}) does not exists as a benchmark and cannot be stopped.",
          obj,
          session);
      throw runtime_error(msg);
    }
    _stopTiming(obj, session);
  }

  void copyTimeDelta(
    string obj1, string session1,
    string obj2, string session2, bool okToReplace)
  {
    _copyTimeDelta(obj1, session1, obj2, session2, okToReplace);
  }

  long getTiming(string obj, string session)
  {
    return _getTiming(obj, session);
  }

  bool alreadyExists(string obj, string session)
  {
    auto obj_session = obj+session;
    return (idmap.count(obj_session) != 0);
  }
};


class Benchmarkable {
protected:
  shared_ptr<Benchmark> pBm;
  string objId; // objId identifies this object to the benchmarker
  // HACK: should not use raw pointers in next line
  vector<pair<Benchmarkable *, string>> slaves;

public:

  void runAndAverage(
	function<void()> bmFn,
	string sessionId, double &out,
	int T)
  {
    double res = 0;
    for (auto i = 1; i <= T; i++) {
      bmFn();
      res += getTimingInMicrosFor(sessionId);
    }
    out = res/T;
  }
  
  void runAndAverage2(
	function<void()> bmFn,
	string sessionId1, double &out1,
	string sessionId2, double &out2,
	int T)
  {
    double res1, res2;
    res1 = res2 = 0;
    for (auto i = 1; i <= T; i++) {
      bmFn();
      res1 += getTimingInMicrosFor(sessionId1);
	  res2 += getTimingInMicrosFor(sessionId2);
    }
    out1 = res1/T;
    out2 = res2/T;
  }
  
   void runAndAverage3(
	function<void()> bmFn,
	string sessionId1, double &out1,
	string sessionId2, double &out2,
	string sessionId3, double &out3,
	int T)
  {
    double res1, res2, res3;
    res1 = res2 = res3 = 0;
    for (auto i = 1; i <= T; i++) {
      bmFn();
      res1 += getTimingInMicrosFor(sessionId1);
	  res2 += getTimingInMicrosFor(sessionId2);
	  res3 += getTimingInMicrosFor(sessionId3);

    }
    out1 = res1/T;
    out2 = res2/T;
	out3 = res3/T;

  }

  void addBenchmarkSlave(Benchmarkable *s, string id)
  {
    slaves.push_back(make_pair(s, id));
  }

  virtual void setBenchmark(shared_ptr<Benchmark> _pBm, string _objId )
  {
    pBm = _pBm;
    objId = _objId;
#if defined(__GNUC__)
    for (auto slave : slaves) {
        auto s = slave.first;
        auto id = slave.second;
#else
    for (auto [s, id] : slaves) {
#endif
        s->setBenchmark(_pBm, id);
    }
  }

  // sessionId is "what we are measuring" (e.g. "verification")
  virtual void startBenchmark(string sessionId, bool okToReplace = true, bool silentSkipIfNoBm = true)
  {

    if(!pBm && !silentSkipIfNoBm) {
      throw runtime_error("Can't start benchmark if it's not set.");
    }
    if (pBm) {
      pBm->startTiming(objId, sessionId, okToReplace);
    }
  }

  virtual void stopBenchmark(string sessionId,  bool silentSkipIfNoBm = true)
  {
    if(!pBm && !silentSkipIfNoBm) {
      throw runtime_error("Can't stop benchmark if it's not set.");
    }
    if (pBm) {
      pBm->stopTiming(objId, sessionId);
    }
  }

  virtual long getTimingInMicrosFor(string sessionId) const
  {
    if(!pBm) {
      throw runtime_error("Can't get timing for benchmark if it's not set.");
    }
    return pBm->getTiming(objId, sessionId);
  }

  // works only if the two objects share the same pBm
  void applyBenchmarkFrom(
      const Benchmarkable &src,
      string srcSess,
      string dstSess,
      bool okToReplace = true,
      bool silentSkipIfNoBm = true)
  {
    string srcId = src.objId;
    if (!pBm && silentSkipIfNoBm) {
      return;
    } else if (!pBm && !silentSkipIfNoBm) {
      throw runtime_error("Copying benchmarks but no benchmark object found.");
    }
    pBm->copyTimeDelta(srcId, srcSess, objId, dstSess, okToReplace);
  }
};

#endif
