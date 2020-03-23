#include "sparsemexp.h"

LG1 simplesparsemexp(const vector<LG1> &bases, const vector<CoeffPos<LFr>> &cps)
{
	#ifdef MULTICORE
    const size_t chunks = omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call omp_set_num_threads()
	#else
    const size_t chunks = 1;
	#endif
	
	return sparsemexpS<LG1, LFr, multi_exp_method_BDLO12>(bases, cps, chunks);
}

LG1 simplesparsemexp(const vector<CoeffPos<LG1>> &cps, const vector<LFr> &exps)
{
	#ifdef MULTICORE
    const size_t chunks = omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call omp_set_num_threads()
	#else
    const size_t chunks = 1;
	#endif
	
	return sparsemexpG<LG1, LFr, multi_exp_method_BDLO12>(cps, exps, chunks);
	
}

LFr sparseinnerproduct(const vector<LFr> &x, const vector<CoeffPos<LFr>> &cps)
{
	LFr res = LFr::zero();
	for (auto cp : cps) {
		res += cp.val*x[cp.pos];
	}
	return res;
}
