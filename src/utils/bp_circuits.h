#ifndef BP_CIRCUITS_H
#define BP_CIRCUITS_H

/* Abstractions for BulletProof style (or Bootle-style) circuits */

#include <string>
#include <cstdlib>
#include <memory>
#include <fstream>
#include <vector>

struct BPRowInfo {
	size_t nr, no;
	BPRowInfo() : nr(0), no(0) {}
};

struct BPCircuit {
	std::string name;
	size_t nrows; // should be n_gates (from the C implementation of BP)
	size_t ncols; // should be n_constraints
	size_t n_bits; // extra: implicit # of bit constraints
	std::vector<BPRowInfo> szWL, szWR, szWO; // number of non-emtpy entries in each row
	
	static std::vector<BPCircuit> readFromFile(std::string fn);
	
	friend BPCircuit parseNextBPCirc(std::ifstream &inFile);

protected:
	BPCircuit() {}
};

#endif
