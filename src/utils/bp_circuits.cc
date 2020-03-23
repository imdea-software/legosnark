#include "bp_circuits.h"
#include "globl.h"

#include <iostream>
#include <stdexcept>
using namespace std;

void parseNextSimpleLine(ifstream &inFile, size_t *x)
{
	string line;
	getline(inFile, line);
	// get substring after colon
	auto start = line.find(": ")+2;
	string fldAsStr = line.substr(start);
	//cout << fldAsStr << endl;
	if (x != nullptr) {
		*x = stoul(fldAsStr); // XXX:using conversion from unsigned long
	}
}

void appendMtxFld(BPCircuit &circ, const string mtxStr, size_t idx, size_t nr, size_t no)
{
	MYREQUIRE(nr <= circ.nrows);
	if (mtxStr == "WL") {
		circ.szWL[idx].nr = nr;
		circ.szWL[idx].no = no;
	} else if (mtxStr == "WR") {
		circ.szWR[idx].nr = nr;
		circ.szWR[idx].no = no;
	} else if (mtxStr == "WO") {
		circ.szWO[idx].nr = nr;
		circ.szWO[idx].no = no;
	}
}

BPCircuit parseNextBPCirc(ifstream &inFile)
{
	BPCircuit circ;
	string line;
	
	// NB: We check only the first line for something wrong
	if (!getline(inFile, line)) {
		throw runtime_error("Problem reading first line");
	}
	// skip "=CIRC= " string
	auto start = string("=CIRC= ").size();
	circ.name = line.substr(start);
	circ.name = line.substr(start);
	cout << circ.name << endl;
	
	// other fields
	parseNextSimpleLine(inFile, &circ.ncols); // n_gates (we basically transpose it)
	parseNextSimpleLine(inFile, nullptr); // n_commits (ignored)
	parseNextSimpleLine(inFile, &circ.nrows); // n_constraints
	parseNextSimpleLine(inFile, &circ.n_bits); // n_bits
	
	circ.szWL = vector<BPRowInfo>(circ.ncols);
	circ.szWR = vector<BPRowInfo>(circ.ncols);
	circ.szWO = vector<BPRowInfo>(circ.ncols);
	
	auto nextTokenPos = [](string &l, size_t bgn) {
		int end = l.size();
		for (auto i = bgn; i < end; i++) {
			if (l[i] == ':')
				return i;
			
		}
		throw runtime_error("Cannot find next token");
	};
	
	// main loop: till empty line or eof
	// bootstrapping loop
	getline(inFile, line);
	while(!inFile.eof() && line != "") {
		//cout << line << endl;
		
		// line format is "Wx:i:nr:no" where:
		// x is a character; i, nr, no are numbers of unknown size
		
		// extract first two characters
		string mtxStr = line.substr(0,2);
		int sndColonIdx = nextTokenPos(line, 3);
		string idxStr = line.substr(3, sndColonIdx-3);
		int trdColonIdx = nextTokenPos(line, sndColonIdx+1);
		//cout << "Colon found at " << sndColonIdx << " and " << trdColonIdx << " in " << line << endl;
		string nrStr = line.substr(sndColonIdx+1, trdColonIdx-sndColonIdx-1);
		string noStr = line.substr(trdColonIdx+1);
		//cout << mtxStr << ":" << idxStr << ":" << nrStr << ":" << noStr << endl;
		
		appendMtxFld(
		  circ, mtxStr, stoul(idxStr), stoul(nrStr), stoul(noStr));
		
		// get next line for next iteration
		getline(inFile, line);
	}
	
	return circ;
}

vector<BPCircuit> BPCircuit::readFromFile(std::string fn)
{
	vector<BPCircuit> ret;
	ifstream inFile(fn);
	if (!inFile) {
		throw runtime_error("Can't open file " + fn);
	}
	
	// NB: We assume there is no empty line at the end of the file
	while (!inFile.eof()) {
		ret.push_back(parseNextBPCirc(inFile));
	}
	
	return ret;
	
}
