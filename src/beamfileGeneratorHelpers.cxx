#include<beamfileGeneratorHelpers.h>

#include<cmath>

#include<TTree.h>
#include<TTreeIndex.h>

antok::beamfileGenerator::fiveDimBin::fiveDimBin(double a0,
                                                 double a1,
                                                 double a2,
                                                 double a3,
                                                 double a4,
                                                 double b0,
                                                 double b1,
                                                 double b2,
                                                 double b3,
                                                 double b4)
	: _a(5, 0.),
	  _b(5, 0.)
{
	_a[0] = a0 < b0 ? a0 : b0;
	_a[1] = a1 < b1 ? a1 : b1;
	_a[2] = a2 < b2 ? a2 : b2;
	_a[3] = a3 < b3 ? a3 : b3;
	_a[4] = a4 < b4 ? a4 : b4;
	_b[0] = a0 < b0 ? b0 : a0;
	_b[1] = a1 < b1 ? b1 : a1;
	_b[2] = a2 < b2 ? b2 : a2;
	_b[3] = a3 < b3 ? b3 : a3;
	_b[4] = a4 < b4 ? b4 : a4;
}

antok::beamfileGenerator::fiveDimBin::fiveDimBin(const std::vector<double>& a, const std::vector<double>& b)
	: _a(5, 0.),
	  _b(5, 0.)
{
	assert(a.size() == 5);
	assert(b.size() == 5);
	for(unsigned int i = 0; i < 5; ++i) {
		_a[i] = a[i] < b[i] ? a[i] : b[i];
		_b[i] = a[i] < b[i] ? b[i] : a[i];
	}
}

void antok::beamfileGenerator::fiveDimBin::set(const std::vector<double>& a, const std::vector<double>& b)
{
	assert(a.size() == 5);
	assert(b.size() == 5);
	for(unsigned int i = 0; i < 5; ++i) {
		_a[i] = a[i] < b[i] ? a[i] : b[i];
		_b[i] = a[i] < b[i] ? b[i] : a[i];
	}
}

bool antok::beamfileGenerator::fiveDimBin::inBin(const std::vector<double>& x) const
{
	assert(x.size() == 5);
	for(unsigned int i = 0; i < 5; ++i) {
		if(x[i] < _a[i] or x[i] >= _b[i]) {
			return false;
		}
	}
	return true;
}

std::ostream& antok::beamfileGenerator::fiveDimBin::print(std::ostream& out) const {
	out << "Five dimensional bin: " << std::endl;
	out << "    lower Corner ....... [" << _a[0];
	for(unsigned int i = 1; i < 5; ++i) {
		out << ", " << _a[i];
	}
	out << "]" << std::endl;
	out << "    upper Corner ....... [" << _b[0];
	for(unsigned int i = 1; i < 5; ++i) {
		out << ", " << _b[i];
	}
	out << "]" << std::endl;
	return out;
}

void antok::beamfileGenerator::getAdaptiveBins(std::list<std::pair<TTree*, antok::beamfileGenerator::fiveDimBin> >& bins,
                                               const antok::beamfileGenerator::fiveDimBin& bin,
                                               TTree* inputTree,
                                               unsigned int dim,
                                               bool debug)
{
	static bool first = true;
	static std::vector<std::string> BRANCH_NAMES;
	if(first) {
		BRANCH_NAMES.resize(5);
		BRANCH_NAMES[0] = "vertex_x_position";
		BRANCH_NAMES[1] = "vertex_y_position";
		BRANCH_NAMES[2] = "beam_momentum_x";
		BRANCH_NAMES[3] = "beam_momentum_y";
		BRANCH_NAMES[4] = "beam_momentum_z";
	}
	long entries = inputTree->GetEntries();
	if(debug) {
		std::cout<<"------------------------------------------------------"<<std::endl;
		std::cout<<"called with " << entries << " entries." << std::endl;
		std::cout<<"input bin:"<<std::endl;
		bin.print(std::cout);
	}
	if(entries < (2 * MIN_ENTRIES)) {
		bins.push_back(std::pair<TTree*, antok::beamfileGenerator::fiveDimBin>(inputTree, bin));
	} else {
		antok::beamfileGenerator::fiveDimBin newBin1;
		antok::beamfileGenerator::fiveDimBin newBin2;
		TTree* inputTree1 = 0;
		TTree* inputTree2 = 0;
		{
			std::vector<double> x(5, 0);
			for(unsigned int i = 0; i < 5; ++i) {
				inputTree->SetBranchAddress(BRANCH_NAMES[i].c_str(), &x[i]);
			}
			TTreeIndex* index = 0;
			switch(dim) {
				case 0:
					inputTree->BuildIndex("0", "vertex_x_position*1000000000000");
					break;
				case 1:
					inputTree->BuildIndex("0", "vertex_y_position*1000000000000");
					break;
				case 2:
					inputTree->BuildIndex("0", "beam_momentum_x*1000000000000");
					break;
				case 3:
					inputTree->BuildIndex("0", "beam_momentum_y*1000000000000");
					break;
				case 4:
					inputTree->BuildIndex("0", "beam_momentum_z*1000000000000");
					break;
				default:
					assert(false);
					break;
			}
			index = (TTreeIndex*)inputTree->GetTreeIndex();
			unsigned int half = (unsigned int)((entries / 2) + 0.5);
			inputTree->GetEntry(index->GetIndex()[half]);
			delete index;
			double middle = x[dim];
			if(debug) {
				std::cout<<"entries: " << entries << ", half: " << half << std::endl;
				std::cout<<"middle coordinates: ["<<x[0];
				for(unsigned int i = 1; i < 4; ++i) {
					std::cout<<", "<<x[i];
				}
				std::cout<<"]"<<std::endl;
				std::cout<<"half entry: " << index->GetIndex()[half] << std::endl;
				std::cout<<"dim: " << dim << ", middle: " << middle << std::endl;
			}
			std::vector<double> upper1 = bin._b;
			std::vector<double> lower2 = bin._a;
			upper1[dim] = middle;
			lower2[dim] = middle;
			newBin1.set(bin._a, upper1);
			newBin2.set(lower2, bin._b);
			if(debug) {
				std::cout<<"bin1: "<<std::endl;
				newBin1.print(std::cout);
				std::cout<<"bin2: "<<std::endl;
				newBin2.print(std::cout);
			}
			inputTree1 = new TTree();
			inputTree2 = new TTree();
			for(unsigned int i = 0; i < 5; ++i) {
				inputTree1->Branch(BRANCH_NAMES[i].c_str(), &x[i]);
				inputTree2->Branch(BRANCH_NAMES[i].c_str(), &x[i]);
			}
			for(unsigned int i = 0; i < entries; ++i) {
				inputTree->GetEntry(i);
				if(x[dim] < middle) {
					inputTree1->Fill();
				} else {
					inputTree2->Fill();
				}
			}
			inputTree1->SetEstimate(-1);
			inputTree2->SetEstimate(-1);
			inputTree1->SetCacheSize(0);
			inputTree2->SetCacheSize(0);
			inputTree->SetBit(kCanDelete);
			inputTree->SetBit(kMustCleanup);
			delete inputTree;
			if(dim == 4) {
				dim = 0;
			} else {
				++dim;
			}
			if(debug) {
				std::cout<<"inputTree1 has " << inputTree1->GetEntries() << " entries." << std::endl;
				std::cout<<"inputTree2 has " << inputTree2->GetEntries() << " entries." << std::endl;
			}
			assert(std::abs((inputTree1->GetEntries() - inputTree2->GetEntries())) <= 1);
		}
		antok::beamfileGenerator::getAdaptiveBins(bins, newBin1, inputTree1, dim, debug);
		antok::beamfileGenerator::getAdaptiveBins(bins, newBin2, inputTree2, dim, debug);
	}
	first = false;
}
