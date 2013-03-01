
#include<iostream>
#include<string>

#include<TFile.h>
#include<TH2D.h>
#include<THnSparse.h>
#include<TLorentzVector.h>
#include<TTree.h>
#include<TTreeIndex.h>

#include<basic_calcs.h>
#include<constants.h>
#include<initializer.h>

const double MIN_ENTRIES = 10.;

struct fiveDimBin {

	fiveDimBin() : _a(5,0), _b(5,0) { }

	fiveDimBin(double a0,
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

	fiveDimBin(double* a, double*b)
		: _a(5, 0.),
		  _b(5, 0.)
	{
		for(unsigned int i = 0; i < 5; ++i) {
			_a[i] = a[i] < b[i] ? a[i] : b[i];
			_b[i] = a[i] < b[i] ? b[i] : a[i];
		}
	}

	fiveDimBin(const std::vector<double>& a, const std::vector<double>& b)
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

	bool inBin(const std::vector<double>& x) {
		assert(x.size() == 5);
		for(unsigned int i = 0; i < 5; ++i) {
			if(x[i] < _a[i] or x[i] >= _b[i]) {
				return false;
			}
		}
		return true;
	}

	std::vector<double> _a;
	std::vector<double> _b;

	std::ostream& print(std::ostream& out) {
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

	bool operator<(const fiveDimBin& right) const {
		return _a[0] < right._a[0];
	}

};

void getAdaptiveBins(std::list<std::pair<TTree*, fiveDimBin> >& bins,
                     fiveDimBin bin,
                     TTree* inputTree,
                     unsigned int dim = 0,
                     bool debug = false)
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
		bins.push_back(std::pair<TTree*, fiveDimBin>(inputTree, bin));
	} else {
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
		}
		index = (TTreeIndex*)inputTree->GetTreeIndex();
		unsigned int half = (unsigned int)((entries / 2) + 0.5);
		inputTree->GetEntry(index->GetIndex()[half]);
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
		fiveDimBin newBin1(bin._a, upper1);
		fiveDimBin newBin2(lower2, bin._b);
		if(debug) {
			std::cout<<"bin1: "<<std::endl;
			newBin1.print(std::cout);
			std::cout<<"bin2: "<<std::endl;
			newBin2.print(std::cout);
		}
		TTree* inputTree1 = new TTree();
		TTree* inputTree2 = new TTree();
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
		if(not first) {
			delete inputTree;
		}
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
		getAdaptiveBins(bins, newBin1, inputTree1, dim, debug);
		getAdaptiveBins(bins, newBin2, inputTree2, dim, debug);
	}
	first = false;
}

void fillFiveDimHist(std::string inFileName) {

	std::string configfilename = "../../antok/config/monte_carlo_run2.yaml";
	antok::Initializer* initializer = antok::Initializer::instance();
	if(not initializer->readConfigFile(configfilename)) {
		std::cerr<<"Could not open config file. Aborting..."<<std::endl;
		exit(1);
	}
	const double& PION_MASS = antok::Constants::charged_pion_mass();

	TFile* inFile = TFile::Open(inFileName.c_str(), "READ");
	TTree* inTree = (TTree*)inFile->Get("Standard Event Selection/USR55");
	TFile* outFile = TFile::Open("bla.root", "RECREATE");
	outFile->cd();

//	Int_t bins[5] = {50, 50, 400, 400, 50};
//	Double_t xMin[5] = {-2, -2, -0.8, -0.8, 187};
//	Double_t xMax[5] = {2, 2, 0.8, 0.8, 196};

//	THnSparseD* hist = new THnSparseD("beamHist", "beamHist", 5, bins, xMin, xMax);
	TH1D* vtxX = new TH1D("vtxX", "vtxX", 10000, -2, 2);
	TH1D* vtxY = new TH1D("vtxY", "vtxY", 10000, -2, 2);
	TH1D* momX = new TH1D("momX", "momX", 50000, -0.8, 0.8);
	TH1D* momY = new TH1D("momY", "momY", 50000, -0.8, 0.8);
	TH1D* momZ = new TH1D("momZ", "momZ", 10000, 187, 196);

//	TH1D* edgeLength = new TH1D("edgeLength", "edgeLength", 501, -0.5, 500.5);

	double px1, py1, pz1;
	double px2, py2, pz2;
	double px3, py3, pz3;
	double px4, py4, pz4;
	double px5, py5, pz5;
	double gradx, grady;
	double primVx, primVy;

	inTree->SetBranchAddress("Mom_x1", &px1);
	inTree->SetBranchAddress("Mom_y1", &py1);
	inTree->SetBranchAddress("Mom_z1", &pz1);
	inTree->SetBranchAddress("Mom_x2", &px2);
	inTree->SetBranchAddress("Mom_y2", &py2);
	inTree->SetBranchAddress("Mom_z2", &pz2);
	inTree->SetBranchAddress("Mom_x3", &px3);
	inTree->SetBranchAddress("Mom_y3", &py3);
	inTree->SetBranchAddress("Mom_z3", &pz3);
	inTree->SetBranchAddress("Mom_x4", &px4);
	inTree->SetBranchAddress("Mom_y4", &py4);
	inTree->SetBranchAddress("Mom_z4", &pz4);
	inTree->SetBranchAddress("Mom_x5", &px5);
	inTree->SetBranchAddress("Mom_y5", &py5);
	inTree->SetBranchAddress("Mom_z5", &pz5);
	inTree->SetBranchAddress("gradx", &gradx);
	inTree->SetBranchAddress("grady", &grady);
	inTree->SetBranchAddress("X_primV", &primVx);
	inTree->SetBranchAddress("Y_primV", &primVy);

	double bx, by, bz;
	TTree* tempTree = new TTree("beamTreeTemp", "beamTreeTemp");
	tempTree->Branch("vertex_x_position", &primVx, "vertex_x_position/D");
	tempTree->Branch("vertex_y_position", &primVy, "vertex_y_position/D");
	tempTree->Branch("beam_momentum_x", &bx,"beam_momentum_x/D");
	tempTree->Branch("beam_momentum_y", &by,"beam_momentum_y/D");
	tempTree->Branch("beam_momentum_z", &bz,"beam_momentum_z/D");
	TTree* outTree = (TTree*)tempTree->Clone("beamTree");

	std::vector<TLorentzVector> particles;
	particles.resize(6);
	unsigned int entries = inTree->GetEntries();

	for(unsigned int i = 0; i < entries; ++i) {

		inTree->GetEntry(i);

		particles[1].SetXYZM(px1, py1, pz1, PION_MASS);
		particles[2].SetXYZM(px2, py2, pz2, PION_MASS);
		particles[3].SetXYZM(px3, py3, pz3, PION_MASS);
		particles[4].SetXYZM(px4, py4, pz4, PION_MASS);
		particles[5].SetXYZM(px5, py5, pz5, PION_MASS);

		TLorentzVector pSum;
		for(unsigned int j = 1; j < 6; ++j) {
			pSum += particles.at(j);
		}
		particles[0] = antok::get_beam_energy(TVector3(gradx, grady, 1.), pSum);
		TVector3 beam = particles[0].Vect();
		bx = beam.X();
		by = beam.Y();
		bz = beam.Z();

//		double values[5] = {primVx, primVy, bx, by, bz};

		vtxX->Fill(primVx);
		vtxY->Fill(primVy);
		momX->Fill(bx);
		momY->Fill(by);
		momZ->Fill(bz);

//		hist->Fill(values);
		tempTree->Fill();

		if(i % 100000 == 0) {
			std::cout<<"Entry "<<i<<" of "<<entries<<std::endl;
		}

	}
/*
	hist->Projection(1, 0)->Write();
	hist->Projection(3, 2)->Write();
	hist->Projection(2, 1)->Write();
	hist->Projection(2, 0)->Write();
	hist->Projection(3, 1)->Write();
	hist->Projection(3, 0)->Write();
	hist->Projection(4, 0)->Write();
	hist->Projection(4, 1)->Write();
	hist->Projection(4, 2)->Write();
	hist->Projection(4, 3)->Write();

	std::cout<<"Histogram has "<<hist->GetNbins()<<" filled bins."<<std::endl;
*/

	std::vector<double> lowerCorner(5, 0);
	std::vector<double> upperCorner(5, 0);
	bool first = true;
	for(unsigned int i = 0; i < entries; ++i) {
		tempTree->GetEntry(i);
		if(first) {
			lowerCorner[0] = primVx;
			lowerCorner[1] = primVy;
			lowerCorner[2] = bx;
			lowerCorner[3] = by;
			lowerCorner[4] = bz;
			upperCorner[0] = primVx;
			upperCorner[1] = primVy;
			upperCorner[2] = bx;
			upperCorner[3] = by;
			upperCorner[4] = bz;
			first = false;
		} else {
			if(primVx < lowerCorner[0]) lowerCorner[0] = primVx;
			if(primVy < lowerCorner[1]) lowerCorner[1] = primVy;
			if(bx < lowerCorner[2]) lowerCorner[2] = bx;
			if(by < lowerCorner[3]) lowerCorner[3] = by;
			if(bz < lowerCorner[4]) lowerCorner[4] = bz;
			if(primVx > upperCorner[0]) upperCorner[0] = primVx;
			if(primVy > upperCorner[1]) upperCorner[1] = primVy;
			if(bx > upperCorner[2]) upperCorner[2] = bx;
			if(by > upperCorner[3]) upperCorner[3] = by;
			if(bz > upperCorner[4]) upperCorner[4] = bz;
		}
	}
	fiveDimBin bin(lowerCorner, upperCorner);
	std::cout<<"Got first bin: "<<std::endl;
	bin.print(std::cout);

	std::list<std::pair<TTree*, fiveDimBin> > adaptiveBins;
	getAdaptiveBins(adaptiveBins, bin, tempTree);
	unsigned int nBins = adaptiveBins.size();
	std::cout<<"Split phase space in "<<nBins<<" bins."<<std::endl;
/*
	double primVx_sigma, primVy_sigma, bx_sigma, by_sigma, bz_sigma;
	int binContent = 0;
	TTree* sigmaTree = new TTree("sigmaTree", "sigmaTree");
	sigmaTree->Branch("vertex_x_position_sigma", &primVx_sigma, "vertex_x_position_sigma/D");
	sigmaTree->Branch("vertex_y_position_sigma", &primVy_sigma, "vertex_y_position_sigma/D");
	sigmaTree->Branch("beam_momentum_x_sigma", &bx_sigma,"beam_momentum_x_sigma/D");
	sigmaTree->Branch("beam_momentum_y_sigma", &by_sigma,"beam_momentum_y_sigma/D");
	sigmaTree->Branch("beam_momentum_z_sigma", &bz_sigma,"beam_momentum_z_sigma/D");
	sigmaTree->Branch("bin_content", &binContent, "bin_content/I");
*/
	double primVx_sigma, primVy_sigma, bx_sigma, by_sigma, bz_sigma;
	int binContent = 0;
	outTree->SetBranchAddress("vertex_x_position", &primVx);
	outTree->SetBranchAddress("vertex_y_position", &primVy);
	outTree->SetBranchAddress("beam_momentum_x", &bx);
	outTree->SetBranchAddress("beam_momentum_y", &by);
	outTree->SetBranchAddress("beam_momentum_z", &bz);
	outTree->Branch("vertex_x_position_sigma", &primVx_sigma, "vertex_x_position_sigma/D");
	outTree->Branch("vertex_y_position_sigma", &primVy_sigma, "vertex_y_position_sigma/D");
	outTree->Branch("beam_momentum_x_sigma", &bx_sigma,"beam_momentum_x_sigma/D");
	outTree->Branch("beam_momentum_y_sigma", &by_sigma,"beam_momentum_y_sigma/D");
	outTree->Branch("beam_momentum_z_sigma", &bz_sigma,"beam_momentum_z_sigma/D");
	outTree->Branch("bin_content", &binContent, "bin_content/I");

	std::cout << "Calculating and saving sigmas." << std:: endl;

	unsigned binNumber = 0;
	for(std::list<std::pair<TTree*, fiveDimBin> >::const_iterator binIt = adaptiveBins.begin(); binIt != adaptiveBins.end(); ++binIt) {
		std::cout<<"Bin "<<++binNumber<<" of "<<nBins<<std::endl;
		const fiveDimBin& currentBin = binIt->second;
		TTree* currentTree = binIt->first;
		currentTree->SetBranchAddress("vertex_x_position", &primVx);
		currentTree->SetBranchAddress("vertex_y_position", &primVy);
		currentTree->SetBranchAddress("beam_momentum_x", &bx);
		currentTree->SetBranchAddress("beam_momentum_y", &by);
		currentTree->SetBranchAddress("beam_momentum_z", &bz);
		binContent = currentTree->GetEntries();
		primVx_sigma = (currentBin._b[0] - currentBin._a[0]) / entries;
		primVy_sigma = (currentBin._b[1] - currentBin._a[1]) / entries;
		bx_sigma = (currentBin._b[2] - currentBin._a[2]) / entries;
		by_sigma = (currentBin._b[3] - currentBin._a[3]) / entries;
		bz_sigma = (currentBin._b[4] - currentBin._a[4]) / entries;
		for(int i = 0; i < binContent; ++i) {
			currentTree->GetEntry(i);
			outTree->Fill();
		}
	}
	tempTree->Delete("all");

/*	for(unsigned int i = 0; i < entries; ++i) {
		tempTree->GetEntry(i);
		std::vector<double> x(5, 0);
		x[0] = primVx;
		x[1] = primVx;
		x[2] = bx;
		x[3] = by;
		x[4] = bz;
		bool found = false;
		for(std::map<fiveDimBin, TTree*>::const_iterator binIt = adaptiveBins.begin(); binIt != adaptiveBins.end(); ++binIt) {
			fiveDimBin currentBin = binIt->first;
			if(currentBin.inBin(x)) {
				binContent = binIt->second->GetEntries();
				primVx_sigma = (currentBin._b[0] - currentBin._a[0]) / entries;
				primVy_sigma = (currentBin._b[1] - currentBin._a[1]) / entries;
				bx_sigma = (currentBin._b[2] - currentBin._a[2]) / entries;
				by_sigma = (currentBin._b[3] - currentBin._a[3]) / entries;
				bz_sigma = (currentBin._b[4] - currentBin._a[4]) / entries;
				sigmaTree->Fill();
				found = true;
				break;
			}
		}
		assert(found);
		if(i % 100000 == 0) {
			std::cout<<"Entry "<<i<<" of "<<entries<<std::endl;
		}
		std::cout<<i<<std::endl;
	}
	tempTree->AddFriend(sigmaTree);
*/
/*
	double primVx_sigma, primVy_sigma, bx_sigma, by_sigma, bz_sigma;
	int cubeSize = 1;
	double binContent = 0;
	TTree* sigmaTree = new TTree("sigmaTree", "sigmaTree");
	sigmaTree->Branch("vertex_x_position_sigma", &primVx_sigma, "vertex_x_position_sigma/D");
	sigmaTree->Branch("vertex_y_position_sigma", &primVy_sigma, "vertex_y_position_sigma/D");
	sigmaTree->Branch("beam_momentum_x_sigma", &bx_sigma,"beam_momentum_x_sigma/D");
	sigmaTree->Branch("beam_momentum_y_sigma", &by_sigma,"beam_momentum_y_sigma/D");
	sigmaTree->Branch("beam_momentum_z_sigma", &bz_sigma,"beam_momentum_z_sigma/D");
	sigmaTree->Branch("cube_edge_length", &cubeSize, "cube_edge_length/I");
	sigmaTree->Branch("bin_content", &binContent, "bin_content/D");

	const double MIN_ENTRIES = 10.;

	for(unsigned int i = 0; i < entries; ++i) {

		outTree->GetEntry(i);
		double values[5] = {primVx, primVy, bx, by, bz};
		int bin = hist->GetBin(values, false);
		int location[5] = {0, 0, 0, 0, 0};
		binContent = hist->GetBinContent(bin, location);
		cubeSize = 1;
		if(binContent < MIN_ENTRIES) {
			std::pair<int, int> verXRange(location[0], location[0]);
			std::pair<int, int> verYRange(location[1], location[1]);
			std::pair<int, int> momXRange(location[2], location[2]);
			std::pair<int, int> momYRange(location[3], location[3]);
			std::pair<int, int> momZRange(location[4], location[4]);

			while(binContent < MIN_ENTRIES) {

				verXRange.first -= 1; verXRange.second += 1;
				verYRange.first -= 1; verYRange.second += 1;
				momXRange.first -= 1; momXRange.second += 1;
				momYRange.first -= 1; momYRange.second += 1;
				momZRange.first -= 1; momZRange.second += 1;
				int step = verXRange.second - verXRange.first;
				cubeSize = step + 1;

				for(int dim = 0; dim < 5; ++dim) {
					for(int l = (dim > 0) ? verXRange.first + 1 : verXRange.first; l <= ((dim > 0) ? verXRange.second - 1 : verXRange.second); l += (dim == 0) ? step : 1) {
						for(int m = (dim > 1) ? verYRange.first + 1 : verYRange.first; m <= ((dim > 1) ? verYRange.second - 1 : verYRange.second); m += (dim == 1) ? step : 1) {
							for(int n = (dim > 2) ? momXRange.first + 1 : momXRange.first; n <= ((dim > 2) ? momXRange.second - 1 : momXRange.second); n += (dim == 2) ? step : 1) {
								for(int o = (dim > 3) ? momYRange.first + 1 : momYRange.first; o <= ((dim > 3) ? momYRange.second - 1 : momYRange.second); o += (dim == 3) ? step : 1) {
									for(int p = momZRange.first; p <= momZRange.second; p += (dim == 4) ? step : 1) {
										int coord[5] = {l, m, n, o, p};
										int binnumber = hist->GetBin(coord, false);
										if(binnumber > 0) {
											binContent += hist->GetBinContent(bin);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		edgeLength->Fill(cubeSize);
		primVx_sigma = (cubeSize * (hist->GetAxis(0)->GetBinWidth(location[0]))) / binContent;
		primVy_sigma = (cubeSize * (hist->GetAxis(1)->GetBinWidth(location[1]))) / binContent;
		bx_sigma = (cubeSize * (hist->GetAxis(2)->GetBinWidth(location[2]))) / binContent;
		by_sigma = (cubeSize * (hist->GetAxis(3)->GetBinWidth(location[3]))) / binContent;
		bz_sigma = (cubeSize * (hist->GetAxis(4)->GetBinWidth(location[4]))) / binContent;
		sigmaTree->Fill();

		if(not (i % 100)) {
			std::cout<<i<<std::endl;
		}

	}

	outTree->AddFriend(sigmaTree);
*/
//	hist->Write();
	outFile->Write();
	outFile->Close();

}

int main(int argc, char* argv[]) {
	if(argc == 2) {
		fillFiveDimHist(argv[1]);
	} else {
		std::cerr<<"Wrong number of arguments, is "<<argc<<", should be 1."<<std::endl;
	}
}
