
#include<iostream>
#include<string>

#include<TFile.h>
#include<TH1D.h>
#include<TLorentzVector.h>
#include<TTree.h>

#include<basic_calcs.h>
#include<beamfile_generator_helpers.h>
#include<constants.h>
#include<initializer.h>


void fillFiveDimHist(std::string inFileName, std::string outFileName) {

	TFile* inFile = TFile::Open(inFileName.c_str(), "READ");
	TTree* inTree = (TTree*)inFile->Get("testBeamTree");
	TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	outFile->cd();

	TH1D* vtxX = new TH1D("vtxX", "vtxX", 10000, -2, 2);
	TH1D* vtxY = new TH1D("vtxY", "vtxY", 10000, -2, 2);
	TH1D* momX = new TH1D("momX", "momX", 50000, -0.8, 0.8);
	TH1D* momY = new TH1D("momY", "momY", 50000, -0.8, 0.8);
	TH1D* momZ = new TH1D("momZ", "momZ", 10000, 187, 196);


	double primVx, primVy;
	double bx, by, bz;

	inTree->SetBranchAddress("vertex_x_position", &primVx);
	inTree->SetBranchAddress("vertex_y_position", &primVy);
	inTree->SetBranchAddress("beam_momentum_x", &bx);
	inTree->SetBranchAddress("beam_momentum_y", &by);
	inTree->SetBranchAddress("beam_momentum_z", &bz);

	std::vector<TLorentzVector> particles;
//	particles.resize(6);
	particles.resize(4);
	unsigned int entries = inTree->GetEntries();

	std::vector<antok::beamfileGenerator::fiveDimCoord*>* tempTree = new std::vector<antok::beamfileGenerator::fiveDimCoord*>();

	for(unsigned int i = 0; i < entries; ++i) {

		inTree->GetEntry(i);

		vtxX->Fill(primVx);
		vtxY->Fill(primVy);
		momX->Fill(bx);
		momY->Fill(by);
		momZ->Fill(bz);

		tempTree->push_back(new antok::beamfileGenerator::fiveDimCoord(primVx, primVy, bx, by, bz));

		if(i % 100000 == 0) {
			std::cout<<"Entry "<<i<<" of "<<entries<<std::endl;
		}
	}

	std::vector<double> lowerCorner(5, 0);
	std::vector<double> upperCorner(5, 0);
	bool first = true;
	for(unsigned int i = 0; i < tempTree->size(); ++i) {
		primVx = (*tempTree)[i]->_coords[0];
		primVy = (*tempTree)[i]->_coords[1];
		bx = (*tempTree)[i]->_coords[2];
		by = (*tempTree)[i]->_coords[3];
		bz = (*tempTree)[i]->_coords[4];
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
	antok::beamfileGenerator::fiveDimBin bin(lowerCorner, upperCorner);
	std::cout<<"Got first bin: "<<std::endl;
	bin.print(std::cout);
	std::vector<double> scalingFactors(5, 0.);
	for(unsigned int i = 0; i < scalingFactors.size(); ++i) {
		scalingFactors[i] = bin._b[i] - bin._a[i];
	}
	std::cout<<"Scaling factors: ["<<scalingFactors[0];
	for(unsigned int i = 1; i < scalingFactors.size(); ++i) {
		std::cout<<", "<<scalingFactors[i];
	}
	std::cout<<"]"<<std::endl;

	std::list<std::pair<std::vector<antok::beamfileGenerator::fiveDimCoord*>*,
	          antok::beamfileGenerator::fiveDimBin> > adaptiveBins;
	antok::beamfileGenerator::getAdaptiveBins(adaptiveBins, bin, tempTree, 4, true);
	unsigned int nBins = adaptiveBins.size();
	std::cout<<"Split phase space in "<<nBins<<" bins."<<std::endl;

	double primVx_sigma, primVy_sigma, bx_sigma, by_sigma, bz_sigma;
	int binContent = 0;
	TTree* outTree = new TTree("beamTree", "beamTree");
	outTree->Branch("vertex_x_position", &primVx, "vertex_x_position/D");
	outTree->Branch("vertex_y_position", &primVy, "vertex_y_position/D");
	outTree->Branch("beam_momentum_x", &bx,"beam_momentum_x/D");
	outTree->Branch("beam_momentum_y", &by,"beam_momentum_y/D");
	outTree->Branch("beam_momentum_z", &bz,"beam_momentum_z/D");
	outTree->Branch("vertex_x_position_sigma", &primVx_sigma, "vertex_x_position_sigma/D");
	outTree->Branch("vertex_y_position_sigma", &primVy_sigma, "vertex_y_position_sigma/D");
	outTree->Branch("beam_momentum_x_sigma", &bx_sigma,"beam_momentum_x_sigma/D");
	outTree->Branch("beam_momentum_y_sigma", &by_sigma,"beam_momentum_y_sigma/D");
	outTree->Branch("beam_momentum_z_sigma", &bz_sigma,"beam_momentum_z_sigma/D");
	outTree->Branch("bin_content", &binContent, "bin_content/I");

	std::vector<double*> coords(5);
	std::vector<double*> sigmas(5);
	coords[0] = &primVx;
	coords[1] = &primVy;
	coords[2] = &bx;
	coords[3] = &by;
	coords[4] = &bz;
	sigmas[0] = &primVx_sigma;
	sigmas[1] = &primVy_sigma;
	sigmas[2] = &bx_sigma;
	sigmas[3] = &by_sigma;
	sigmas[4] = &bz_sigma;


	std::cout << "Calculating and saving sigmas." << std:: endl;

	unsigned binNumber = 0;
	for(
		std::list<std::pair<std::vector<antok::beamfileGenerator::fiveDimCoord*>*,
		          antok::beamfileGenerator::fiveDimBin> >::const_iterator binIt = adaptiveBins.begin();
		binIt != adaptiveBins.end();
		++binIt
	)
	{
		std::cout<<"Bin "<<++binNumber<<" of "<<nBins<<std::endl;
		const antok::beamfileGenerator::fiveDimBin& currentBin = binIt->second;
		std::vector<antok::beamfileGenerator::fiveDimCoord*>* currentTree = binIt->first;
		binContent = currentTree->size();
		const double binContentAsDouble = (double)binContent;

		double meanVolume = 1.;
		for(unsigned int i = 0; i < 5; ++i) {
			meanVolume *= (currentBin._b[i] - currentBin._a[i]) / scalingFactors[i];
		}
		meanVolume /= binContentAsDouble;
		const double meanVolumeEdgeLength = pow(meanVolume, 0.2);

		for(unsigned int i = 0; i < 5; ++i) {
			*(sigmas[i]) = meanVolumeEdgeLength * scalingFactors[i];
		}
		for(int i = 0; i < binContent; ++i) {
			for(unsigned int j = 0; j < 5; ++j) {
				*(coords[j]) = (*currentTree)[i]->_coords[j];
			}
			outTree->Fill();
		}
	}

	outFile->Write();
	outFile->Close();

}

int main(int argc, char* argv[]) {
	if(argc == 3) {
		fillFiveDimHist(argv[1], argv[2]);
	} else {
		std::cerr<<"Wrong number of arguments, is "<<argc<<", should be 3."<<std::endl;
	}
}
