
#include<iomanip>
#include<iostream>
#include<string>

#include<boost/shared_ptr.hpp>

#include<TFile.h>
#include<TH1D.h>
#include<TLorentzVector.h>
#include<TTree.h>

#include<basic_calcs.h>
#include<beamfile_generator_helpers.h>
#include<beamfile_generator_5dBin.h>
#include<beamfile_generator_5dCoord.h>
#include<constants.h>
#include<initializer.h>


void fillFiveDimHist(std::string inFileName, std::string outFileName, int startCoord = 0) {

	if(startCoord < 0 or startCoord > 4) {
		std::cout<<"Invalid starting coordinate: "<<startCoord<<". Aborting..."<<std::endl;
		exit(1);
	}

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

		tempTree->push_back(new antok::beamfileGenerator::fiveDimCoord(primVx, primVy, bx, by, bz, i));

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
	std::list<boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> > adaptiveBins;
	{
		boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> bin(new antok::beamfileGenerator::fiveDimBin (lowerCorner,
		                                                                                                      upperCorner,
		                                                                                                      tempTree));
		bin->setOnLowerEdge(std::vector<bool>(5, true));
		bin->setOnUpperEdge(std::vector<bool>(5, true));
		std::cout<<"Got first bin: "<<std::endl;
		bin->print(std::cout);

		std::cout<<"starting binning with coordinate "<<startCoord<<"."<<std::endl;
		antok::beamfileGenerator::getAdaptiveBins(adaptiveBins, bin, startCoord, false);
	}
	const unsigned int nBins = adaptiveBins.size();
	std::cout<<"Split phase space in "<<nBins<<" bins."<<std::endl;
	std::cout<<"(self-reporting of the bin class gives "
	         <<antok::beamfileGenerator::fiveDimBin::getNExistingBins()<<" bins)" <<std::endl;

	double primVx_sigma, primVy_sigma, bx_sigma, by_sigma, bz_sigma;
	int binContent, nNeighbors, edgeity, sigmaCalculationMethod;
	long eventNumber;
	double binVolume;
	TTree* outTree = new TTree("beamTree", "beamTree");
	outTree->Branch("vertex_x_position", &primVx, "vertex_x_position/D");
	outTree->Branch("vertex_y_position", &primVy, "vertex_y_position/D");
	outTree->Branch("beam_momentum_x", &bx,"beam_momentum_x/D");
	outTree->Branch("beam_momentum_y", &by,"beam_momentum_y/D");
	outTree->Branch("beam_momentum_z", &bz,"beam_momentum_z/D");
	outTree->Branch("vertex_x_position_sigma", &primVx_sigma, "vertex_x_position_sigma/D");
	outTree->Branch("vertex_y_position_sigma", &primVy_sigma, "vertex_y_position_sigma/D");
	outTree->Branch("beam_momentum_x_sigma", &bx_sigma, "beam_momentum_x_sigma/D");
	outTree->Branch("beam_momentum_y_sigma", &by_sigma, "beam_momentum_y_sigma/D");
	outTree->Branch("beam_momentum_z_sigma", &bz_sigma, "beam_momentum_z_sigma/D");
	outTree->Branch("event_number", &eventNumber, "event_number/L");
	outTree->Branch("bin_content", &binContent, "bin_content/I");
	outTree->Branch("number_of_neighbors", &nNeighbors, "number_of_neighbors/I");
	outTree->Branch("edgeity", &edgeity, "edgeity/I");
	outTree->Branch("bin_volume", &binVolume,"bin_volume/D");
	outTree->Branch("sigma_calculation_method", &sigmaCalculationMethod, "sigma_calculation_method/I");

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

	std::cout << "Calculating sigmas." << std:: endl;

	unsigned int binNumber = 0;
	const unsigned int roundingNumber = int(std::pow(10., (unsigned int)(log10((double)nBins / 100.) + 0.5)) + 0.5);

	antok::beamfileGenerator::fiveDimBin::setDebug(false);
	antok::beamfileGenerator::fiveDimBin::setPrintNeighbors(false);
	antok::beamfileGenerator::fiveDimBin::setDifferentSigmaCalculationForEdges(true);

	std::vector<antok::beamfileGenerator::eventBookkeeper> eventsToSave;

	for(
		std::list<boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> >::const_iterator binIt = adaptiveBins.begin();
		binIt != adaptiveBins.end();
		++binIt
	)
	{
		if(not (binNumber % roundingNumber)) {
			std::cout<<"Bin "<<binNumber<<" of "<<nBins<<" ("<<std::setprecision(2)
			         <<(binNumber/(double)nBins*100)<<"%)"<<std::endl;
		}
		++binNumber;
		const antok::beamfileGenerator::fiveDimBin& currentBin = *(*binIt);
		const std::vector<antok::beamfileGenerator::fiveDimCoord*>* currentTree = currentBin.getEvents();

		const boost::shared_ptr<const std::vector<std::vector<double> > > sigmasFromBin = currentBin.getSigmas();
		for(unsigned int i = 0; i < currentTree->size(); ++i) {
			assert((*currentTree)[i]->_eventNumber >= 0);
			antok::beamfileGenerator::eventBookkeeper currentEvent;
			currentEvent.binContent = currentTree->size();
			currentEvent.binVolume = currentBin.getVolume();
			currentEvent.nNeighbors = currentBin.getNeighbors().size();
			currentEvent.edgeity = currentBin.getEdgeity();
			currentEvent.sigmaCalculationMethod = currentBin.getSigmaCalculationMethod();
			currentEvent.sigmas = sigmasFromBin;
			currentEvent.sigmaIndex = i;
			currentEvent.coords = (*currentTree)[i];
			eventsToSave.push_back(currentEvent);
		}
	}
	for(
		std::list<boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> >::const_iterator binIt = adaptiveBins.begin();
		binIt != adaptiveBins.end();
		++binIt
	)
	{
		(*binIt)->clearSigmaCache();
	}

	std::cout << "Sorting events." << std::endl;
	std::sort(eventsToSave.begin(), eventsToSave.end(), antok::beamfileGenerator::compareEventBookkeepers());

	std::cout << "Saving events." << std:: endl;
	const unsigned int eventRoundingNumber = int(std::pow(10., (unsigned int)(log10((double)eventsToSave.size() / 100.) + 0.5)) + 0.5);
	for(unsigned int i = 0; i < eventsToSave.size(); ++i) {
		if(not (i % eventRoundingNumber)) {
			std::cout<<"Bin "<<i<<" of "<<eventsToSave.size()<<" ("<<std::setprecision(2)
			         <<(i/(double)eventsToSave.size()*100)<<"%)"<<std::endl;
		}
		const antok::beamfileGenerator::eventBookkeeper& currentEvent = eventsToSave[i];
		binContent = currentEvent.binContent;
		binVolume = currentEvent.binVolume;
		nNeighbors = currentEvent.nNeighbors;
		edgeity = currentEvent.edgeity;
		sigmaCalculationMethod = currentEvent.sigmaCalculationMethod;
		for(unsigned int j = 0; j < 5; ++j) {
			*(sigmas[j]) = (*currentEvent.sigmas)[currentEvent.sigmaIndex][j];
			*(coords[j]) = currentEvent.coords->_coords[j];
		}
		eventNumber = currentEvent.coords->_eventNumber;
		outTree->Fill();
	}

	outFile->Write();
	outFile->Close();

}

int main(int argc, char* argv[]) {
	if(argc == 3) {
		fillFiveDimHist(argv[1], argv[2]);
	} else if (argc == 4) {
		fillFiveDimHist(argv[1], argv[2], std::atoi(argv[3]));
	} else {
		std::cerr<<"Wrong number of arguments, is "<<argc<<", should be 2 or 3."<<std::endl;
	}
}
