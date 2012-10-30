
#include<iostream>
#include<signal.h>

#include<TApplication.h>
#include<TFile.h>
#include<TTree.h>
#include<TStyle.h>

#include<constants.h>
#include<cutter.h>
#include<event.h>
#include<initializer.h>
#include<object_manager.h>
#include<plotter.h>

#include<assert.h>

static bool ABORT = false;

void signal_handler(int signum) {
	ABORT = true;
}

void treereader(char* infilename=0, char* outfilename=0, std::string configfilename = "../config/default.yaml") {

	new TApplication("app", 0, 0);

	gStyle->SetPalette(1);
	gStyle->SetCanvasColor(10);
	gStyle->SetPadColor(10);

	TFile* infile;
	if(infilename != 0) {
		infile = TFile::Open(infilename, "READ");
	} else {
		infile = TFile::Open("/afs/cern.ch/user/k/kbicker/w0/analysis/phast/5Pi_fhaasUE_10chunks.root", "READ");
// 		infile = TFile::Open("/afs/cern.ch/user/k/kbicker/scratch0/prefiltering_run1_merged/files_H_2008_26.root");
	}
	if(infile == 0) {
		std::cerr<<"Could not open input file. Aborting..."<<std::endl;
		return;
	}

	TFile* outfile;
	if(outfilename != 0) {
		outfile = TFile::Open(outfilename, "NEW");
	} else {
		outfile = TFile::Open("out_tree.root", "RECREATE");
	}
	if(outfile == 0) {
		return;
	}

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	assert(objectManager->setInFile(infile));
	assert(objectManager->setOutFile(outfile));

	antok::Initializer* initializer = antok::Initializer::instance();
	if(not initializer->readConfigFile(configfilename)) {
		std::cerr<<"Could not open config file. Aborting..."<<std::endl;
		exit(1);
	}

	if(not initializer->initAll()) {
		std::cerr<<"Error while initializing. Aborting..."<<std::endl;
		exit(1);
	}
	TTree* inTree = objectManager->getInTree();

	for(unsigned int i = 0; i < inTree->GetEntries(); ++i) {

		if(ABORT) {
			double percent = 100. * ((double)i / (double)(inTree->GetEntries()));
			std::cout<<"At event "<<i<<" of "<<inTree->GetEntries()<<" ("<<percent<<"%)."<<std::endl;
			std::cout<<"Caught CTRL-C, aborting..."<<std::endl;
			break;
		}

		inTree->GetEntry(i);

		if(not objectManager->magic()) {
			std::cerr<<"Could not process event "<<i<<". Aborting..."<<std::endl;
			exit(1);
		}

	}

	if(not objectManager->finish()) {
		std::cerr<<"Problem when writing TObjects and/or closing output file."<<std::endl;
	}

}

int main(int argc, char* argv[]) {

	signal(SIGINT, signal_handler);

	if(argc == 1) {
		treereader();
	} else if (argc == 3) {
		treereader(argv[1], argv[2]);
	} else if (argc == 4) {
		treereader(argv[1], argv[2], argv[3]);
	} else {
		std::cerr<<"Wrong number of arguments, is "<<argc<<", should be in [0, 2, 3]."<<std::endl;
	}
}
