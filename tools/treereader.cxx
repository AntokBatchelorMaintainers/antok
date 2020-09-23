
#include<iostream>
#include<vector>
#include<signal.h>

#include <boost/version.hpp>
// progress_display was moved to different header file and name space in version 1.72.0
// in 1.71.0 deprecation message is printed but there is no progress_display in <boost/timer/timer.hpp> (sigh)
// see https://github.com/boostorg/timer/issues/12
#if (BOOST_VERSION >= 107200)
#include <boost/timer/progress_display.hpp>
using boost::timer::progress_display;
#else
#if (BOOST_VERSION == 107100)
// suppress deprecation message
#define BOOST_ALLOW_DEPRECATED_HEADERS
#endif
#include <boost/progress.hpp>
using boost::progress_display;
#endif

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

int treereader(std::vector<const char*> infilenames, char* outfilename=0, std::string configfilename = "../config/default.yaml") {

	new TApplication("app", 0, 0);

	gStyle->SetPalette(1);
	gStyle->SetCanvasColor(10);
	gStyle->SetPadColor(10);

	TFile* outfile;
	if(outfilename != 0) {
		outfile = TFile::Open(outfilename, "NEW");
	} else {
		outfile = TFile::Open("out_tree.root", "RECREATE");
	}
	if(outfile == 0) {
		return 1;
	}

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	assert(objectManager->setOutFile(outfile));

	antok::Initializer* initializer = antok::Initializer::instance();
	if(not initializer->readConfigFile(configfilename)) {
		std::cerr<<"Could not open config file. Aborting..."<<std::endl;
		exit(1);
	}

	for(std::vector<const char*>::const_iterator infilename = infilenames.begin(); infilename != infilenames.end(); ++infilename ){
		if(ABORT) break;

		TFile* infile;
		if(*infilename != 0) {
			infile = TFile::Open(*infilename, "READ");
		} else {
			infile = TFile::Open("/afs/cern.ch/user/k/kbicker/w0/analysis/phast/5Pi_fhaasUE_10chunks.root", "READ");
		}
		if(infile == 0) {
			std::cerr<<"Could not open input file. Aborting..."<<std::endl;
			return 1;
		}


		if( infilename == infilenames.begin() ){ // first initialization
			assert(objectManager->setInFile(infile));
			if(not initializer->initAll()) {
				std::cerr<<"Error while initializing. Aborting..."<<std::endl;
				exit(1);
			}
		} else {
			assert( objectManager->changeInFile(infile) );
			if(not initializer->updateInput()) {
				std::cerr<<"Error while initializing input tree. Aborting..."<<std::endl;
				exit(1);
			}
		}

		std::cout << "Processing input file '" << *infilename << "'" << std::endl;
		TTree* inTree = objectManager->getInTree();

		progress_display* progressIndicator = new progress_display(inTree->GetEntries(), std::cout, "");

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

			++(*progressIndicator);

		}


	}


	if(not objectManager->finish()) {
		std::cerr<<"Problem when writing TObjects and/or closing output file."<<std::endl;
		return 1;
	}
	return 0;

}
int treereader(char* infilename=0, char* outfilename=0, std::string configfilename = "../config/default.yaml") {
	std::vector<const char*> inputfiles(1, infilename);
	return treereader(inputfiles, outfilename, configfilename);
}

int main(int argc, char* argv[]) {

	signal(SIGINT, signal_handler);

	if(argc == 1) {
		return treereader();
	} else if (argc == 3) {
		return treereader(argv[1], argv[2]);
	} else if (argc == 4) {
		return treereader(argv[1], argv[2], argv[3]);
	} else { // multiple input files
		std::vector<const char*> input_files( argc - 3 );
		for( int i = 1; i < argc - 2; ++i ) // loop over input files
			input_files[i-1] = argv[i];
		return treereader(input_files, argv[argc-2], argv[argc-1]);
	}
}
