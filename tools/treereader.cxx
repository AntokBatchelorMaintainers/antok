
#include <cassert>
#include <csignal>
#include <iostream>
#include <vector>

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

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"

#include "constants.h"
#include "cutter.h"
#include "event.h"
#include "initializer.h"
#include "object_manager.h"
#include "plotter.h"


static bool ABORT = false;


void
signal_handler(int signum)
{
	ABORT = true;
}


int
treereader(const std::vector<const char*>& inFileNames,
           const char*                     outFileName    = nullptr,
           const std::string&              configFileName = "../config/default.yaml")
{
	TApplication* app = new TApplication("app", 0, 0);
	gStyle->SetPalette(1);
	gStyle->SetCanvasColor(10);
	gStyle->SetPadColor(10);

	TFile* outFile = nullptr;
	if (outFileName == nullptr) {
		outFile = TFile::Open("out_tree.root", "RECREATE");
	} else {
		outFile = TFile::Open(outFileName, "NEW");
	}
	if (outFile == nullptr) {
		std::cerr << "Could not open output file. Aborting..." << std::endl;
		exit(1);
	}
	std::cout << "Successfully opened output file '" << outFile->GetPath() << "' for writing." << std::endl;

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	assert(objectManager->setOutFile(outFile));
	antok::Initializer* initializer = antok::Initializer::instance();
	if (not initializer->readConfigFile(configFileName)) {
		std::cerr << "Could not read config file. Aborting..." << std::endl;
		exit(1);
	}
	std::cout << "Successfully read config file '" << configFileName << "'." << std::endl;

	bool firstFile = true;
	for (const auto& inFileName : inFileNames) {
		if (ABORT) {
			break;
		}

		// open input files
		TFile* inFile = nullptr;
		if (inFileName == nullptr) {
			//TODO add meaningful default; best to add test case(s) with root and config files to Antok repo
			inFile = TFile::Open("/afs/cern.ch/user/k/kbicker/w0/analysis/phast/5Pi_fhaasUE_10chunks.root", "READ");
		} else {
			inFile = TFile::Open(inFileName, "READ");
		}
		if (inFile == nullptr) {
			std::cerr << "Could not open input file";
			if (inFileName != nullptr) {
				std::cerr << " '" << inFileName << "'";
			}
			std::cerr << ". Aborting..." << std::endl;
			exit(1);
		}
		std::cout << "Processing input file '" << inFile->GetPath() << "'" << std::endl;

		if (firstFile) {
			firstFile = false;
			// first initialization
			assert(objectManager->setInFile(inFile));
			if (not initializer->initAll()) {
				std::cerr << "Error while initializing. Aborting..." << std::endl;
				exit(1);
			}
			std::cout << "Initialization successful." << std::endl;
		} else {
			assert(objectManager->changeInFile(inFile));
			if (not initializer->updateInput()) {
				std::cerr << "Error while initializing input tree. Aborting..." << std::endl;
				exit(1);
			}
			std::cout << "Initialization of input tree successful." << std::endl;
		}

		TTree* inTree = objectManager->getInTree();
		if (inTree == nullptr) {
			std::cerr << "Error cannot get input tree. Aborting..." << std::endl;
			exit(1);
		}
		std::cout << "Reading from input tree '" << inTree->GetName() << "'." << std::endl;
		progress_display progressIndicator(inTree->GetEntries(), std::cout, "");
		for (unsigned int i = 0; i < inTree->GetEntries(); ++i) {
			if (ABORT) {
				const double percent = 100. * ((double)i / (double)(inTree->GetEntries()));
				std::cout << "At event " << i << " of " << inTree->GetEntries() << " (" << percent << "%)." << std::endl;
				std::cout << "Caught CTRL-C, aborting..." << std::endl;
				break;
			}

			inTree->GetEntry(i);
			if (not objectManager->magic()) {
				std::cerr << "Could not process event " << i << ". Aborting..." << std::endl;
				exit(1);
			}
			++progressIndicator;
		}
		std::cout << "Successfully read " << inTree->GetEntries() << " events from input tree '" << inTree->GetName() << "'." << std::endl;

	}  // loop over input files

	if (not objectManager->finish()) {
		std::cerr << "Problem when writing TObjects and/or closing output file." << std::endl;
		exit(1);
	}
	delete app;
	return 0;
}


int
treereader(const char*        inFileName     = nullptr,
           const char*        outFileName    = nullptr,
           const std::string& configFileName = "../config/default.yaml")
{
	const std::vector<const char*> inFileNames(1, inFileName);
	return treereader(inFileNames, outFileName, configFileName);
}


int
main(int argc,
     char* argv[])
{
	signal(SIGINT, signal_handler);

	if (argc == 1) {
		return treereader();
	} else if (argc == 3) {
		return treereader(argv[1], argv[2]);
	} else if (argc == 4) {
		return treereader(argv[1], argv[2], argv[3]);
	} else {
		// multiple input files
		std::vector<const char*> inFileNames(argc - 3);
		for (int i = 1; i < argc - 2; ++i) {  // loop over input files
			inFileNames[i - 1] = argv[i];
		}
		return treereader(inFileNames, argv[argc - 2], argv[argc - 1]);
	}
}
