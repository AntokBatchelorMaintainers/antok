
#include<iostream>

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
	objectManager->setInFile(infile);
	antok::Initializer* initializer = antok::Initializer::instance();
	if(not initializer->readConfigFile(configfilename)) {
		std::cerr<<"Could not open config file. Aborting..."<<std::endl;
		exit(1);
	}

	if(not initializer->initAll()) {
		std::cerr<<"Error while initializing. Aborting..."<<std::endl;
		exit(1);
	}
	antok::Cutter& cutter = objectManager->getCutter();
	antok::Event& event = objectManager->getEvent();
	antok::Plotter& plotter = objectManager->getPlotter();
	TTree* inTree = objectManager->getInTree();

	TTree* out_tree = inTree->CloneTree(0);

	TH1D* stats_pre = (TH1D*)infile->Get("kbicker/statistic");
	TH1D* stats = (TH1D*)stats_pre->Clone("statistics");

	assert(cutter.set_stats_histogram(stats));

	for(unsigned int i = 0; i < inTree->GetEntries(); ++i) {

		inTree->GetEntry(i);

		event.update();

		int cutmask = cutter.get_cutmask();
		if(cutmask == 0) {
			out_tree->Fill();
		}

		plotter.fill(event, cutmask);

	}

	plotter.save(outfile);
	outfile->cd();
	out_tree->Write();
	(cutter.get_stats_histogram())->Write();

	infile->Close();
	outfile->Close();

}

int main(int argc, char* argv[]) {
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
