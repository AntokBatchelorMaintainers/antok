
#include<iostream>

#include<TApplication.h>
#include<TFile.h>
#include<TTree.h>
#include<TStyle.h>

#include<constants.h>
#include<cutter.h>
#include<data.hpp>
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
	antok::Data& data = objectManager->getData();
	antok::Event& event = objectManager->getEvent();
	antok::Plotter& plotter = objectManager->getPlotter();
	TTree* inTree = objectManager->getInTree();

	TTree* out_tree = inTree->CloneTree(0);

	const unsigned int& N_PARTICLES = antok::Constants::n_particles();

	inTree->SetBranchAddress("Run", &data.Run);
	inTree->SetBranchAddress("TrigMask", &data.TrigMask);
	inTree->SetBranchAddress("EvNbr", &data.EvNbr);
	inTree->SetBranchAddress("SpillNbr", &data.SpillNbr);

	inTree->SetBranchAddress("X_primV", &data.X_primV);
	inTree->SetBranchAddress("Y_primV", &data.Y_primV);
	inTree->SetBranchAddress("Z_primV", &data.Z_primV);

	inTree->SetBranchAddress("gradx", &data.gradx);
	inTree->SetBranchAddress("grady", &data.grady);

	inTree->SetBranchAddress("beam_time", &data.beam_time);

/*	for(unsigned int i = 0; i < N_PARTICLES; ++i) {
		std::stringstream sstr;
		sstr<<"Mom_x"<<(i + 1);
		inTree->SetBranchAddress(sstr.str().c_str(), &data.Mom_x.at(i));
		sstr.str("");
		sstr<<"Mom_y"<<(i + 1);
		inTree->SetBranchAddress(sstr.str().c_str(), &data.Mom_y.at(i));
		sstr.str("");
		sstr<<"Mom_z"<<(i + 1);
		inTree->SetBranchAddress(sstr.str().c_str(), &data.Mom_z.at(i));
		sstr.str("");
		sstr<<"theta_RICH_"<<(i + 1);
		inTree->SetBranchAddress(sstr.str().c_str(), &data.theta_RICH.at(i));
		sstr.str("");
		sstr<<"PID_RICH_"<<(i + 1);
		inTree->SetBranchAddress(sstr.str().c_str(), &data.PID_RICH.at(i));
		sstr.str("");
		sstr<<"zmax"<<(i + 1);
		inTree->SetBranchAddress(sstr.str().c_str(), &data.z_max.at(i));
	}
*/
	inTree->SetBranchAddress("chi2PV", &data.chi2PV);

	inTree->SetBranchAddress("RPD_Px", &data.RPD_Px);
	inTree->SetBranchAddress("RPD_Py", &data.RPD_Py);
	inTree->SetBranchAddress("RPD_Pz", &data.RPD_Pz);
	inTree->SetBranchAddress("RPD_E", &data.RPD_E);
	inTree->SetBranchAddress("RPD_Tz", &data.RPD_Tz);
	inTree->SetBranchAddress("RPD_z", &data.RPD_z);
	inTree->SetBranchAddress("RPD_beta", &data.RPD_beta);
	inTree->SetBranchAddress("RPD_Phi", &data.RPD_Phi);
	inTree->SetBranchAddress("RPD_dEA", &data.RPD_dEA);
	inTree->SetBranchAddress("RPD_dEB", &data.RPD_dEB);
	inTree->SetBranchAddress("nbrRPDTracks", &data.nbrRPDTracks);

	inTree->SetBranchAddress("cedarID_bayes", &data.cedarID_bayes);
	inTree->SetBranchAddress("cedarTheta_X", &data.cedarTheta_X);
	inTree->SetBranchAddress("cedarTheta_Y", &data.cedarTheta_Y);
	inTree->SetBranchAddress("cedarProbK1", &data.cedarProbK1);
	inTree->SetBranchAddress("cedarProbK2", &data.cedarProbK2);
	inTree->SetBranchAddress("cedarProbK3", &data.cedarProbK3);

	TH1D* stats_pre = (TH1D*)infile->Get("kbicker/statistic");
	TH1D* stats = (TH1D*)stats_pre->Clone("statistics");

	assert(cutter.set_stats_histogram(stats));

	for(unsigned int i = 0; i < inTree->GetEntries(); ++i) {

		inTree->GetEntry(i);

		event.update(data);

		int cutmask = cutter.get_cutmask(event);
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
