
#include<iostream>
#include<sstream>

#include<TApplication.h>
#include<TChain.h>
#include<TFile.h>
#include<TMath.h>
#include<TStyle.h>
#include<TH1D.h>
#include<TH2D.h>
#include<TCanvas.h>
#include<TRotation.h>
#include<TLorentzVector.h>
#include<TLorentzRotation.h>
#include<TVector3.h>

#include<basic_calcs.h>
#include<constants.hpp>
#include<cutter.h>
#include<data.hpp>
#include<event.h>
#include<plotter.h>

#include<assert.h>

void treereader(char* infilename=NULL, char* outfilename=NULL) {

	using hlib::PION_MASS;
	using hlib::PROTON_MASS;

	new TApplication("app", 0, NULL);

	gStyle->SetPalette(1);
	gStyle->SetCanvasColor(10);
	gStyle->SetPadColor(10);

	TFile* infile;
	if(infilename != NULL) {
		infile = TFile::Open(infilename, "READ");
	} else {
		infile = TFile::Open("/afs/cern.ch/user/k/kbicker/w0/analysis/phast/5Pi_fhaasUE_10chunks.root", "READ");
// 		infile = TFile::Open("/afs/cern.ch/user/k/kbicker/scratch0/prefiltering_run1_merged/files_H_2008_26.root");
	}
	if(infile == NULL) {
		return;
	}

	TTree* tree_chain = (TTree*)infile->Get("kbicker/USR55");
//	TTree* tree_chain = (TTree*)infile->Get("fhaas/USR52");
	if(tree_chain == NULL) {
		std::cout<<"Error opening in-TTree."<<std::endl;
		return;
	}

	TFile* outfile;
	if(outfilename != NULL) {
		outfile = TFile::Open(outfilename, "NEW");
	} else {
		outfile = TFile::Open("out_tree.root", "RECREATE");
	}
	if(outfile == NULL) {
		return;
	}

	TTree* out_tree = tree_chain->CloneTree(0);

	hlib::Data data;
	hlib::Event* event = hlib::Event::instance();
	hlib::Cutter* cutter = hlib::Cutter::instance();
	hlib::Plotter* plotter = hlib::Plotter::instance();

	tree_chain->SetBranchAddress("Run", &data.Run);
	tree_chain->SetBranchAddress("TrigMask", &data.TrigMask);
	tree_chain->SetBranchAddress("EvNbr", &data.EvNbr);
	tree_chain->SetBranchAddress("SpillNbr", &data.SpillNbr);

	tree_chain->SetBranchAddress("X_primV", &data.X_primV);
	tree_chain->SetBranchAddress("Y_primV", &data.Y_primV);
	tree_chain->SetBranchAddress("Z_primV", &data.Z_primV);

	tree_chain->SetBranchAddress("gradx", &data.gradx);
	tree_chain->SetBranchAddress("grady", &data.grady);

	for(unsigned int i = 0; i < hlib::N_PARTICLES; ++i) {
		std::stringstream sstr;
		sstr<<"Mom_x"<<(i + 1);
		tree_chain->SetBranchAddress(sstr.str().c_str(), &data.Mom_x.at(i));
		sstr.str("");
		sstr<<"Mom_y"<<(i + 1);
		tree_chain->SetBranchAddress(sstr.str().c_str(), &data.Mom_y.at(i));
		sstr.str("");
		sstr<<"Mom_z"<<(i + 1);
		tree_chain->SetBranchAddress(sstr.str().c_str(), &data.Mom_z.at(i));
		sstr.str("");
		sstr<<"theta_RICH_"<<(i + 1);
		tree_chain->SetBranchAddress(sstr.str().c_str(), &data.theta_RICH.at(i));
		sstr.str("");
		sstr<<"PID_RICH_"<<(i + 1);
		tree_chain->SetBranchAddress(sstr.str().c_str(), &data.PID_RICH.at(i));
		sstr.str("");
		sstr<<"zmax"<<(i + 1);
		tree_chain->SetBranchAddress(sstr.str().c_str(), &data.z_max.at(i));
	}

	tree_chain->SetBranchAddress("chi2PV", &data.chi2PV);

	tree_chain->SetBranchAddress("RPD_Px", &data.RPD_Px);
	tree_chain->SetBranchAddress("RPD_Py", &data.RPD_Py);
	tree_chain->SetBranchAddress("RPD_Pz", &data.RPD_Pz);
	tree_chain->SetBranchAddress("RPD_E", &data.RPD_E);
	tree_chain->SetBranchAddress("RPD_Tz", &data.RPD_Tz);
	tree_chain->SetBranchAddress("RPD_z", &data.RPD_z);
	tree_chain->SetBranchAddress("RPD_beta", &data.RPD_beta);
	tree_chain->SetBranchAddress("RPD_Phi", &data.RPD_Phi);
	tree_chain->SetBranchAddress("RPD_dEA", &data.RPD_dEA);
	tree_chain->SetBranchAddress("RPD_dEB", &data.RPD_dEB);
	tree_chain->SetBranchAddress("nbrRPDTracks", &data.nbrRPDTracks);

	tree_chain->SetBranchAddress("isKaon", &data.isKaon);

	TH1D* stats_pre = (TH1D*)infile->Get("kbicker/statistic");
	TH1D* stats = (TH1D*)stats_pre->Clone("statistics");

	assert(cutter->set_stats_histogram(stats));

	for(unsigned int i = 0; i < tree_chain->GetEntries(); ++i) {

		tree_chain->GetEntry(i);

		event->update(data);

		int cutmask = cutter->get_cutmask(*event);
		if(cutmask == 0) {
			out_tree->Fill();
		}

		plotter->fill(*event, cutmask);

	}

	plotter->save(outfile);
	outfile->cd();
	out_tree->Write();
	(cutter->get_stats_histogram())->Write();

	infile->Close();
	outfile->Close();

}

int main(int argc, char* argv[]) {
	if(argc != 3) {
		treereader();
	} else {
		treereader(argv[1], argv[2]);
	}
}
