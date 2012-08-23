
#include<iostream>

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
#include<constants.h>
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
	hlib::Event event;
	hlib::Cutter* cutter = hlib::Cutter::instance();
	hlib::Plotter plotter;

	tree_chain->SetBranchAddress("Run", &data.Run);
	tree_chain->SetBranchAddress("TrigMask", &data.TrigMask);
	tree_chain->SetBranchAddress("EvNbr", &data.EvNbr);
	tree_chain->SetBranchAddress("SpillNbr", &data.SpillNbr);

	tree_chain->SetBranchAddress("X_primV", &data.X_primV);
	tree_chain->SetBranchAddress("Y_primV", &data.Y_primV);
	tree_chain->SetBranchAddress("Z_primV", &data.Z_primV);

	tree_chain->SetBranchAddress("gradx", &data.gradx);
	tree_chain->SetBranchAddress("grady", &data.grady);

	tree_chain->SetBranchAddress("Mom_x1", &data.Mom_x1);
	tree_chain->SetBranchAddress("Mom_x2", &data.Mom_x2);
	tree_chain->SetBranchAddress("Mom_x3", &data.Mom_x3);
	tree_chain->SetBranchAddress("Mom_x4", &data.Mom_x4);
	tree_chain->SetBranchAddress("Mom_x5", &data.Mom_x5);

	tree_chain->SetBranchAddress("Mom_y1", &data.Mom_y1);
	tree_chain->SetBranchAddress("Mom_y2", &data.Mom_y2);
	tree_chain->SetBranchAddress("Mom_y3", &data.Mom_y3);
	tree_chain->SetBranchAddress("Mom_y4", &data.Mom_y4);
	tree_chain->SetBranchAddress("Mom_y5", &data.Mom_y5);

	tree_chain->SetBranchAddress("Mom_z1", &data.Mom_z1);
	tree_chain->SetBranchAddress("Mom_z2", &data.Mom_z2);
	tree_chain->SetBranchAddress("Mom_z3", &data.Mom_z3);
	tree_chain->SetBranchAddress("Mom_z4", &data.Mom_z4);
	tree_chain->SetBranchAddress("Mom_z5", &data.Mom_z5);

	tree_chain->SetBranchAddress("chi2PV", &data.chi2PV);

	tree_chain->SetBranchAddress("theta_RICH_1", &data.theta_RICH_1);
	tree_chain->SetBranchAddress("theta_RICH_2", &data.theta_RICH_2);
	tree_chain->SetBranchAddress("theta_RICH_3", &data.theta_RICH_3);
	tree_chain->SetBranchAddress("theta_RICH_4", &data.theta_RICH_4);
	tree_chain->SetBranchAddress("theta_RICH_5", &data.theta_RICH_5);

	tree_chain->SetBranchAddress("PID_RICH_1", &data.PID_RICH_1);
	tree_chain->SetBranchAddress("PID_RICH_2", &data.PID_RICH_2);
	tree_chain->SetBranchAddress("PID_RICH_3", &data.PID_RICH_3);
	tree_chain->SetBranchAddress("PID_RICH_4", &data.PID_RICH_4);
	tree_chain->SetBranchAddress("PID_RICH_5", &data.PID_RICH_5);

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

	tree_chain->SetBranchAddress("zmax1", &data.zmax1);
	tree_chain->SetBranchAddress("zmax2", &data.zmax2);
	tree_chain->SetBranchAddress("zmax3", &data.zmax3);
	tree_chain->SetBranchAddress("zmax4", &data.zmax4);
	tree_chain->SetBranchAddress("zmax5", &data.zmax5);

	TH1D* stats_pre = (TH1D*)infile->Get("kbicker/statistic");
	TH1D* stats = (TH1D*)stats_pre->Clone("statistics");

	assert(cutter->set_stats_histogram(stats));

	for(unsigned int i = 0; i < tree_chain->GetEntries(); ++i) {

		tree_chain->GetEntry(i);

		event.update(data);

		int cutmask = cutter->get_cutmask(event);
		if(cutmask == 0) {
			out_tree->Fill();
		}

		plotter.fill(event, cutmask);

	}

	plotter.save(outfile);
	outfile->cd();
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
