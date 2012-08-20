
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

#include<assert.h>

void treereader(char* infilename=NULL, char* outfilename=NULL) {

	using hlib::PION_MASS;
	using hlib::PROTON_MASS;

	TApplication* app = new TApplication("app", 0, NULL);

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
	hlib::Cutter cutter;

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


	std::vector<TH1*> hists;

	TH1D* stats_pre = (TH1D*)infile->Get("kbicker/statistic");
//	TH1D* stats_pre = (TH1D*)infile->Get("fhaas/statistic");
	TH1D* stats = (TH1D*)stats_pre->Clone("statistics");
	hists.push_back(stats);
	TH1D* mass_5pi = new TH1D("mass_5pi", "mass_5Pi", 500, 0, 7);
	hists.push_back(mass_5pi);
	TH1D* mom_5pi = new TH1D("mom_5pi", "mom_5Pi", 500, 0, 250);
	hists.push_back(mom_5pi);
	TH1D* mom_5pi_raw = new TH1D("mom_5pi_raw", "mom_5Pi_raw", 500, 0, 250);
	hists.push_back(mom_5pi_raw);
	TH1D* calc_beam_E = new TH1D("calc_beam_E", "calc_beam_E", 500, 0, 250);
	hists.push_back(calc_beam_E);
	TH1D* rpd_mult = new TH1D("rpd_mult", "rpd_mult", 10, 0, 10);
	hists.push_back(rpd_mult);
	TH1D* rpd_Pxh = new TH1D("rpd_Px", "rpd_Px", 1000, 0, 10);
	hists.push_back(rpd_Pxh);
	TH2D* vtx_pos = new TH2D("vtx_pos", "vtx_pos", 1000, -5, 5, 1000, -5, 5);
	vtx_pos->SetDrawOption("colz");
	hists.push_back(vtx_pos);
	TH1D* vtx_zh = new TH1D("vtx_z", "vtx_z", 2000, -200, 200);
	hists.push_back(vtx_zh);
	TH1D* t_primh = new TH1D("t_prime", "t_prime", 1000, -5, 5);
	hists.push_back(t_primh);
	TH1D* delta_phih = new TH1D("delta_phi", "delta_phi", 500, -7, 7);
	hists.push_back(delta_phih);
	TH1D* trig_maskh = new TH1D("trigger_mask", "trigger_mask", 15, 0, 15);
	hists.push_back(trig_maskh);

int test = 0;

	for(unsigned int i = 0; i < tree_chain->GetEntries(); ++i) {

		tree_chain->GetEntry(i);

		event.update(data);

		int cutmask = cutter.get_cutmask(event);
		if(cutmask == 0) {
			out_tree->Fill();
		}

assert(test == 0);
if(cutter.get_cutmask(event) == 0) {
	test = 1;
}

		stats->Fill("All events", 1);

		trig_maskh->Fill(data.TrigMask);
		if(!(data.TrigMask&0x1)) {
			continue;
		}
		stats->Fill("Trigger Mask = 1", 1);

		vtx_zh->Fill(data.Z_primV);
		if((data.Z_primV > -28.4) || (data.Z_primV < -68.4)) {
			continue;
		}
		stats->Fill("Vertex z in ]-28.4,-68.4[", 1);

		vtx_pos->Fill(data.X_primV, data.Y_primV);
		if(std::pow(data.X_primV, 2) + std::pow(data.Y_primV, 2) > 3.0625) {
			continue;
		}
		stats->Fill("Vertex.R() < 1.75", 1);

		rpd_mult->Fill(data.nbrRPDTracks);
		if(data.nbrRPDTracks != 1) {
			continue;
		}
		stats->Fill("1 RPD track", 1);

		rpd_Pxh->Fill(event.get_pProton().M());
		if(event.get_pProton().M() < 0.2) {
			continue;
		}
		stats->Fill("Proton mass > 0.2", 1);

		if(data.isKaon != 0) {
			continue;
		}
		stats->Fill("data.isKaon = 0", 1);

		t_primh->Fill(event.get_tPrime());

		if(event.get_tPrime() < 0.1) {
			continue;
		}
		stats->Fill("T-prime > 0.1", 1);

		mom_5pi_raw->Fill(event.get_pSum().Energy());

		delta_phih->Fill(event.get_RpdDeltaPhi());
		if(std::fabs(event.get_RpdDeltaPhi()) > event.get_RpdPhiRes()) {
			continue;
		}
		stats->Fill("RPD planarity cut", 1);

		mom_5pi->Fill(event.get_pSum().Energy());
		calc_beam_E->Fill(event.get_pBeam().E());
		if(std::fabs(event.get_pSum().Energy()-191.) > 3.28) {
			continue;
		}
		stats->Fill("Exclusivity 191+-3.28GeV", 1);

		mass_5pi->Fill(event.get_pSum().M());

assert(cutter.get_cutmask(event) == 0);
assert(test == 1);
test = 0;

//		out_tree->Fill();

	}

	outfile->cd();
	out_tree->Write();
	outfile->mkdir("Histograms");
	outfile->Cd("Histograms");
	for(unsigned int i = 0; i < hists.size(); ++i) {
		hists.at(i)->Write();
	}

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
