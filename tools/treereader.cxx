
#include<iostream>

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

#include<constants.h>
#include<basic_calcs.h>

void treereader(char* infilename=NULL, char* outfilename=NULL) {

	using hlib::PION_MASS;
	using hlib::PROTON_MASS;

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

	double px1, py1, pz1;
	TLorentzVector p1;
	double px2, py2, pz2;
	TLorentzVector p2;
	double px3, py3, pz3;
	TLorentzVector p3;
	double px4, py4, pz4;
	TLorentzVector p4;
	double px5, py5, pz5;
	TLorentzVector p5;

	TLorentzVector pTot;
	double gradx, grady;

	int nbrRPDTracks;
	double rpd_E, rpd_Px, rpd_Py, rpd_Pz, rpd_Phi;
	TLorentzVector proton;

	double vtx_x, vtx_y, vtx_z;
	int trigMask;

	int isKaon;

	tree_chain->SetBranchAddress("Mom_x1", &px1);
	tree_chain->SetBranchAddress("Mom_y1", &py1);
	tree_chain->SetBranchAddress("Mom_z1", &pz1);
	tree_chain->SetBranchAddress("Mom_x2", &px2);
	tree_chain->SetBranchAddress("Mom_y2", &py2);
	tree_chain->SetBranchAddress("Mom_z2", &pz2);
	tree_chain->SetBranchAddress("Mom_x3", &px3);
	tree_chain->SetBranchAddress("Mom_y3", &py3);
	tree_chain->SetBranchAddress("Mom_z3", &pz3);
	tree_chain->SetBranchAddress("Mom_x4", &px4);
	tree_chain->SetBranchAddress("Mom_y4", &py4);
	tree_chain->SetBranchAddress("Mom_z4", &pz4);
	tree_chain->SetBranchAddress("Mom_x5", &px5);
	tree_chain->SetBranchAddress("Mom_y5", &py5);
	tree_chain->SetBranchAddress("Mom_z5", &pz5);
	tree_chain->SetBranchAddress("nbrRPDTracks", &nbrRPDTracks);
	tree_chain->SetBranchAddress("RPD_Px", &rpd_Px);
	tree_chain->SetBranchAddress("RPD_Py", &rpd_Py);
	tree_chain->SetBranchAddress("RPD_Pz", &rpd_Pz);
	tree_chain->SetBranchAddress("RPD_E", &rpd_E);
	tree_chain->SetBranchAddress("RPD_Phi", &rpd_Phi);
	tree_chain->SetBranchAddress("X_primV", &vtx_x);
	tree_chain->SetBranchAddress("Y_primV", &vtx_y);
	tree_chain->SetBranchAddress("Z_primV", &vtx_z);
	tree_chain->SetBranchAddress("TrigMask", &trigMask);
	tree_chain->SetBranchAddress("gradx", &gradx);
	tree_chain->SetBranchAddress("grady", &grady);
	tree_chain->SetBranchAddress("isKaon", &isKaon);

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
	TH1D* t_prim_alth = new TH1D("t_prime_alt", "t_prime_alt", 1000, -5, 5);
	hists.push_back(t_prim_alth);
	TH1D* delta_phih = new TH1D("delta_phi", "delta_phi", 500, -7, 7);
	hists.push_back(delta_phih);
	TH1D* trig_maskh = new TH1D("trigger_mask", "trigger_mask", 15, 0, 15);
	hists.push_back(trig_maskh);

	for(unsigned int i = 0; i < tree_chain->GetEntries(); ++i) {

		tree_chain->GetEntry(i);
		stats->Fill("All events", 1);

/*std::cout<<"-------------"<<std::endl;
std::cout<<trigMask<<std::endl;
std::cout<<(trigMask&0x1)<<std::endl;
std::cout<<"-------------"<<std::endl;*/
		trig_maskh->Fill(trigMask);
		if(!(trigMask&0x1)) {
			continue;
		}
		stats->Fill("Trigger Mask = 1", 1);

		vtx_zh->Fill(vtx_z);
		if((vtx_z > -28.4) || (vtx_z < -68.4)) {
			continue;
		}
		stats->Fill("Vertex z in ]-28.4,-68.4[", 1);

		vtx_pos->Fill(vtx_x, vtx_y);
		if(std::pow(vtx_x, 2) + std::pow(vtx_y, 2) > 3.0625) {
			continue;
		}
		stats->Fill("Vertex.R() < 1.75", 1);

		rpd_mult->Fill(nbrRPDTracks);
		if(nbrRPDTracks != 1) {
			continue;
		}
		stats->Fill("1 RPD track", 1);

		proton.SetPxPyPzE(rpd_Px, rpd_Py, rpd_Pz, rpd_E);
		rpd_Pxh->Fill(proton.M());
		if(proton.M() < 0.2) {
			continue;
		}
		stats->Fill("Proton mass > 0.2", 1);

		if(isKaon != 0) {
			continue;
		}
		stats->Fill("isKaon = 0", 1);

		p1.SetXYZM(px1, py1, pz1, PION_MASS);
		p2.SetXYZM(px2, py2, pz2, PION_MASS);
		p3.SetXYZM(px3, py3, pz3, PION_MASS);
		p4.SetXYZM(px4, py4, pz4, PION_MASS);
		p5.SetXYZM(px5, py5, pz5, PION_MASS);

		pTot = p1+p2+p3+p4+p5;

		// Get beam 4-mom
		TVector3 p3_beam(gradx, grady, 1.);
		TLorentzVector pBeam = hlib::get_beam_energy(p3_beam, pTot);

		// Get the t's
		double t = std::fabs((pBeam - pTot).Mag2());
		// fhaas' method
		double t_min = std::fabs((std::pow(pTot.M2() - pBeam.M2(), 2)) / (4. * (pBeam.Vect()).Mag2()));
		double t_prim = t - t_min;
		// from suchung paper, seems to be wrong
		double t_min_alt = (pBeam - pTot).Mag2() - std::pow((pBeam.E() - pTot.E()), 2);
		double t_prim_alt = t - t_min_alt;

		t_primh->Fill(t_prim);
		t_prim_alth->Fill(t_prim_alt);

		if(t_prim < 0.1) {
			continue;
		}
		stats->Fill("T-prime > 0.1", 1);

		mom_5pi_raw->Fill(pTot.Energy());

		// Planarity cut
		double delta_phi, res;
		hlib::get_RPD_delta_phi_res(pBeam, proton, pTot, delta_phi, res);

		delta_phih->Fill(delta_phi);
		if(std::fabs(delta_phi) > res) {
			continue;
		}
		stats->Fill("RPD planarity cut", 1);

		mom_5pi->Fill(pTot.Energy());
		calc_beam_E->Fill(pBeam.E());
		if(std::fabs(pTot.Energy()-191.) > 3.28) {
			continue;
		}
		stats->Fill("Exclusivity 191+-3.28GeV", 1);

		mass_5pi->Fill(pTot.M());

		out_tree->Fill();

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
