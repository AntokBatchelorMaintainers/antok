
#include<iostream>
#include<sstream>
#include<string>

#include<TApplication.h>
#include<TCanvas.h>
#include<TFile.h>
#include<TTree.h>
#include<TH1D.h>
#include<TH2D.h>
#include<TLorentzVector.h>

#include<basic_calcs.h>
#include<constants.h>
#include<initializer.h>

void gen_kin(char* infile_name = 0, char* outfile_name = 0, std::string configfilename="../config/default.yaml") {

	antok::Initializer* initializer = antok::Initializer::instance();
	if(not initializer->readConfigFile(configfilename)) {
		std::cerr<<"Could not open config file. Aborting..."<<std::endl;
		exit(1);
	}

	const double& CHARGED_KAON_MASS = antok::Constants::chargedKaonMass();
	const double& PION_MASS = antok::Constants::chargedPionMass();

	new TApplication("app", 0, 0);

	TFile* infile;
	if(infile_name == 0) {
		infile = TFile::Open("/afs/cern.ch/user/k/kbicker/scratch0/filtered_run1/files_H_2008_10.root");
	} else {
		infile = TFile::Open(infile_name);
	}
	if(infile == 0) {
		std::cerr<<"Could not open input file. Aborting..."<<std::endl;
		return;
	}

	TTree* intree = (TTree*)infile->Get("Standard Event Selection/USR55");
	if(intree == 0) {
		std::cerr<<"Could not get input TTree. Aborting..."<<std::endl;
		return;
	}

	TFile* outfile;
	if(outfile_name == 0) {
		outfile = TFile::Open("outfile.root", "RECREATE");
	} else {
		outfile = TFile::Open(outfile_name, "NEW");
	}
	if(outfile == 0) {
		return;
	}

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
	std::vector<TLorentzVector*> allMomenta;
	allMomenta.push_back(&p1);
	allMomenta.push_back(&p2);
	allMomenta.push_back(&p3);
	allMomenta.push_back(&p4);
	allMomenta.push_back(&p5);

	double beamx, beamy, beamz;

	TLorentzVector k1;
	TLorentzVector k2;
	TLorentzVector k3;
	TLorentzVector k4;
	TLorentzVector k5;

	TLorentzVector pTot;

	double gradx, grady;

	const bool GENERATED_MONTECARLO = false;

	if(not GENERATED_MONTECARLO) {
		intree->SetBranchAddress("Mom_x1", &px1);
		intree->SetBranchAddress("Mom_y1", &py1);
		intree->SetBranchAddress("Mom_z1", &pz1);
		intree->SetBranchAddress("Mom_x2", &px2);
		intree->SetBranchAddress("Mom_y2", &py2);
		intree->SetBranchAddress("Mom_z2", &pz2);
		intree->SetBranchAddress("Mom_x3", &px3);
		intree->SetBranchAddress("Mom_y3", &py3);
		intree->SetBranchAddress("Mom_z3", &pz3);
		intree->SetBranchAddress("Mom_x4", &px4);
		intree->SetBranchAddress("Mom_y4", &py4);
		intree->SetBranchAddress("Mom_z4", &pz4);
		intree->SetBranchAddress("Mom_x5", &px5);
		intree->SetBranchAddress("Mom_y5", &py5);
		intree->SetBranchAddress("Mom_z5", &pz5);
		intree->SetBranchAddress("gradx", &gradx);
		intree->SetBranchAddress("grady", &grady);
	} else {
		intree->SetBranchAddress("Mom_MCTruth_x1", &px1);
		intree->SetBranchAddress("Mom_MCTruth_y1", &py1);
		intree->SetBranchAddress("Mom_MCTruth_z1", &pz1);
		intree->SetBranchAddress("Mom_MCTruth_x2", &px2);
		intree->SetBranchAddress("Mom_MCTruth_y2", &py2);
		intree->SetBranchAddress("Mom_MCTruth_z2", &pz2);
		intree->SetBranchAddress("Mom_MCTruth_x3", &px3);
		intree->SetBranchAddress("Mom_MCTruth_y3", &py3);
		intree->SetBranchAddress("Mom_MCTruth_z3", &pz3);
		intree->SetBranchAddress("Mom_MCTruth_x4", &px4);
		intree->SetBranchAddress("Mom_MCTruth_y4", &py4);
		intree->SetBranchAddress("Mom_MCTruth_z4", &pz4);
		intree->SetBranchAddress("Mom_MCTruth_x5", &px5);
		intree->SetBranchAddress("Mom_MCTruth_y5", &py5);
		intree->SetBranchAddress("Mom_MCTruth_z5", &pz5);
		intree->SetBranchAddress("Mom_x0_MCTruth", &beamx);
		intree->SetBranchAddress("Mom_y0_MCTruth", &beamy);
		intree->SetBranchAddress("Mom_z0_MCTruth", &beamz);
	}

	std::vector<TH1*> hists;
	hists.push_back(new TH1D("5_pi_mass", "5 Pion Mass", 1000, 0., 7.));
	hists.push_back(new TH1D("4_pi_mass", "4 Pion Subsystem", 500, 0.5, 3.5));
	hists.push_back(new TH1D("3_pi_mass", "3 Pion Subsystem", 1500, 0.5, 2.5));
	hists.push_back(new TH2D("3_dalitz", "3 Pion Subsystem", 1000, 0., 2.5, 1000, 0., 2.5));
	hists.push_back(new TH2D("3_dalitz_a2", "3 Pion Subsystem in a2 region", 700, 0., 2.5, 700, 0., 2.5));
	hists.push_back(new TH1D("2_pi_mass", "2 Pion Subsystem", 1000, 0.2, 2));
	hists.push_back(new TH2D("2_pi_4_pi_dalitz", "4 Pion Subsystem", 600, 0.5, 3., 600, 0.2, 1.));

	std::vector<double> bounds;
	bounds.resize(10, 0.784);
	bounds.at(1) = 0.961;
	bounds.at(2) = 1.137;
	bounds.at(3) = 1.312;
	bounds.at(4) = 1.488;
	bounds.at(5) = 1.664;
	bounds.at(6) = 1.840;
	bounds.at(7) = 2.016;
	bounds.at(8) = 2.192;
	bounds.at(9) = 2.368;

	for(unsigned int i = 1; i < bounds.size(); ++i) {
		std::ostringstream strs;
		strs << bounds.at(i-1);
		std::string lb = strs.str();
		std::ostringstream strs2;
		strs2 << bounds.at(i);
		std::string ub = strs2.str();
		std::string name("2_pi_4_pi_dalitz_" + lb + "_" + ub);
		std::string title("2 Pion [" + lb + ", " + ub + "[");
		hists.push_back(new TH2D(name.c_str(), title.c_str(), 200, 0.3, 1.3, 200, 0.3, 1.4));
	}

	hists.push_back(new TH1D("3_pi-_rapidity", "3 negative Pion Rapidity", 2000, 0., 10.));
	hists.push_back(new TH1D("Max_pi-_rapidity", "negative Pion with max. Rapidity", 2000, 0., 10.));
	hists.push_back(new TH1D("4_pi_mass_---+", "4 Pion Subsystem (---+)", 500, 0.5, 3.5));
	hists.push_back(new TH2D("4_pi_mass_ag_rapidity", "4 Pion Mass against Rapidity", 250, 0.5, 3.5, 1000, 0., 10.));
	hists.push_back(new TH2D("5_pi_mass_ag_rapidity", "5 Pion Mass against Rapidity", 250, 0.5, 3.5, 1000, 0., 10.));

	hists.push_back(new TH2D("4_pi_mass_ag_tp", "4 Pion Mass against t prime", 250, 0.5, 3.5, 1000, 0.1, 1.));
	hists.push_back(new TH2D("5_pi_mass_ag_tp", "5 Pion Mass against t prime", 1000, 0., 8., 1000, 0., 1.1));

	hists.push_back(new TH1D("2_kaon_mass", "2 Kaon (assumption) Subsystem", 1000, 0.9, 2));
	hists.push_back(new TH1D("5_pi_mass_pwa_bins", "5 Pion Mass", 233, 0.01, 7.));

	hists.push_back(new TH1D("pi_minus_momentum", "Momentum of negative Pions", 1000, 0, 210));
	hists.push_back(new TH1D("pi_plus_momentum", "Momentum of positive Pions", 1000, 0, 210));

	hists.push_back(new TH2D("xF_fast_pion_ag_rapidity_gap", "x_F of fast Pion against Rapidity Gap", 1000, -10, 10, 1000, 0, 1));

	hists.push_back(new TH1D("fast_pi-_momentum", "Momentum of fastest negative Pion", 1000, 0, 210));

	for(unsigned int i = 0; i < intree->GetEntries(); ++i) {

		intree->GetEntry(i);

		p1.SetXYZM(px1, py1, pz1, PION_MASS);
		p2.SetXYZM(px2, py2, pz2, PION_MASS);
		p3.SetXYZM(px3, py3, pz3, PION_MASS);
		p4.SetXYZM(px4, py4, pz4, PION_MASS);
		p5.SetXYZM(px5, py5, pz5, PION_MASS);

		k1.SetXYZM(px1, py1, pz1, CHARGED_KAON_MASS);
		k2.SetXYZM(px2, py2, pz2, CHARGED_KAON_MASS);
		k3.SetXYZM(px3, py3, pz3, CHARGED_KAON_MASS);
		k4.SetXYZM(px4, py4, pz4, CHARGED_KAON_MASS);
		k5.SetXYZM(px5, py5, pz5, CHARGED_KAON_MASS);

		pTot = p1+p2+p3+p4+p5;

		TVector3 p3Beam;
		TLorentzVector pBeam;
		if(not GENERATED_MONTECARLO) {
			p3Beam.SetXYZ(gradx, grady, 1.);
			pBeam = antok::getBeamEnergy(p3Beam, pTot);
			p3Beam = pBeam.Vect();
		} else {
			p3Beam.SetXYZ(beamx, beamy, beamz);
			pBeam.SetXYZM(p3Beam.X(), p3Beam.Y(), p3Beam.Z(), PION_MASS);
		}

		double t = std::fabs((pBeam - pTot).Mag2());
		double tMin = std::fabs((std::pow(pTot.M2() - pBeam.M2(), 2)) / (4. * p3Beam.Mag2()));
		double tPrime = t - tMin;

		hists.at(0)->Fill(pTot.M());
		hists.at(1)->Fill((p1+p2+p4+p5).M());
		hists.at(1)->Fill((p1+p3+p4+p5).M());
		hists.at(1)->Fill((p2+p3+p4+p5).M());
		hists.at(2)->Fill((p1+p2+p5).M());
		hists.at(2)->Fill((p1+p3+p5).M());
		hists.at(2)->Fill((p2+p3+p5).M());
		hists.at(2)->Fill((p1+p2+p4).M());
		hists.at(2)->Fill((p1+p3+p4).M());
		hists.at(2)->Fill((p2+p3+p4).M());

		hists.at(3)->Fill((p1+p5).M(), (p2+p5).M());
		hists.at(3)->Fill((p1+p5).M(), (p3+p5).M());
		hists.at(3)->Fill((p2+p5).M(), (p3+p5).M());
		hists.at(3)->Fill((p1+p4).M(), (p2+p4).M());
		hists.at(3)->Fill((p1+p4).M(), (p3+p4).M());
		hists.at(3)->Fill((p2+p4).M(), (p3+p4).M());

		const double a2_mass = 1.318;
		if(std::fabs((p1+p2+p5).M()) - a2_mass < 0.107) {
			hists.at(4)->Fill((p1+p5).M(), (p2+p5).M());
		}
		if(std::fabs((p1+p3+p5).M()) - a2_mass < 0.107) {
			hists.at(4)->Fill((p1+p5).M(), (p3+p5).M());
		}
		if(std::fabs((p2+p3+p5).M()) - a2_mass < 0.107) {
			hists.at(4)->Fill((p2+p5).M(), (p3+p5).M());
		}
		if(std::fabs((p1+p2+p4).M()) - a2_mass < 0.107) {
			hists.at(4)->Fill((p1+p4).M(), (p2+p4).M());
		}
		if(std::fabs((p1+p3+p4).M()) - a2_mass < 0.107) {
			hists.at(4)->Fill((p1+p4).M(), (p3+p4).M());
		}
		if(std::fabs((p2+p3+p4).M()) - a2_mass < 0.107) {
			hists.at(4)->Fill((p2+p4).M(), (p3+p4).M());
		}

		hists.at(5)->Fill((p1+p5).M());
		hists.at(5)->Fill((p2+p5).M());
		hists.at(5)->Fill((p3+p5).M());
		hists.at(5)->Fill((p1+p4).M());
		hists.at(5)->Fill((p2+p4).M());
		hists.at(5)->Fill((p3+p4).M());

		hists.at(6)->Fill((p1+p2+p4+p5).M(), (p1+p4).M());
		hists.at(6)->Fill((p1+p2+p4+p5).M(), (p1+p5).M());
		hists.at(6)->Fill((p1+p2+p4+p5).M(), (p2+p4).M());
		hists.at(6)->Fill((p1+p2+p4+p5).M(), (p2+p5).M());
		hists.at(6)->Fill((p1+p3+p4+p5).M(), (p1+p4).M());
		hists.at(6)->Fill((p1+p3+p4+p5).M(), (p1+p5).M());
		hists.at(6)->Fill((p1+p3+p4+p5).M(), (p3+p4).M());
		hists.at(6)->Fill((p1+p3+p4+p5).M(), (p3+p5).M());
		hists.at(6)->Fill((p2+p3+p4+p5).M(), (p2+p4).M());
		hists.at(6)->Fill((p2+p3+p4+p5).M(), (p2+p5).M());
		hists.at(6)->Fill((p2+p3+p4+p5).M(), (p3+p4).M());
		hists.at(6)->Fill((p2+p3+p4+p5).M(), (p3+p5).M());

		for(unsigned int i = 1; i < bounds.size(); ++i) {
			if(((p1+p2+p4+p5).M() >= bounds.at(i-1)) && ((p1+p2+p4+p5).M() < bounds.at(i))) {
				hists.at(6+i)->Fill((p1+p4).M(), (p2+p5).M());
				hists.at(6+i)->Fill((p2+p4).M(), (p1+p5).M());
			}
			if(((p1+p3+p4+p5).M() >= bounds.at(i-1)) && ((p1+p3+p4+p5).M() < bounds.at(i))) {
				hists.at(6+i)->Fill((p1+p4).M(), (p3+p5).M());
				hists.at(6+i)->Fill((p3+p4).M(), (p1+p5).M());
			}
			if(((p2+p3+p4+p5).M() >= bounds.at(i-1)) && ((p2+p3+p4+p5).M() < bounds.at(i))) {
				hists.at(6+i)->Fill((p2+p4).M(), (p3+p5).M());
				hists.at(6+i)->Fill((p3+p4).M(), (p2+p5).M());
			}
		}

		hists.at(6 + bounds.size())->Fill(p1.Rapidity());
		hists.at(6 + bounds.size())->Fill(p2.Rapidity());
		hists.at(6 + bounds.size())->Fill(p3.Rapidity());

		TLorentzVector* p_max_rap = &p1;
		unsigned int p_max_rap_index = 1;
		if(p_max_rap->Rapidity() < p2.Rapidity()) {
			p_max_rap = &p2;
			p_max_rap_index = 2;
		}
		if(p_max_rap->Rapidity() < p3.Rapidity()) {
			p_max_rap = &p3;
			p_max_rap_index = 3;
		}
		TLorentzVector p_4pi_central(0., 0., 0., 0.);
		for(unsigned int j = 1; j <= 5; ++j) {
			if(j == p_max_rap_index) {
				continue;
			}
			p_4pi_central += *(allMomenta[j-1]);
		}

		hists.at(7 + bounds.size())->Fill(p_max_rap->Rapidity());
		hists.at(8 + bounds.size())->Fill((p1+p2+p3+p4).M());
		hists.at(8 + bounds.size())->Fill((p1+p2+p3+p5).M());

		hists.at(9 + bounds.size())->Fill((p1+p2+p4+p5).M(), p3.Rapidity());
		hists.at(9 + bounds.size())->Fill((p1+p3+p4+p5).M(), p2.Rapidity());
		hists.at(9 + bounds.size())->Fill((p2+p3+p4+p5).M(), p1.Rapidity());

		hists.at(10 + bounds.size())->Fill((p1+p2+p3+p4+p5).M(), p_max_rap->Rapidity());

		hists.at(11 + bounds.size())->Fill((p1+p2+p4+p5).M(), tPrime);
		hists.at(11 + bounds.size())->Fill((p1+p3+p4+p5).M(), tPrime);
		hists.at(11 + bounds.size())->Fill((p2+p3+p4+p5).M(), tPrime);

		hists.at(12 + bounds.size())->Fill((p1+p2+p3+p4+p5).M(), tPrime);

		hists.at(13 + bounds.size())->Fill((k1+k5).M());
		hists.at(13 + bounds.size())->Fill((k2+k5).M());
		hists.at(13 + bounds.size())->Fill((k3+k5).M());
		hists.at(13 + bounds.size())->Fill((k1+k4).M());
		hists.at(13 + bounds.size())->Fill((k2+k4).M());
		hists.at(13 + bounds.size())->Fill((k3+k4).M());

		hists.at(14 + bounds.size())->Fill(pTot.M());

		hists.at(15 + bounds.size())->Fill(p1.Vect().Mag());
		hists.at(15 + bounds.size())->Fill(p2.Vect().Mag());
		hists.at(15 + bounds.size())->Fill(p3.Vect().Mag());

		hists.at(16 + bounds.size())->Fill(p4.Vect().Mag());
		hists.at(16 + bounds.size())->Fill(p5.Vect().Mag());

		double centerOfMassEnergy;
		TVector3 boostToCenterOfMassSystem;
		antok::getBoostToCenterOfMassSystem(pBeam, centerOfMassEnergy, boostToCenterOfMassSystem);
		{
			TLorentzVector p_max_rap_CM = *p_max_rap;
			p_max_rap_CM.Boost(-boostToCenterOfMassSystem);
			double xF_maxRap = (2 * p_max_rap_CM.Pz()) / centerOfMassEnergy;
			hists.at(17 + bounds.size())->Fill(p_max_rap->Rapidity() - p_4pi_central.Rapidity(), xF_maxRap);
		}

		hists.at(18 + bounds.size())->Fill(p_max_rap->Vect().Mag());

	}

	outfile->cd();
	outfile->mkdir("Canvases");
	TCanvas* c2 = new TCanvas("2_pions_combined", "2_pions_combined");
	unsigned int needed = bounds.size() - 1;
	double root = std::sqrt(needed);
	int spl1 = (int)root;
	double spl2d = needed / spl1;
	int spl2;
	if((spl2d - std::floor(spl2d)) > 0.) {
		spl2 = (int)std::floor(spl2d + 1.);
	} else {
		spl2 = (int)spl2d;
	}
	c2->Divide(spl2, spl1);

	for(unsigned int i = 0; i < hists.size(); ++i) {
		hists.at(i)->Write();
		std::string name(std::string(hists.at(i)->GetName()) + "_c");
		outfile->Cd("Canvases");
		TCanvas* c;
		c = new TCanvas(name.c_str(), hists.at(i)->GetTitle());
		hists.at(i)->Draw();
		if(dynamic_cast<TH2D*>(hists.at(i)) != 0) {
			hists.at(i)->SetDrawOption("colz");
		}
		c->Write();
		if((i >= 7) && (i < (6+bounds.size()))) {
			c2->cd(i-6);
			hists.at(i)->Draw();
			hists.at(i)->SetDrawOption("colz");
		}
		outfile->cd();
	}
	outfile->Cd("Canvases");
	c2->Write();
	outfile->Close();

}

int main(int argc, char* argv[]) {
	if(argc == 1) {
		gen_kin();
	} else if (argc == 3) {
		gen_kin(argv[1], argv[2]);
	} else if (argc == 4) {
		gen_kin(argv[1], argv[2], argv[3]);
	} else {
		std::cerr<<"Wrong number of arguments, is "<<argc<<", should be in [0, 2, 3]."<<std::endl;
	}
}
