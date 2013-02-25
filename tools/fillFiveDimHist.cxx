
#include<iostream>
#include<string>
#include<time.h>

#include<TFile.h>
#include<TH2D.h>
#include<THnSparse.h>
#include<TLorentzVector.h>
#include<TTree.h>

#include<basic_calcs.h>
#include<constants.h>
#include<initializer.h>

void fillFiveDimHist(std::string inFileName) {

	std::string configfilename = "../../antok/config/monte_carlo_run2.yaml";
	antok::Initializer* initializer = antok::Initializer::instance();
	if(not initializer->readConfigFile(configfilename)) {
		std::cerr<<"Could not open config file. Aborting..."<<std::endl;
		exit(1);
	}
	const double& PION_MASS = antok::Constants::charged_pion_mass();

	TFile* inFile = TFile::Open(inFileName.c_str(), "READ");
	TTree* inTree = (TTree*)inFile->Get("Standard Event Selection/USR55");
	TFile* outFile = TFile::Open("bla.root", "RECREATE");
	outFile->cd();

	double px1, py1, pz1;
	double px2, py2, pz2;
	double px3, py3, pz3;
	double px4, py4, pz4;
	double px5, py5, pz5;
	double gradx, grady;
	double primVx, primVy;

	inTree->SetBranchAddress("Mom_x1", &px1);
	inTree->SetBranchAddress("Mom_y1", &py1);
	inTree->SetBranchAddress("Mom_z1", &pz1);
	inTree->SetBranchAddress("Mom_x2", &px2);
	inTree->SetBranchAddress("Mom_y2", &py2);
	inTree->SetBranchAddress("Mom_z2", &pz2);
	inTree->SetBranchAddress("Mom_x3", &px3);
	inTree->SetBranchAddress("Mom_y3", &py3);
	inTree->SetBranchAddress("Mom_z3", &pz3);
	inTree->SetBranchAddress("Mom_x4", &px4);
	inTree->SetBranchAddress("Mom_y4", &py4);
	inTree->SetBranchAddress("Mom_z4", &pz4);
	inTree->SetBranchAddress("Mom_x5", &px5);
	inTree->SetBranchAddress("Mom_y5", &py5);
	inTree->SetBranchAddress("Mom_z5", &pz5);
	inTree->SetBranchAddress("gradx", &gradx);
	inTree->SetBranchAddress("grady", &grady);
	inTree->SetBranchAddress("X_primV", &primVx);
	inTree->SetBranchAddress("Y_primV", &primVy);

	std::vector<TLorentzVector> particles;
	particles.resize(6);
	unsigned int entries = inTree->GetEntries();

	double bx, by, bz;
	TTree* outTree = new TTree("beamTree", "beamTree");
	outTree->Branch("vertex_x_position", &primVx, "vertex_x_position/D");
	outTree->Branch("vertex_y_position", &primVy, "vertex_y_position/D");
	outTree->Branch("beam_momentum_x", &bx,"beam_momentum_x/D");
	outTree->Branch("beam_momentum_y", &by,"beam_momentum_y/D");
	outTree->Branch("beam_momentum_z", &bz,"beam_momentum_z/D");

	for(unsigned int i = 0; i < entries; ++i) {

		inTree->GetEntry(i);

		particles[1].SetXYZM(px1, py1, pz1, PION_MASS);
		particles[2].SetXYZM(px2, py2, pz2, PION_MASS);
		particles[3].SetXYZM(px3, py3, pz3, PION_MASS);
		particles[4].SetXYZM(px4, py4, pz4, PION_MASS);
		particles[5].SetXYZM(px5, py5, pz5, PION_MASS);

		TLorentzVector pSum;
		for(unsigned int j = 1; j < 6; ++j) {
			pSum += particles.at(j);
		}
		particles[0] = antok::get_beam_energy(TVector3(gradx, grady, 1.), pSum);
		TVector3 beam = particles[0].Vect();
		bx = beam.X();
		by = beam.Y();
		bz = beam.Z();

		outTree->Fill();

		if(i % 100000 == 0) {
			std::cout<<"Filling beamTree: Entry "<<i<<" of "<<entries<<std::endl;
		}

	}

	double primVx_sigma, primVy_sigma, bx_sigma, by_sigma, bz_sigma;
	TTree* sigmaTree = new TTree("sigmaTree", "sigmaTree");
	outTree->Branch("vertex_x_position_sigma", &primVx_sigma, "vertex_x_position_sigma/D");
	outTree->Branch("vertex_y_position_sigma", &primVy_sigma, "vertex_y_position_sigma/D");
	outTree->Branch("beam_momentum_x_sigma", &bx_sigma,"beam_momentum_x_sigma/D");
	outTree->Branch("beam_momentum_y_sigma", &by_sigma,"beam_momentum_y_sigma/D");
	outTree->Branch("beam_momentum_z_sigma", &bz_sigma,"beam_momentum_z_sigma/D");

	for(unsigned int i = 0; i < entries; ++i) {

		outTree->GetEntry(i);
		double refVX = primVx;
		double refVY = primVy;
		double refbx = bx;
		double refby = by;
		double refbz = bz;

		time_t time_first = time(NULL);

		const unsigned int NVALUES = 10;

		std::vector<double> dists;
		std::vector<int> closest_events;

		for(unsigned int j = 0; j < entries; ++j) {

			outTree->GetEntry(j);

			double dist = std::sqrt( (primVx-refVX)*(primVx-refVX) +
			                         (primVy-refVY)*(primVy-refVY) +
			                         (bx-refbx)*(bx-refbx) +
			                         (by-refby)*(by-refby) +
			                         (bz-refbz)*(bz-refbz) );

			if(dists.size() < NVALUES) {
				dists.push_back(dist);
				closest_events.push_back(j);
			} else {
				unsigned int max_k = 0;
				double max = 0;
				bool set = false;
				for(unsigned int k = 0; k < NVALUES; ++k) {
					if(not set && dist < dists[k]) {
						set = true;
						max = dists[k];
						max_k = k;
					}
					if(set) {
						if(dists[k] > max) {
							max = dists[k];
							max_k = k;
						}
					}
				}
				if(set) {
					dists[max_k] = dist;
					closest_events[max_k] = j;
				}
			}
		}

		primVx_sigma = 0;
		primVy_sigma = 0;
		bx_sigma = 0;
		by_sigma = 0;
		bz_sigma = 0;

		for(unsigned int j = 0; j < NVALUES; ++j) {
			int index = closest_events[j];
			outTree->GetEntry(index);
			primVx_sigma += std::fabs(primVx - refVX);
			primVy_sigma += std::fabs(primVy - refVY);
			bx_sigma += std::fabs(bx - refbx);
			by_sigma += std::fabs(by - refby);
			bz_sigma += std::fabs(bz - refbz);
		}
		primVx_sigma /= NVALUES;
		primVy_sigma /= NVALUES;
		bx_sigma /= NVALUES;
		by_sigma /= NVALUES;
		bz_sigma /= NVALUES;
		sigmaTree->Fill();

		time_t time_last = time(NULL);
		double seconds = difftime(time_last, time_first);

		std::cout<<"Filling sigmaTree: Entry "<<i<<" of "<<entries<<" (took "<<seconds<<" seconds)"<<std::endl;

	}

	outTree->AddFriend(sigmaTree);

	outFile->Write();
	outFile->Close();

}

int main(int argc, char* argv[]) {
	if(argc == 2) {
		fillFiveDimHist(argv[1]);
	} else {
		std::cerr<<"Wrong number of arguments, is "<<argc<<", should be 1."<<std::endl;
	}
}
