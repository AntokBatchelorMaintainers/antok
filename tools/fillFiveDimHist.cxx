
#include<iostream>
#include<string>

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

	Int_t bins[5] = {50, 50, 400, 400, 50};
	double_t xMin[5] = {-2, -2, -0.8, -0.8, 187};
	double_t xMax[5] = {2, 2, 0.8, 0.8, 196};

	THnSparseD* hist = new THnSparseD("beamHist", "beamHist", 5, bins, xMin, xMax);
	TH1D* vtxX = new TH1D("vtxX", "vtxX", 10000, -2, 2);
	TH1D* vtxY = new TH1D("vtxY", "vtxY", 10000, -2, 2);
	TH1D* momX = new TH1D("momX", "momX", 50000, -0.8, 0.8);
	TH1D* momY = new TH1D("momY", "momY", 50000, -0.8, 0.8);
	TH1D* momZ = new TH1D("momZ", "momZ", 10000, 187, 196);

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
		double bx = beam.X();
		double by = beam.Y();
		double bz = beam.Z();

		double values[5] = {primVx, primVy, bx, by, bz};

		vtxX->Fill(primVx);
		vtxY->Fill(primVy);
		momX->Fill(bx);
		momY->Fill(by);
		momZ->Fill(bz);

		hist->Fill(values);

		if(i % 100000 == 0) {
			std::cout<<"Entry "<<i<<" of "<<entries<<std::endl;
		}

	}

	hist->Projection(1, 0)->Write();
	hist->Projection(3, 2)->Write();
	hist->Projection(2, 1)->Write();
	hist->Projection(2, 0)->Write();
	hist->Projection(3, 1)->Write();
	hist->Projection(3, 0)->Write();
	hist->Projection(4, 0)->Write();
	hist->Projection(4, 1)->Write();
	hist->Projection(4, 2)->Write();
	hist->Projection(4, 3)->Write();

	std::cout<<"Histogram has "<<hist->GetNbins()<<" filled bins."<<std::endl;

	hist->Write();
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
