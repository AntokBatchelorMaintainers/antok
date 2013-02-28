
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
	Double_t xMin[5] = {-2, -2, -0.8, -0.8, 187};
	Double_t xMax[5] = {2, 2, 0.8, 0.8, 196};

	THnSparseD* hist = new THnSparseD("beamHist", "beamHist", 5, bins, xMin, xMax);
	TH1D* vtxX = new TH1D("vtxX", "vtxX", 10000, -2, 2);
	TH1D* vtxY = new TH1D("vtxY", "vtxY", 10000, -2, 2);
	TH1D* momX = new TH1D("momX", "momX", 50000, -0.8, 0.8);
	TH1D* momY = new TH1D("momY", "momY", 50000, -0.8, 0.8);
	TH1D* momZ = new TH1D("momZ", "momZ", 10000, 187, 196);

	TH1D* edgeLength = new TH1D("edgeLength", "edgeLength", 501, -0.5, 500.5);

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

	double bx, by, bz;
	TTree* outTree = new TTree("beamTree", "beamTree");
	outTree->Branch("vertex_x_position", &primVx, "vertex_x_position/D");
	outTree->Branch("vertex_y_position", &primVy, "vertex_y_position/D");
	outTree->Branch("beam_momentum_x", &bx,"beam_momentum_x/D");
	outTree->Branch("beam_momentum_y", &by,"beam_momentum_y/D");
	outTree->Branch("beam_momentum_z", &bz,"beam_momentum_z/D");

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
		bx = beam.X();
		by = beam.Y();
		bz = beam.Z();

		double values[5] = {primVx, primVy, bx, by, bz};

		vtxX->Fill(primVx);
		vtxY->Fill(primVy);
		momX->Fill(bx);
		momY->Fill(by);
		momZ->Fill(bz);

		hist->Fill(values);
		outTree->Fill();

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

	double primVx_sigma, primVy_sigma, bx_sigma, by_sigma, bz_sigma;
	int cubeSize = 1;
	double binContent = 0;
	TTree* sigmaTree = new TTree("sigmaTree", "sigmaTree");
	sigmaTree->Branch("vertex_x_position_sigma", &primVx_sigma, "vertex_x_position_sigma/D");
	sigmaTree->Branch("vertex_y_position_sigma", &primVy_sigma, "vertex_y_position_sigma/D");
	sigmaTree->Branch("beam_momentum_x_sigma", &bx_sigma,"beam_momentum_x_sigma/D");
	sigmaTree->Branch("beam_momentum_y_sigma", &by_sigma,"beam_momentum_y_sigma/D");
	sigmaTree->Branch("beam_momentum_z_sigma", &bz_sigma,"beam_momentum_z_sigma/D");
	sigmaTree->Branch("cube_edge_length", &cubeSize, "cube_edge_length/I");
	sigmaTree->Branch("bin_content", &binContent, "bin_content/D");

	const double MIN_ENTRIES = 10.;

	for(unsigned int i = 0; i < entries; ++i) {

		outTree->GetEntry(i);
		double values[5] = {primVx, primVy, bx, by, bz};
		int bin = hist->GetBin(values, false);
		int location[5] = {0, 0, 0, 0, 0};
		binContent = hist->GetBinContent(bin, location);
		cubeSize = 1;
		if(binContent < MIN_ENTRIES) {
			std::pair<int, int> verXRange(location[0], location[0]);
			std::pair<int, int> verYRange(location[1], location[1]);
			std::pair<int, int> momXRange(location[2], location[2]);
			std::pair<int, int> momYRange(location[3], location[3]);
			std::pair<int, int> momZRange(location[4], location[4]);

			while(binContent < MIN_ENTRIES) {

				verXRange.first -= 1; verXRange.second += 1;
				verYRange.first -= 1; verYRange.second += 1;
				momXRange.first -= 1; momXRange.second += 1;
				momYRange.first -= 1; momYRange.second += 1;
				momZRange.first -= 1; momZRange.second += 1;
				int step = verXRange.second - verXRange.first;
				cubeSize = step + 1;

				for(int dim = 0; dim < 5; ++dim) {
					for(int l = (dim > 0) ? verXRange.first + 1 : verXRange.first; l <= ((dim > 0) ? verXRange.second - 1 : verXRange.second); l += (dim == 0) ? step : 1) {
						for(int m = (dim > 1) ? verYRange.first + 1 : verYRange.first; m <= ((dim > 1) ? verYRange.second - 1 : verYRange.second); m += (dim == 1) ? step : 1) {
							for(int n = (dim > 2) ? momXRange.first + 1 : momXRange.first; n <= ((dim > 2) ? momXRange.second - 1 : momXRange.second); n += (dim == 2) ? step : 1) {
								for(int o = (dim > 3) ? momYRange.first + 1 : momYRange.first; o <= ((dim > 3) ? momYRange.second - 1 : momYRange.second); o += (dim == 3) ? step : 1) {
									for(int p = momZRange.first; p <= momZRange.second; p += (dim == 4) ? step : 1) {
										int coord[5] = {l, m, n, o, p};
										int binnumber = hist->GetBin(coord, false);
										if(binnumber > 0) {
											binContent += hist->GetBinContent(bin);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		edgeLength->Fill(cubeSize);
		primVx_sigma = (cubeSize * (hist->GetAxis(0)->GetBinWidth(location[0]))) / binContent;
		primVy_sigma = (cubeSize * (hist->GetAxis(1)->GetBinWidth(location[1]))) / binContent;
		bx_sigma = (cubeSize * (hist->GetAxis(2)->GetBinWidth(location[2]))) / binContent;
		by_sigma = (cubeSize * (hist->GetAxis(3)->GetBinWidth(location[3]))) / binContent;
		bz_sigma = (cubeSize * (hist->GetAxis(4)->GetBinWidth(location[4]))) / binContent;
		sigmaTree->Fill();

		if(not (i % 100)) {
			std::cout<<i<<std::endl;
		}

	}

	outTree->AddFriend(sigmaTree);
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
