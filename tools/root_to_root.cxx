
#include<cmath>
#include<iostream>
#include<sstream>
#include<sys/stat.h>
#include<sys/types.h>

#include<TClonesArray.h>
#include<TFile.h>
#include<TLorentzVector.h>
#include<TObjString.h>
#include<TTree.h>
#include<TVector3.h>

#include<basic_calcs.h>
#include<constants.h>
#include<initializer.h>

void convert_root_to_txt(char* infile_name,
                         char* outfile_name,
                         std::string configfilename,
                         std::string inputType)
{

	if(not (inputType == "data" or inputType == "gen" or inputType == "acc")) {
		std::cerr<<"Invalid input type '"<<inputType<<"', must be in {'data', 'gen', 'acc'}"<<std::endl;
		return;
	}

	antok::Initializer* initializer = antok::Initializer::instance();
	if(not initializer->readConfigFile(configfilename)) {
		std::cerr<<"Could not open config file. Aborting..."<<std::endl;
		exit(1);
	}

	const double& PION_MASS = antok::Constants::chargedPionMass();

	const std::string PHASE_SPACE_FILENAME_POSTFIX = "genPS";
	const std::string RECO_FILENAME_POSTFIX = "accPS";

	double bin_width = 0.03;		//GeV
	double mass_range_lb = 1.3;		//GeV
	double mass_range_ub = 4.;		//GeV

	if((int)(floor(mass_range_ub*1000. - mass_range_lb*1000.) + 0.5) % (int)(floor(bin_width*1000. + 0.5))) {
		std::cout<<"Mass bins don't fit into range"<<std::endl;
		return;
	}

	// Get all the bounds.
	std::vector<double> bounds;
	bounds.resize(((int)((mass_range_ub - mass_range_lb) / bin_width)) + 1);
	for(unsigned int i = 0; i < bounds.size(); ++i) {
		bounds.at(i) = mass_range_lb + i*bin_width;
	}

	// Make dirs and open output files.
	std::vector<TFile*> tfiles;
	tfiles.resize(bounds.size()-1, 0);
	for(unsigned int i = 0; i < tfiles.size(); ++i) {
		std::ostringstream strs;
		strs<<outfile_name<<"/"<<(int)(floor(bounds.at(i)*1000. + 0.5))<<"."<<(int)(floor(bounds.at(i+1)*1000. + 0.5));
		mkdir(strs.str().c_str(), S_IRWXU | S_IRWXG);
		strs<<"/"<<(int)(floor(bounds.at(i)*1000. + 0.5))<<"."<<(int)(floor(bounds.at(i+1)*1000. + 0.5));
		if(inputType == "gen") {
			strs<<"."<<PHASE_SPACE_FILENAME_POSTFIX;
		} else if(inputType == "acc") {
			strs<<"."<<RECO_FILENAME_POSTFIX;
		}
		strs<<".root";
		tfiles.at(i) = TFile::Open(strs.str().c_str(), "NEW");
		if(tfiles.at(i) == nullptr) {
			std::cout<<"Error opening file for writing."<<std::endl;
			return;
		}
	}

	// Open input file and do all the tree stuff
	TFile* infile = TFile::Open(infile_name, "READ");
	if(infile == nullptr) {
		return;
	}
	std::string inFileName = "Standard Event Selection/USR55";
	if(inputType == "gen") {
		inFileName = "kbicker_5pic/USR55";
	}
	TTree* tree = (TTree*)infile->Get(inFileName.c_str());
	if(tree == nullptr) {
		std::cout<<"Error opening in-TTree."<<std::endl;
		return;
	}
	double px0, py0, pz0;
	double px1, py1, pz1;
	double px2, py2, pz2;
	double px3, py3, pz3;
	double px4, py4, pz4;
	double px5, py5, pz5;
	double gradx, grady;

	if(inputType == "data") {
		tree->SetBranchAddress("Mom_x1", &px1);
		tree->SetBranchAddress("Mom_y1", &py1);
		tree->SetBranchAddress("Mom_z1", &pz1);
		tree->SetBranchAddress("Mom_x2", &px2);
		tree->SetBranchAddress("Mom_y2", &py2);
		tree->SetBranchAddress("Mom_z2", &pz2);
		tree->SetBranchAddress("Mom_x3", &px3);
		tree->SetBranchAddress("Mom_y3", &py3);
		tree->SetBranchAddress("Mom_z3", &pz3);
		tree->SetBranchAddress("Mom_x4", &px4);
		tree->SetBranchAddress("Mom_y4", &py4);
		tree->SetBranchAddress("Mom_z4", &pz4);
		tree->SetBranchAddress("Mom_x5", &px5);
		tree->SetBranchAddress("Mom_y5", &py5);
		tree->SetBranchAddress("Mom_z5", &pz5);
		tree->SetBranchAddress("gradx", &gradx);
		tree->SetBranchAddress("grady", &grady);
	} else if (inputType == "gen" or inputType == "acc") {
		tree->SetBranchAddress("Mom_MCTruth_x1", &px1);
		tree->SetBranchAddress("Mom_MCTruth_y1", &py1);
		tree->SetBranchAddress("Mom_MCTruth_z1", &pz1);
		tree->SetBranchAddress("Mom_MCTruth_x2", &px2);
		tree->SetBranchAddress("Mom_MCTruth_y2", &py2);
		tree->SetBranchAddress("Mom_MCTruth_z2", &pz2);
		tree->SetBranchAddress("Mom_MCTruth_x3", &px3);
		tree->SetBranchAddress("Mom_MCTruth_y3", &py3);
		tree->SetBranchAddress("Mom_MCTruth_z3", &pz3);
		tree->SetBranchAddress("Mom_MCTruth_x4", &px4);
		tree->SetBranchAddress("Mom_MCTruth_y4", &py4);
		tree->SetBranchAddress("Mom_MCTruth_z4", &pz4);
		tree->SetBranchAddress("Mom_MCTruth_x5", &px5);
		tree->SetBranchAddress("Mom_MCTruth_y5", &py5);
		tree->SetBranchAddress("Mom_MCTruth_z5", &pz5);
		tree->SetBranchAddress("Mom_x0_MCTruth", &px0);
		tree->SetBranchAddress("Mom_y0_MCTruth", &py0);
		tree->SetBranchAddress("Mom_z0_MCTruth", &pz0);
	}

	std::vector<TLorentzVector> particles;
	particles.resize(6);

	TClonesArray prodMomName ("TObjString", 1);
	TClonesArray decayMomName("TObjString", 5);

	new (prodMomName  [0]) TObjString("pi-");
	new (decayMomName [0]) TObjString("pi-");
	new (decayMomName [1]) TObjString("pi-");
	new (decayMomName [2]) TObjString("pi-");
	new (decayMomName [3]) TObjString("pi+");
	new (decayMomName [4]) TObjString("pi+");

	TClonesArray* prodMom = new TClonesArray("TVector3");
	TClonesArray* decayMom = new TClonesArray("TVector3");

	std::vector<TTree*> trees;
	trees.resize(tfiles.size(), 0);

	const int splitLevel = 99;
	const int buffsize = 256000;

	// Create all the trees with their branches.
	for(unsigned int i = 0; i < tfiles.size(); ++i) {
		tfiles.at(i)->cd();
		trees.at(i) = new TTree("rootPwaEvtTree", "rootPwaEvtTree");
		trees.at(i)->Branch("prodKinMomenta", "TClonesArray", &prodMom, buffsize, splitLevel);
		trees.at(i)->Branch("decayKinMomenta", "TClonesArray", &decayMom, buffsize, splitLevel);
		prodMomName.Write("prodKinParticles", TObject::kSingleKey);
		decayMomName.Write("decayKinParticles", TObject::kSingleKey);
	}

	// Loop over events.
	for(unsigned int i = 0; i < tree->GetEntries(); ++i) {

		// Get stuff from the tree
		tree->GetEntry(i);

		particles.at(1).SetXYZM(px1, py1, pz1, PION_MASS);
		particles.at(2).SetXYZM(px2, py2, pz2, PION_MASS);
		particles.at(3).SetXYZM(px3, py3, pz3, PION_MASS);
		particles.at(4).SetXYZM(px4, py4, pz4, PION_MASS);
		particles.at(5).SetXYZM(px5, py5, pz5, PION_MASS);

		TLorentzVector pSum;
		for(unsigned int i = 1; i < 6; ++i) {
			pSum += particles.at(i);
		}
		double mass = pSum.M();
		if(inputType == "data") {
			particles.at(0) = antok::getBeamEnergy(TVector3(gradx, grady, 1.), pSum);
		} else {
			particles.at(0).SetXYZM(px0, py0, pz0, PION_MASS);
		}

		// If out of bounds, go to next event
		if(mass < bounds.at(0) || mass > bounds.at(bounds.size()-1)) {
			continue;
		}

		// Fill the TClonesArrays
		new ((*prodMom)[0]) TVector3(particles.at(0).Vect());
		for(unsigned int i = 1; i < 6; ++i) {
			new ((*decayMom)[i-1]) TVector3(particles.at(i).Vect());
		}

		// Fill the corresponding tree
		for(unsigned int i = 0; i < tfiles.size(); ++i) {
			if((mass > bounds.at(i)) && (mass < bounds.at(i+1))) {
				trees.at(i)->Fill();
			}
		}

	} // End loop over events

	for(unsigned int i = 0; i < tfiles.size(); ++i) {
		std::cout<<"Mass bin "<<i<<" has "<<trees.at(i)->GetEntries()<<" events."<<std::endl;
		tfiles.at(i)->cd();
		trees.at(i)->Write();
		tfiles.at(i)->Close();
	}

}

int main(int argc, char* argv[]) {
	if(argc == 5) {
		convert_root_to_txt(argv[1], argv[2], argv[3], argv[4]);
	} else {
		std::cerr<<"Wrong number of arguments, is "<<argc<<", should be 4."<<std::endl;
	}
}
