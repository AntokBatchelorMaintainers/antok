
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

void convert_root_to_txt(char* infile_name, char* outfile_name) {

	double bin_width = 0.05;		//GeV
	double mass_range_lb = 1.3;		//GeV
	double mass_range_ub = 4.;		//GeV

	using hlib::PION_MASS;

	if((int)(mass_range_ub*1000. - mass_range_lb*1000.) % (int)(bin_width*1000.)) {
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
	tfiles.resize(bounds.size()-1, NULL);
	for(unsigned int i = 0; i < tfiles.size(); ++i) {
		std::string dir_to_make = outfile_name;
		std::ostringstream strs;
		strs<<outfile_name<<"/"<<(int)(bounds.at(i)*1000.)<<"."<<(int)((bounds.at(i+1))*1000.);
		mkdir(strs.str().c_str(), S_IRWXU | S_IRWXG);
		strs<<"/"<<(int)(bounds.at(i)*1000.)<<"."<<(int)((bounds.at(i+1))*1000.)<<".root";
		tfiles.at(i) = TFile::Open(strs.str().c_str(), "NEW");
		if(tfiles.at(i) == NULL) {
			std::cout<<"Error opening file for writing."<<std::endl;
			return;
		}
	}

	// Open input file and do all the tree stuff
	TFile* infile = TFile::Open(infile_name, "READ");
	if(infile == NULL) {
		return;
	}
	TTree* tree = (TTree*)infile->Get("USR55");
	if(tree == NULL) {
		std::cout<<"Error opening in-TTree."<<std::endl;
		return;
	}
	double px1, py1, pz1;
	double px2, py2, pz2;
	double px3, py3, pz3;
	double px4, py4, pz4;
	double px5, py5, pz5;
	double gradx, grady;
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

	std::vector<TLorentzVector> particles;
	particles.resize(6);

	TClonesArray prodMomName("TObjString", 1);
	TClonesArray decayMomName("TObjString", 5);

	new (prodMomName [0]) TObjString("pi-");
	new (decayMomName [0]) TObjString("pi-");
	new (decayMomName [1]) TObjString("pi-");
	new (decayMomName [2]) TObjString("pi-");
	new (decayMomName [3]) TObjString("pi+");
	new (decayMomName [4]) TObjString("pi+");

	TClonesArray* prodMom = new TClonesArray("TVector3");
	TClonesArray* decayMom = new TClonesArray("TVector3");

	std::vector<TTree*> trees;
	trees.resize(tfiles.size(), NULL);

	const int splitLevel = 0;
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
		particles.at(0) = hlib::get_beam_energy(TVector3(gradx, grady, 1.), pSum);
		double mass = pSum.M();

		if(mass < bounds.at(0) || mass > bounds.at(bounds.size()-1)) {
			continue;
		}

		new ((*prodMom)[0]) TVector3(particles.at(0).Vect());
		for(unsigned int i = 1; i < 6; ++i) {
			new ((*decayMom)[i-1]) TVector3(particles.at(i).Vect());
		}

		for(unsigned int i = 0; i < tfiles.size(); ++i) {
			if((mass > bounds.at(i)) && (mass < bounds.at(i+1))) {
				trees.at(i)->Fill();
			}
		}

	} // End loop over evnts.

	for(unsigned int i = 0; i < tfiles.size(); ++i) {
		tfiles.at(i)->cd();
		trees.at(i)->Write();
		tfiles.at(i)->Close();
	}

}

int main(int argc, char* argv[]) {
	if(argc == 3) {
		convert_root_to_txt(argv[1], argv[2]);
	}
}
