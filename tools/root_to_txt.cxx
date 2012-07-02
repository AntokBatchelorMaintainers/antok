
#include<fstream>
#include<iostream>
#include<iomanip>

#include<TFile.h>
#include<TLorentzVector.h>
#include<TTree.h>

#include<basic_calcs.h>
#include<constants.h>

void convert_root_to_txt(char* infile_name, char* outfile_name) {

	using hlib::PION_MASS;

	TFile* infile = TFile::Open(infile_name, "READ");
	if(infile == NULL) {
		return;
	}

	ofstream outfile;
	outfile.open(outfile_name);

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

		outfile<<"6\n";
		for(unsigned int i = 0; i < particles.size(); ++i) {
			std::string header;
			if(i < 4) {
				header = "9 -1 ";
			} else {
				header = "8 1 ";
			}
			outfile<<header;
			outfile<<std::setprecision(20)<<particles.at(i).X()<<" ";
			outfile<<std::setprecision(20)<<particles.at(i).Y()<<" ";
			outfile<<std::setprecision(20)<<particles.at(i).Z()<<" ";
			outfile<<std::setprecision(20)<<particles.at(i).T()<<" \n";
		}

	}

	outfile.close();

}

int main(int argc, char* argv[]) {
	if(argc == 3) {
		convert_root_to_txt(argv[1], argv[2]);
	}
}
