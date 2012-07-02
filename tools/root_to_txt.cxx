
#include<fstream>

#include<TFile.h>

void convert_root_to_txt(char* infile_name, char* outfile_name) {

	TFile* infile = TFile::Open(infile_name, "READ");
	if(infile == NULL) {
		return;
	}

	ofstream outfile;
	outfile.open(outfile_name);

	TTree* tree = (TTree*)infile->Get("kbicker/USR55");
	if(tree == NULL) {
		std::cout<<"Error opening in-TTree."<<std::endl;
		return;
	}

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
	tree_chain->SetBranchAddress("gradx", &gradx);
	tree_chain->SetBranchAddress("grady", &grady);

	for(unsigned int i = 0; i < tree_chain->GetEntries(); ++i) {


	}


}

int main(int argc, char* argv[]) {
	if(argc == 3) {
		convert_root_to_txt(argv[1], argv[2]);
	}
}
