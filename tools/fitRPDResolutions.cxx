
#include<iomanip>
#include<iostream>
#include<string>
#include<sstream>
#include<cmath>

#include<TApplication.h>
#include<TF1.h>
#include<TFile.h>
#include<TFitResult.h>
#include<TH1D.h>
#include<TLorentzVector.h>
#include<TTree.h>
#include<TFitter.h>

#include<basic_calcs.h>
#include<initializer.h>
#include<constants.h>

static bool useTransform = true;

class bin0r {

  public:

	bin0r() {
		_radiusLimits.push_back(0);
		for(unsigned int i = 1; i <= NUMBER_OF_RINGS; ++i) {
			_radiusLimits.push_back((MAXIMUM_RADIUS/NUMBER_OF_RINGS)*i);
		}
		for(unsigned int i = 0; i < _radiusLimits.size(); ++i) {
			_phiLimits.push_back(std::vector<double>());
			const double adjustifier = -TMath::Pi();
			_phiLimits[i].push_back(adjustifier);
			unsigned int limit = (i+1)*(i+1);
			for(unsigned int j = 1; j <= limit; ++j) {
				_phiLimits[i].push_back(((2*TMath::Pi()/limit)*j)+adjustifier);
			}
		}

		unsigned int binNo = 0;
		for(unsigned int i = 0; i < _radiusLimits.size()-1; ++i) {
			for(unsigned int j = 0; j < _phiLimits[i].size()-1; ++j) {
				std::stringstream sstr;
				sstr<<"Bin_"<<binNo<<"_r"<<_radiusLimits[i+1]<<"_phi"<<_phiLimits[i][j+1];
				_binNames.push_back(sstr.str());
				std::cout<<"Bin: "<<binNo++<<": r["<<_radiusLimits[i]
				         <<", "<<_radiusLimits[i+1]<<"], phi["<<_phiLimits[i][j]
				         <<", "<<_phiLimits[i][j+1]<<"]"<<std::endl;
			}
		}

		if(useTransform) {
			RPD_SLAB_PHI_ANGLES.push_back(-0.196);
			RPD_SLAB_PHI_ANGLES.push_back(0.002);
			RPD_SLAB_PHI_ANGLES.push_back(0.194);
			if(RPD_SLAB_PHI_ANGLES.size() != 3) {
				throw;
			}
		} else {
			RPD_SLAB_PHI_ANGLES.push_back(-3.140);
			RPD_SLAB_PHI_ANGLES.push_back(-2.945);
			RPD_SLAB_PHI_ANGLES.push_back(-2.817);
			RPD_SLAB_PHI_ANGLES.push_back(-2.619);
			RPD_SLAB_PHI_ANGLES.push_back(-2.422);
			RPD_SLAB_PHI_ANGLES.push_back(-2.290);
			RPD_SLAB_PHI_ANGLES.push_back(-2.098);
			RPD_SLAB_PHI_ANGLES.push_back(-1.900);
			RPD_SLAB_PHI_ANGLES.push_back(-1.768);
			RPD_SLAB_PHI_ANGLES.push_back(-1.570);
			RPD_SLAB_PHI_ANGLES.push_back(-1.378);
			RPD_SLAB_PHI_ANGLES.push_back(-1.246);
			RPD_SLAB_PHI_ANGLES.push_back(-1.048);
			RPD_SLAB_PHI_ANGLES.push_back(-0.850);
			RPD_SLAB_PHI_ANGLES.push_back(-0.718);
			RPD_SLAB_PHI_ANGLES.push_back(-0.526);
			RPD_SLAB_PHI_ANGLES.push_back(-0.329);
			RPD_SLAB_PHI_ANGLES.push_back(-0.196);
			RPD_SLAB_PHI_ANGLES.push_back(0.002);
			RPD_SLAB_PHI_ANGLES.push_back(0.194);
			RPD_SLAB_PHI_ANGLES.push_back(0.326);
			RPD_SLAB_PHI_ANGLES.push_back(0.524);
			RPD_SLAB_PHI_ANGLES.push_back(0.722);
			RPD_SLAB_PHI_ANGLES.push_back(0.848);
			RPD_SLAB_PHI_ANGLES.push_back(1.046);
			RPD_SLAB_PHI_ANGLES.push_back(1.244);
			RPD_SLAB_PHI_ANGLES.push_back(1.376);
			RPD_SLAB_PHI_ANGLES.push_back(1.568);
			RPD_SLAB_PHI_ANGLES.push_back(1.766);
			RPD_SLAB_PHI_ANGLES.push_back(1.898);
			RPD_SLAB_PHI_ANGLES.push_back(2.096);
			RPD_SLAB_PHI_ANGLES.push_back(2.288);
			RPD_SLAB_PHI_ANGLES.push_back(2.420);
			RPD_SLAB_PHI_ANGLES.push_back(2.618);
			RPD_SLAB_PHI_ANGLES.push_back(2.816);
			RPD_SLAB_PHI_ANGLES.push_back(2.943);
			if(RPD_SLAB_PHI_ANGLES.size() != 36) {
				throw;
			}
		}
	}

	void print() const {
		unsigned int binNo = 0;
		for(unsigned int i = 0; i < _radiusLimits.size()-1; ++i) {
			for(unsigned int j = 0; j < _phiLimits[i].size()-1; ++j) {
				std::cout<<"Bin: "<<binNo++<<": r["<<_radiusLimits[i]
				         <<", "<<_radiusLimits[i+1]<<"], phi["<<_phiLimits[i][j]
				         <<", "<<_phiLimits[i][j+1]<<"]"<<std::endl;
			}
		}
	}

	bool getBinNo(const double& r, const double& phi, unsigned int& binNo) const {
		if(r > _radiusLimits[_radiusLimits.size()-1]) {
			return false;
		}
		binNo = 0;
		unsigned int radiusI = 1;
		for(; r > _radiusLimits[radiusI]; ++radiusI) {
			binNo += (radiusI)*(radiusI);
		}
		for(unsigned int phiI = 1; phi > _phiLimits[radiusI-1][phiI]; ++phiI) {
			++binNo;
		}
		return true;
	}

	TH1D* getHist(const double& r, const double& phi, const double& rpdPhi) {
		unsigned int binNo;
		bool success = getBinNo(r, phi, binNo);
		if(not success) {
			return 0;
		}
		unsigned int histIndex = 0;
		for(; rpdPhi > RPD_SLAB_PHI_ANGLES[histIndex] + RPD_SLAB_STEP; ++histIndex);
		return _histMap[binNo][histIndex];
	}

	std::vector<TH1D*> getHistsInBin(const int& binNo) {
		return _histMap[binNo];
	}

	unsigned int getNHistsInBin() const {
		std::map<unsigned int, std::vector<TH1D*> >::const_iterator it = _histMap.find(0);
		return it->second.size();
	}

	unsigned int getNBins() const {
		const unsigned int& n = _radiusLimits.size() - 1;
		return ((2*n*n*n + 3*n*n + n)/6);
	}

	void prepareOutFile(TDirectory* dir) {
		assert(_binNames.size() == getNBins());
		for(unsigned int i = 0; i < getNBins(); ++i) {
			dir->cd();
			dir->mkdir(_binNames[i].c_str());
			dir->cd(_binNames[i].c_str());
			for(unsigned int j = 0; j < RPD_SLAB_PHI_ANGLES.size(); ++j) {
				std::stringstream sstr;
				sstr<<"hist_bin"<<i<<"_phi"<<RPD_SLAB_PHI_ANGLES[j];
				_histMap[i].push_back(new TH1D(sstr.str().c_str(), sstr.str().c_str(), 1000, -0.5, 0.5));
			}
		}
	}

  private:

	static const unsigned int NUMBER_OF_RINGS = 6;
	static const double MAXIMUM_RADIUS = 1.60;

	static const double RPD_SLAB_STEP = 0.025;

	std::vector<double> _radiusLimits;
	std::vector<std::vector<double> > _phiLimits;
	std::map<unsigned int, std::vector<TH1D*> > _histMap;
	std::vector<double> RPD_SLAB_PHI_ANGLES;

	std::vector<std::string> _binNames;

};

void fitErrorFunctionForRPD(std::string inFileName, std::string outFileName, std::string configFileName) {

	antok::Initializer* initializer = antok::Initializer::instance();
	if(not initializer->readConfigFile(configFileName)) {
		std::cerr<<"Could not open config file. Aborting..."<<std::endl;
		exit(1);
	}
	const double& PION_MASS = antok::Constants::chargedPionMass();

	new TApplication("app", 0, 0);

	TFile* inFile = TFile::Open(inFileName.c_str(), "READ");

	TTree* inTree = (TTree*)inFile->Get("RPDFitsTree/USR55");

	double gradx, grady;
	inTree->SetBranchAddress("gradx", &gradx);
	inTree->SetBranchAddress("grady", &grady);

	double rpdProtonX, rpdProtonY, rpdProtonZ, rpdProtonE;
	inTree->SetBranchAddress("RPD_Px", &rpdProtonX);
	inTree->SetBranchAddress("RPD_Py", &rpdProtonY);
	inTree->SetBranchAddress("RPD_Pz", &rpdProtonZ);
	inTree->SetBranchAddress("RPD_E", &rpdProtonE);

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

	double vertexX, vertexY, vertexZ;
	inTree->SetBranchAddress("X_primV", &vertexX);
	inTree->SetBranchAddress("Y_primV", &vertexY);
	inTree->SetBranchAddress("Z_primV", &vertexZ);

	TFile* outFile = TFile::Open(outFileName.c_str(), "NEW");
	if(not outFile) {
		std::cerr<<"Could not open outfile. Aborting..."<<std::endl;
		return;
	}

	bin0r bins;

	bins.prepareOutFile(outFile);

	long nBins = inTree->GetEntries();

	unsigned int roundingNumber = int(std::pow(10., (unsigned int)(log10((double)nBins / 100.) + 0.5)) + 0.5);

	const double rpdPeriod = TMath::Pi() / 6.;

	for(long i = 0; i < nBins; ++i) {

		if(not (i % roundingNumber)) {
			std::cout<<"Event "<<i<<" of "<<nBins<<" ("<<std::setprecision(2)
			         <<(i/(double)nBins*100)<<"%)"<<std::endl;
		}

		inTree->GetEntry(i);

		p1.SetXYZM(px1, py1, pz1, PION_MASS);
		p2.SetXYZM(px2, py2, pz2, PION_MASS);
		p3.SetXYZM(px3, py3, pz3, PION_MASS);
		p4.SetXYZM(px4, py4, pz4, PION_MASS);
		p5.SetXYZM(px5, py5, pz5, PION_MASS);

		TLorentzVector xVector = *(allMomenta[0]);
		for(unsigned int j = 1; j < allMomenta.size(); ++j) {
			xVector += *(allMomenta[j]);
		}

		TVector3 p3Beam;
		TLorentzVector pBeam;
		p3Beam.SetXYZ(gradx, grady, 1.);
		pBeam = antok::getBeamEnergy(p3Beam, xVector);
		p3Beam = pBeam.Vect();

/*		double t = std::fabs((pBeam - xVector).Mag2());
		double tMin = std::fabs((std::pow(xVector.M2() - pBeam.M2(), 2)) / (4. * p3Beam.Mag2()));
		double tPrime = t - tMin;
*/
		TLorentzVector rpdProton(rpdProtonX, rpdProtonY, rpdProtonZ, rpdProtonE);

		double deltaPhi, res;
		double protonPhi, xPhi;
		antok::getRPDDeltaPhiResPrediction(pBeam, rpdProton, xVector, deltaPhi, res, protonPhi, xPhi);

		TVector3 planarVertex(vertexX, vertexY, 0.);

		double correctedVertexPhi = planarVertex.Phi();
		double correctedProtonPhi = protonPhi;
		double correctedDeltaPhi = deltaPhi;
		if(useTransform) {
			while(correctedProtonPhi > 0.25) {
				correctedProtonPhi -= rpdPeriod;
				correctedVertexPhi -= rpdPeriod;

			}
			while(correctedProtonPhi < -0.25) {
				correctedProtonPhi += rpdPeriod;
				correctedVertexPhi += rpdPeriod;
			}
			if(correctedVertexPhi < -TMath::Pi()) {
				correctedVertexPhi += TMath::TwoPi();
			}
			if(correctedVertexPhi > TMath::Pi()) {
				correctedVertexPhi -= TMath::TwoPi();
			}
			if(correctedDeltaPhi < -TMath::Pi()) {
				correctedDeltaPhi += TMath::TwoPi();
			}
			if(correctedDeltaPhi > TMath::Pi()) {
				correctedDeltaPhi -= TMath::TwoPi();
			}

			/*
			while(correctedVertexPhi < -TMath::Pi()) {
				correctedVertexPhi = TMath::Pi() - correctedVertexPhi;
			}
			while(correctedVertexPhi > TMath::Pi()) {
				correctedVertexPhi = correctedVertexPhi - TMath::Pi();
			}
			while(correctedDeltaPhi < -TMath::Pi()) {
				correctedDeltaPhi = TMath::Pi() - correctedDeltaPhi;
			}
			while(correctedDeltaPhi > TMath::Pi()) {
				correctedDeltaPhi = correctedDeltaPhi - TMath::Pi();
			}
			*/
		}
/*
		std::cout<<"gradx="<<gradx<<" grady="<<grady<<std::endl;
		pBeam.Print();
		rpdProton.Print();
		xVector.Print();
		std::cout<<"--------"<<std::endl;
		std::cout<<"r(v)="<<planarVertex.Mag()<<std::endl;
		std::cout<<"phi(v)="<<planarVertex.Phi()<<std::endl;
		std::cout<<"phi(p+)="<<protonPhi<<std::endl;
		std::cout<<"deltaPhi="<<deltaPhi<<std::endl;
		std::cout<<"deltaPhi'="<<correctedDeltaPhi<<std::endl;
		std::cout<<"phi'(v)="<<correctedVertexPhi<<std::endl;
		std::cout<<"phi'(p+)="<<correctedProtonPhi<<std::endl;
		std::cout<<"########"<<std::endl;
*/

		TH1D* hist = bins.getHist(planarVertex.Mag(), correctedVertexPhi, correctedProtonPhi);
		if(not hist) {
			continue;
		}
		hist->Fill(correctedDeltaPhi);

	}

	std::cout<<std::endl;
	std::cout<<"#---------------------#"<<std::endl;
	std::cout<<"| Starting fitting... |"<<std::endl;
	std::cout<<"#---------------------#"<<std::endl;
	std::cout<<std::endl;

	const unsigned int nFits = bins.getNBins() * bins.getNHistsInBin();
	unsigned int currentFit = 0;

	roundingNumber = int(std::pow(10., (unsigned int)(log10((double)nFits / 100.) + 0.5)) + 0.5);

	for(unsigned int i = 0; i < bins.getNBins(); ++i) {
		std::vector<TH1D*> histsInBin = bins.getHistsInBin(i);
		for(unsigned int j = 0; j < histsInBin.size(); ++j) {

			if(not (currentFit % roundingNumber)) {
				std::cout<<"Fit "<<currentFit<<" of "<<nFits<<" ("<<std::setprecision(2)
				         <<(currentFit/(double)nFits*100)<<"%)"<<std::endl;
			}

			TH1D* fitHist = histsInBin[j];
			TF1* fitFunction = new TF1("fitFunction", "(-[0])*(TMath::Erf((x-[1])/[3])-TMath::Erf((x+[2])/[4])+[5]*(TMath::Erf((x-[2])/[6])-TMath::Erf((x+[4])/[7])))", -0.3, 0.3);
			fitFunction->SetParameter(0, fitHist->GetBinContent(fitHist->GetMaximumBin()) / 2.);        // overall strength
			TFitter::SetMaxIterations(20000);
			fitFunction->SetParLimits(0, fitHist->GetBinContent(fitHist->GetMaximumBin()) / 20., 10000);
			fitFunction->SetParameter(1, 0.130899694); // position both erf pairs right
			fitFunction->SetParLimits(1, 0, 0.5);
			fitFunction->SetParameter(2, 0.130899694); // position both erf pairs left
			fitFunction->SetParLimits(2, 0, 0.5);
			fitFunction->SetParameter(3, 0.025);       // sigma 1st erf pair right
			fitFunction->SetParLimits(3, 0, 5);
			fitFunction->SetParameter(4, 0.025);       // sigma 1st erf pair left
			fitFunction->SetParLimits(4, 0, 5);
			fitFunction->SetParameter(5, 0.02);        // relative strength 2nd erf pair
			fitFunction->SetParLimits(5, 0., 0.5);
			fitFunction->SetParameter(6, 0.1);         // sigma 2nd erf pair left
			fitFunction->SetParLimits(6, 0, 10);
			fitFunction->SetParameter(7, 0.1);        // sigma 2nd erf pair left
			fitFunction->SetParLimits(7, 0, 10);

			TFitResultPtr result = fitHist->Fit(fitFunction, "RSEMI");
			int status = result->Status();
			std::cout<<status<<std::endl;
			currentFit++;
		}
	}

	outFile->Write();
	outFile->Close();

}

int main(int argc, char** argv) {

	fitErrorFunctionForRPD(argv[1], argv[2], argv[3]);

}
