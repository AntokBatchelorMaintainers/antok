
#include<fstream>
#include<iomanip>
#include<iostream>
#include<string>
#include<sstream>
#include<cmath>

#include<TApplication.h>
#include<TF1.h>
#include<TFile.h>
#include<TFitResult.h>
#include<TFitter.h>
#include<TGraph2D.h>
#include<TH1D.h>
#include<TLorentzVector.h>
#include<TTree.h>
#include<TVector2.h>
#include<TVector3.h>

#include<basic_calcs.h>
#include<constants.h>
#include<initializer.h>
#include<rpd_helper_helper.h>

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
				_binNoToIndexMap[binNo] = std::pair<unsigned int, unsigned int>(i, j);
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

	std::vector<TH1D*> getHistsInBin(const unsigned int& binNo) {
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

	TVector2 getBinCenter(const unsigned int& binNo) {
		TVector2 retval;
		std::pair<unsigned int, unsigned int> indices = _binNoToIndexMap[binNo];
		double r = 0;
		double phi = 0;
		if(binNo != 0) {
			r = (_radiusLimits[indices.first] + _radiusLimits[indices.first+1]) / 2.;
			phi = (_phiLimits[indices.first][indices.second] + _phiLimits[indices.first][indices.second+1]) / 2.;
		}
		retval.SetMagPhi(r, phi);
		return retval;
	}

	void writeFitResultsToRoot(TDirectory* dir, std::vector<TFitResultPtr> fitResults) {
		unsigned int fitResultIndex = 0;
		for(unsigned int i = 0; i < getNBins(); ++i) {
			dir->cd(_binNames[i].c_str());
			for(unsigned int j = 0; j < RPD_SLAB_PHI_ANGLES.size(); ++j) {
				std::stringstream sstr;
				sstr<<"fitResult_"<<i<<"_phi"<<RPD_SLAB_PHI_ANGLES[j];
				fitResults[fitResultIndex++]->Write(sstr.str().c_str());
			}
		}
	}

	void writeFitResultsAsGraphs(TDirectory* dir, const std::vector<TFitResultPtr>& fitResults) {
		unsigned int fitResultIndex = 0;
		unsigned int nEntries = getNBins();
		double x[nEntries];
		double y[nEntries];
		double z_ll[nEntries];
		double z_lu[nEntries];
		double z_ml[nEntries];
		double z_mu[nEntries];
		double z_rl[nEntries];
		double z_ru[nEntries];

		double z_lp[nEntries];
		double z_lw[nEntries];
		double z_mp[nEntries];
		double z_mw[nEntries];
		double z_rp[nEntries];
		double z_rw[nEntries];
		for(unsigned int i = 0; i < getNBins(); ++i) {
			TVector2 binCenter = getBinCenter(i);
			x[i] = binCenter.X();
			y[i] = binCenter.Y();
			for(unsigned int j = 0; j < RPD_SLAB_PHI_ANGLES.size(); ++j) {
				const TFitResultPtr fitResult = fitResults[fitResultIndex];
				double lowerBound = -fitResult->Parameter(2);
				double upperBound = fitResult->Parameter(1);
				double width = upperBound - lowerBound;
				double position = (upperBound + lowerBound) / 2.;
				if(j == 0) {
					z_ll[i] = lowerBound;
					z_lu[i] = upperBound;
					z_lp[i] = position;
					z_lw[i] = width;
				} else if(j == 1) {
					z_ml[i] = lowerBound;
					z_mu[i] = upperBound;
					z_mp[i] = position;
					z_mw[i] = width;
				} else if(j == 2) {
					z_rl[i] = lowerBound;
					z_ru[i] = upperBound;
					z_rp[i] = position;
					z_rw[i] = width;
				}
				fitResultIndex++;
			}
		}
/*		for(unsigned int i = 0; i < getNBins(); ++i) {
			std::cout<<"z_ll["<<i<<
		}
*/		std::vector<TGraph2D*> graphs;
		dir->cd();
		graphs.push_back(new TGraph2D("left_slab_lower_bound", "left_slab_lower_bound", nEntries, x, y, z_ll));
		graphs.push_back(new TGraph2D("left_slab_upper_bound", "left_slab_upper_bound", nEntries, x, y, z_lu));
		graphs.push_back(new TGraph2D("middle_slab_lower_bound", "middle_slab_lower_bound", nEntries, x, y, z_ml));
		graphs.push_back(new TGraph2D("middle_slab_upper_bound", "middle_slab_upper_bound", nEntries, x, y, z_mu));
		graphs.push_back(new TGraph2D("right_slab_lower_bound", "right_slab_lower_bound", nEntries, x, y, z_rl));
		graphs.push_back(new TGraph2D("right_slab_upper_bound", "right_slab_upper_bound", nEntries, x, y, z_ru));
		graphs.push_back(new TGraph2D("left_slab_position", "left_slab_position", nEntries, x, y, z_lp));
		graphs.push_back(new TGraph2D("left_slab_width", "left_slab_width", nEntries, x, y, z_lw));
		graphs.push_back(new TGraph2D("middle_slab_position", "middle_slab_position", nEntries, x, y, z_mp));
		graphs.push_back(new TGraph2D("middle_slab_width", "middle_slab_width", nEntries, x, y, z_mw));
		graphs.push_back(new TGraph2D("right_slab_position", "right_slab_position", nEntries, x, y, z_rp));
		graphs.push_back(new TGraph2D("right_slab_width", "right_slab_width", nEntries, x, y, z_rw));
	}

	void writeFitResultsToAscii(std::string fileNamePrefix, const std::vector<TFitResultPtr>& fitResults) {
		unsigned int fitResultIndex = 0;
		std::string fileName_ll = fileNamePrefix + "left_lower";
		std::string fileName_lu = fileNamePrefix + "left_upper";
		std::string fileName_ml = fileNamePrefix + "middle_lower";
		std::string fileName_mu = fileNamePrefix + "middle_upper";
		std::string fileName_rl = fileNamePrefix + "right_lower";
		std::string fileName_ru = fileNamePrefix + "right_upper";
		std::string fileName_lp = fileNamePrefix + "left_position";
		std::string fileName_lw = fileNamePrefix + "left_width";
		std::string fileName_mp = fileNamePrefix + "middle_position";
		std::string fileName_mw = fileNamePrefix + "middle_width";
		std::string fileName_rp = fileNamePrefix + "right_position";
		std::string fileName_rw = fileNamePrefix + "right_width";
		ofstream file_ll;
		ofstream file_lu;
		ofstream file_ml;
		ofstream file_mu;
		ofstream file_rl;
		ofstream file_ru;
		ofstream file_lp;
		ofstream file_lw;
		ofstream file_mp;
		ofstream file_mw;
		ofstream file_rp;
		ofstream file_rw;
		file_ll.open(fileName_ll.c_str());
		file_lu.open(fileName_lu.c_str());
		file_ml.open(fileName_ml.c_str());
		file_mu.open(fileName_mu.c_str());
		file_rl.open(fileName_rl.c_str());
		file_ru.open(fileName_ru.c_str());
		file_lp.open(fileName_lp.c_str());
		file_lw.open(fileName_lw.c_str());
		file_mp.open(fileName_mp.c_str());
		file_mw.open(fileName_mw.c_str());
		file_rp.open(fileName_rp.c_str());
		file_rw.open(fileName_rw.c_str());
		for(unsigned int i = 0; i < getNBins(); ++i) {
			TVector2 binCenter = getBinCenter(i);
			double x = binCenter.X();
			double y = binCenter.Y();
			for(unsigned int j = 0; j < RPD_SLAB_PHI_ANGLES.size(); ++j) {
				const TFitResultPtr fitResult = fitResults[fitResultIndex++];
				double lowerBound = -fitResult->Parameter(2);
				double upperBound = fitResult->Parameter(1);
				double width = upperBound - lowerBound;
				double position = (upperBound + lowerBound) / 2.;
				if(j == 0) {
					file_ll<<x<<" "<<y<<" "<<lowerBound<<"\n";
					file_lu<<x<<" "<<y<<" "<<upperBound<<"\n";
					file_lp<<x<<" "<<y<<" "<<position<<"\n";
					file_lw<<x<<" "<<y<<" "<<width<<"\n";
				} else if(j == 1) {
					file_ml<<x<<" "<<y<<" "<<lowerBound<<"\n";
					file_mu<<x<<" "<<y<<" "<<upperBound<<"\n";
					file_mp<<x<<" "<<y<<" "<<position<<"\n";
					file_mw<<x<<" "<<y<<" "<<width<<"\n";
				} else if(j == 2) {
					file_rl<<x<<" "<<y<<" "<<lowerBound<<"\n";
					file_ru<<x<<" "<<y<<" "<<upperBound<<"\n";
					file_rp<<x<<" "<<y<<" "<<position<<"\n";
					file_rw<<x<<" "<<y<<" "<<width<<"\n";
				}
			}
		}
		file_ll.close();
		file_lu.close();
		file_ml.close();
		file_mu.close();
		file_rl.close();
		file_ru.close();
		file_lp.close();
		file_lw.close();
		file_mp.close();
		file_mw.close();
		file_rp.close();
		file_rw.close();
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

	std::map<unsigned int, std::pair<unsigned int, unsigned int> > _binNoToIndexMap;

	std::vector<std::string> _binNames;

};

void fitErrorFunctionForRPD(std::string inFileName,
                            std::string outFileName,
                            std::string configFileName,
                            std::string graphDataFileNamePrefix)
{

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

	outFile->cd();
	TH1D* predictedProtonMag = new TH1D("predictedProtonMag", "predictedProtonMag", 1000, 0, 50);

	long nBins = inTree->GetEntries();

	unsigned int roundingNumber = int(std::pow(10., (unsigned int)(log10((double)nBins / 100.) + 0.5)) + 0.5);

	const double rpdPeriod = antok::Constants::RPDPeriod();

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

		TVector3 proton = p3Beam - xVector.Vect();
		predictedProtonMag->Fill(proton.Mag());
		if(proton.Mag() < 0.425) {
			continue;
		}

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

	std::vector<TFitResultPtr> fitResults;

	std::vector<std::pair<unsigned int, unsigned int> > failedFitBins;

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
/*			fitFunction->SetParameter(1, 0.130899694); // position both erf pairs right
			fitFunction->SetParLimits(1, -0.1, 0.5);
			fitFunction->SetParameter(2, 0.130899694); // position both erf pairs left
			fitFunction->SetParLimits(2, -0.1, 0.5);
*/
			antok::RpdHelperHelper* rpdHelperHelper = antok::RpdHelperHelper::getInstance();
			TVector2 vertexXY = bins.getBinCenter(i);
			double rpdPhi = 0.;
			switch(j) {
				case(0):
						std::cout<<"left"<<std::endl;
						rpdPhi = -0.196;
						break;
				case(1):
						std::cout<<"middle"<<std::endl;
						rpdPhi = 0.002;
						break;
				case(2):
						std::cout<<"right"<<std::endl;
						rpdPhi = 0.194;
						break;

			}
			std::pair<double, double> limits = rpdHelperHelper->getLimits(rpdPhi, vertexXY);
			std::cout<<"limits = "<<limits.first<<", "<<limits.second<<std::endl;
			const double& LIMIT_WINDOW = 0.01;
			std::cout<<"lower limit = "<<-limits.first<<" ["<<(-limits.first-LIMIT_WINDOW)<<", "<<(-limits.first+LIMIT_WINDOW)<<"]"<<std::endl;
			std::cout<<"upper limit = "<<limits.second<<" ["<<(limits.second-LIMIT_WINDOW)<<", "<<(limits.second+LIMIT_WINDOW)<<"]"<<std::endl;

			fitFunction->FixParameter(1, limits.second);
			fitFunction->FixParameter(2, -limits.first);

/*
			fitFunction->SetParameter(1, limits.second); // position both erf pairs right
			fitFunction->SetParLimits(1, limits.second-LIMIT_WINDOW, limits.second+LIMIT_WINDOW);
			fitFunction->SetParameter(2, -limits.first); // position both erf pairs left
			fitFunction->SetParLimits(2, -limits.first-LIMIT_WINDOW, -limits.first+LIMIT_WINDOW);
*/
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

//			TFitResultPtr result = fitHist->Fit(fitFunction, "RSEMIL");
			TFitResultPtr result = fitHist->Fit(fitFunction, "RSL");
			fitResults.push_back(result);
			int status = result->Status();
			std::cout<<status<<std::endl;
			if(status != 0) {
				std::cout<<"FIT FAILED FIT FAILED FIT FAILED FIT FAILED FIT FAILED FIT FAILED FIT FAILED FIT FAILED"<<std::endl;
				failedFitBins.push_back(std::pair<unsigned int, unsigned int>(i, j));
			}
			currentFit++;
		}
	}

	if(failedFitBins.size() > 0) {
		std::cout<<"Fit failed in bins: ["<<failedFitBins[0].first<<"."<<failedFitBins[0].second;
		for(unsigned int i = 1; i < failedFitBins.size(); ++i) {
			std::cout<<", "<<failedFitBins[i].first<<"."<<failedFitBins[i].second;
		}
		std::cout<<"]"<<std::endl;
	}

	bins.writeFitResultsAsGraphs(outFile, fitResults);
	bins.writeFitResultsToAscii(graphDataFileNamePrefix, fitResults);
	bins.writeFitResultsToRoot(outFile, fitResults);

	outFile->Write();
	outFile->Close();

}

int main(int argc, char** argv) {

	fitErrorFunctionForRPD(argv[1], argv[2], argv[3], argv[4]);

}
