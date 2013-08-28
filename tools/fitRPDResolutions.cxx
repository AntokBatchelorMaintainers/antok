
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

	bin0r(TDirectory* dir) : _dir(dir) {
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

	TVector2 getBinCenter(const TH1D* hist) {
		return getBinCenter(getBinNumber(hist));
	}

	unsigned int getBinNumber(const TH1D* hist) {
		return _reverseHistMap[hist];
	}

	std::pair<double, double> getLimitsForFit(const unsigned int& binNumber, const TH1D* hist) {
		antok::RpdHelperHelper* rpdHelperHelper = antok::RpdHelperHelper::getInstance();
		TVector2 vertexXY = getBinCenter(hist);
		double rpdPhi;
		std::vector<TH1D*> histsInBin = getHistsInBin(getBinNumber(hist));
		assert(histsInBin.size() == 3);
		if(hist == histsInBin[0]) {
			rpdPhi = -0.196;
		} else if(hist == histsInBin[1]) {
			rpdPhi = 0.002;
		} else if(hist == histsInBin[2]) {
			rpdPhi = 0.194;
		} else {
			std::cerr<<"Could not find histogram. Aborting..."<<std::endl;
			throw;
		}
		return rpdHelperHelper->getLimits(rpdPhi, vertexXY);
	}

	void addFitResults(const unsigned int& binNumber, std::vector<TFitResultPtr> fitResults) {
		if(binNumber != _fitResults.size()) {
			std::cerr<<"Fit results have to be added sequentially. Aborting..."<<std::endl;
			throw;
		}
		_fitResults.push_back(fitResults);
	}

	void writeFitResultsToRoot(double* median1 = 0,
	                           double* median2 = 0,
	                           double* median3 = 0,
	                           double* median4 = 0,
	                           double* median5 = 0) {
		bool medians = (median1 != 0) and (median2 != 0) and
		               (median3 != 0) and (median4 != 0) and
		               (median5 != 0);
		std::vector<TFitResultPtr> fitResults = serializeFitResult();
		_dir->cd();
		TH1D* sigma1Hist[3] = { new TH1D("sigma11Hist", "sigma11Hist", 1000, 0, 0.5),
				                new TH1D("sigma12Hist", "sigma12Hist", 1000, 0, 0.5),
				                new TH1D("sigma13Hist", "sigma13Hist", 1000, 0, 0.5) };
		TH1D* sigma2Hist[3] = { new TH1D("sigma21Hist", "sigma21Hist", 1000, 0, 0.5),
				                new TH1D("sigma22Hist", "sigma22Hist", 1000, 0, 0.5),
				                new TH1D("sigma23Hist", "sigma23Hist", 1000, 0, 0.5) };
		TH1D* sigma3Hist = new TH1D("sigma3Hist", "sigma3Hist", 1000, 0, 0.5);
		TH1D* sigma4Hist = new TH1D("sigma4Hist", "sigma4Hist", 1000, 0, 0.5);
		TH1D* relativeContributionHist = new TH1D("relativeContributionHist", "relativeContributionHist", 5000, 0, 1.);
		double results1[getNBins()*getNHistsInBin()]; // sharp sigmas left
		double results2[getNBins()*getNHistsInBin()]; // sharp sigmas right
		double results3[getNBins()*getNHistsInBin()]; // broad sigmas left
		double results4[getNBins()*getNHistsInBin()]; // broad sigmas right
		double results5[getNBins()*getNHistsInBin()]; // relative contribution broad sigmas
		unsigned int fitResultIndex = 0;
		for(unsigned int i = 0; i < getNBins(); ++i) {
			_dir->cd(_binNames[i].c_str());
			for(unsigned int j = 0; j < RPD_SLAB_PHI_ANGLES.size(); ++j) {
				std::stringstream sstr;
				sstr<<"fitResult_"<<i<<"_phi"<<RPD_SLAB_PHI_ANGLES[j];
				fitResults[fitResultIndex]->Write(sstr.str().c_str());
				sigma1Hist[j]->Fill(fitResults[fitResultIndex]->Parameter(4));
				sigma2Hist[j]->Fill(fitResults[fitResultIndex]->Parameter(3));
				sigma3Hist->Fill(fitResults[fitResultIndex]->Parameter(7));
				sigma4Hist->Fill(fitResults[fitResultIndex]->Parameter(6));
				relativeContributionHist->Fill(fitResults[fitResultIndex]->Parameter(5));
				if(medians) {
					results1[fitResultIndex] = fitResults[fitResultIndex]->Parameter(4);
					results2[fitResultIndex] = fitResults[fitResultIndex]->Parameter(3);
					results3[fitResultIndex] = fitResults[fitResultIndex]->Parameter(7);
					results4[fitResultIndex] = fitResults[fitResultIndex]->Parameter(6);
					results5[fitResultIndex] = fitResults[fitResultIndex]->Parameter(5);
				}
				++fitResultIndex;
			}
		}
		if(medians) {
			(*median1) = TMath::Median(getNBins()*getNHistsInBin(), results1);
			(*median2) = TMath::Median(getNBins()*getNHistsInBin(), results2);
			(*median3) = TMath::Median(getNBins()*getNHistsInBin(), results3);
			(*median4) = TMath::Median(getNBins()*getNHistsInBin(), results4);
			(*median5) = TMath::Median(getNBins()*getNHistsInBin(), results5);
		}
	}

	void writeFitResultsAsGraphs() {
		std::vector<TFitResultPtr> fitResults = serializeFitResult();
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

		std::vector<TGraph2D*> graphs;
		_dir->cd();
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

	void writeFitResultsToAscii(std::string fileNamePrefix) {
		std::vector<TFitResultPtr> fitResults = serializeFitResult();
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

	void prepareOutFile() {
		assert(_binNames.size() == getNBins());
		for(unsigned int i = 0; i < getNBins(); ++i) {
			_dir->cd();
			_dir->mkdir(_binNames[i].c_str());
			_dir->cd(_binNames[i].c_str());
			for(unsigned int j = 0; j < RPD_SLAB_PHI_ANGLES.size(); ++j) {
				std::stringstream sstr;
				sstr<<"hist_bin"<<i<<"_phi"<<RPD_SLAB_PHI_ANGLES[j];
				TH1D* hist = new TH1D(sstr.str().c_str(), sstr.str().c_str(), 1000, -0.5, 0.5);
				_histMap[i].push_back(hist);
				_reverseHistMap[hist] = i;
			}
		}
	}

  private:

	std::vector<TFitResultPtr> serializeFitResult() {
		std::vector<TFitResultPtr> retval;
		for(unsigned int i = 0; i < _fitResults.size(); ++i) {
			for(unsigned int j = 0; j < _fitResults[i].size(); ++j) {
				retval.push_back(_fitResults[i][j]);
			}
		}
		return retval;
	}

	static const unsigned int NUMBER_OF_RINGS = 6;
	static const double MAXIMUM_RADIUS = 1.60;

	static const double RPD_SLAB_STEP = 0.025;

	std::vector<double> _radiusLimits;
	std::vector<std::vector<double> > _phiLimits;
	std::map<unsigned int, std::vector<TH1D*> > _histMap;
	std::vector<double> RPD_SLAB_PHI_ANGLES;

	std::map<unsigned int, std::pair<unsigned int, unsigned int> > _binNoToIndexMap;

	std::map<const TH1D*, unsigned int> _reverseHistMap;

	std::vector<std::string> _binNames;

	TDirectory* _dir;

	std::vector<std::vector<TFitResultPtr> > _fitResults;

};

class uberBin0r {

  public:

	uberBin0r(TDirectory* dir) : _dir(dir) {
		double tBinWidth = (LIMITS_T_UPPER - LIMITS_T_LOWER) / NBINS_T;
		double mBinWidth = (LIMITS_M_UPPER - LIMITS_M_LOWER) / NBINS_M;
		for(unsigned int i = 0; i <= NBINS_T; ++i) {
			_tLimits.push_back(LIMITS_T_LOWER + i*tBinWidth);
		}
		for(unsigned int i = 0; i <= NBINS_M; ++i) {
			_mLimits.push_back(LIMITS_M_LOWER + i*mBinWidth);
		}
		unsigned int binNumber = 0;
		for(unsigned int i = 1; i <= NBINS_T; ++i) {
			for(unsigned int j = 1; j <= NBINS_M; ++j) {
				std::stringstream sstr;
				sstr<<"Bin_"<<binNumber++;
				sstr<<"_t-"<<_tLimits[i];
				sstr<<"_m-"<<_mLimits[j];
				_binNames.push_back(sstr.str());
			}
		}
	};

	unsigned int getNBins() const {
		return NBINS_T * NBINS_M;
	}

	unsigned int getNHistsInBin() const {
		return _bins[0].getNBins() * _bins[0].getNHistsInBin();
	}

	std::vector<TH1D*> getHistsInBin(const unsigned int& binNo) {
		bin0r& bin = _bins[binNo];
		std::vector<TH1D*> retval;
		for(unsigned int i = 0; i < bin.getNBins(); ++i) {
			std::vector<TH1D*> histsInBin = bin.getHistsInBin(i);
			retval.insert(retval.end(), histsInBin.begin(), histsInBin.end());
		}
		return retval;
	}

	void prepareOutFile() {
		for(unsigned int i = 0; i < getNBins(); ++i) {
			TDirectory* subDir = _dir->mkdir(_binNames[i].c_str());
			_bins.push_back(bin0r(subDir));
			_bins[i].prepareOutFile();
		}
	}

	unsigned int getBinNumber(double t,
	                          double m)
	{
//		std::cout<<"original = ("<<t<<", "<<m<<", "<<l<<")"<<std::endl;
		if(t < _tLimits[0]) {
			t = (_tLimits[0] + _tLimits[1]) / 2.;
		} else if (t > _tLimits[_tLimits.size()-1]) {
			t = (_tLimits[_tLimits.size()-2] + _tLimits[_tLimits.size()-1]) / 2.;
		}
		if(m < _mLimits[0]) {
			m = (_mLimits[0] + _mLimits[1]) / 2.;
		} else if (m > _mLimits[_mLimits.size()-1]) {
			m = (_mLimits[_mLimits.size()-2] + _mLimits[_mLimits.size()-1]) / 2.;
		}
//		std::cout<<"adjusted = ("<<t<<", "<<m<<", "<<l<<")"<<std::endl;
//		std::stringstream sstr;
		unsigned int binNumber = 0;
		bool found = false;
		for(unsigned int i = 0; i < _tLimits.size()-1; ++i) {
			for(unsigned int j = 0; j < _tLimits.size()-1; ++j) {
/*					sstr<<"####### Bin "<<binNumber<<"#########"<<std::endl;
					sstr<<_tLimits[i]<<" <? "<<t<<" <? "<<_tLimits[i+1]<<std::endl;
					sstr<<_mLimits[j]<<" <? "<<m<<" <? "<<_mLimits[j+1]<<std::endl;
					sstr<<_lLimits[k]<<" <? "<<l<<" <? "<<_lLimits[k+1]<<std::endl;
					sstr<<"########################"<<std::endl;
*/				if((_tLimits[i] <= t and t < _tLimits[i+1]) and
				   (_mLimits[j] <= m and m < _mLimits[j+1]))
				{
					found = true;
					break;
				}
				++binNumber;
			}
			if(found) {
				break;
			}
		}
//		std::cout<<"Bin number = "<<binNumber<<std::endl;
		if(binNumber >= _bins.size()) {
			std::cerr<<"Bin number found out of range."<<std::endl;
//			std::cout<<sstr.str();
			throw;
		}
		return binNumber;
	}

	TVector2 getBinCenter(const unsigned int& binNumber, const TH1D* hist) {
		return _bins[binNumber].getBinCenter(hist);
	}

	TH1D* getHist(const double& r,
	              const double& phi,
	              const double& rpdPhi,
	              const double& t,
	              const double& m) {
//		std::cout<<"_bins.size()="<<_bins.size()<<std::endl;
//		std::cout<<"binNumber = "<<getBinNumber(t, m, l)<<std::endl;
		return _bins[getBinNumber(t, m)].getHist(r, phi, rpdPhi);
	}

	void addFitResults(const unsigned int& binNumber, std::vector<TFitResultPtr> fitResults) {
		assert(fitResults.size() == _bins[binNumber].getNBins() * _bins[binNumber].getNHistsInBin());
		const unsigned int histsInBin = _bins[binNumber].getNHistsInBin();
		std::vector<TFitResultPtr> fitResultsInner;
		unsigned int innerBinNumber = 0;
		for(unsigned int i = 0; i < fitResults.size(); ++i) {
			fitResultsInner.push_back(fitResults[i]);
			if(fitResultsInner.size() == histsInBin) {
				_bins[binNumber].addFitResults(innerBinNumber, fitResultsInner);
				fitResultsInner.clear();
				++innerBinNumber;
			}
		}
	}

	TVector2 getBinCoordinates(const unsigned int& binNumber) {
		unsigned int binIndex = 0;
		for(unsigned int i = 0; i < _tLimits.size()-1; ++i) {
			for(unsigned int j = 0; j < _mLimits.size()-1; ++j) {
				if(binIndex == binNumber) {
					double t = (_tLimits[i] + _tLimits[i+1]) / 2.;
					double m = (_mLimits[j] + _mLimits[j+1]) / 2.;
					return TVector2(t, m);
				}
				binIndex++;
			}
		}
		throw;
	}

	void writeFitResultsToRoot() {
		unsigned int nEntries = getNBins();
		double m[nEntries];
		double t[nEntries];
		double m1[nEntries];
		double m2[nEntries];
		double m3[nEntries];
		double m4[nEntries];
		double m5[nEntries];
		for(unsigned int i = 0; i < _bins.size(); ++i) {
			TVector2 binCoordinates = getBinCoordinates(i);
			m[i] = binCoordinates.X();
			t[i] = binCoordinates.Y();
			_bins[i].writeFitResultsToRoot(&(m1[i]), &(m2[i]), &(m3[i]), &(m4[i]), &(m5[i]));
		}
		_dir->cd();
		new TGraph2D("median_sharp_sigmas_left", "median_sharp_sigmas_left", nEntries, m, t, m1);
		new TGraph2D("median_sharp_sigmas_right", "median_sharp_sigmas_right", nEntries, m, t, m2);
		new TGraph2D("median_broad_sigmas_left", "median_broad_sigmas_left", nEntries, m, t, m3);
		new TGraph2D("median_broad_sigmas_right", "median_broad_sigmas_right", nEntries, m, t, m4);
		new TGraph2D("relative_contribution_broad_sigmas", "relative_contribution_broad_sigmas", nEntries, m, t, m5);
	}

	void writeFitResultsAsGraphs() {
		for(unsigned int i = 0; i < _bins.size(); ++i) {
			_bins[i].writeFitResultsAsGraphs();
		}
	}

	void writeFitResultsToAscii(std::string fileNamePrefix) {
		for(unsigned int i = 0; i < _bins.size(); ++i) {
			std::stringstream sstr;
			sstr<<fileNamePrefix<<_binNames[i]<<"_";
			_bins[i].writeFitResultsToAscii(sstr.str());
		}
	}

	std::pair<double, double> getLimitsForFit(const unsigned int& binNumber, const TH1D* hist) {
		return _bins[binNumber].getLimitsForFit(binNumber, hist);
	}

	std::pair<double, double> getSharpSigmasForBin(const unsigned int& binNumber) {
		TVector2 tm = getBinCoordinates(binNumber);
		antok::RpdHelperHelper* rpdHelperHelper = antok::RpdHelperHelper::getInstance();
		return rpdHelperHelper->getSharpSigmas(tm.Y(), tm.X());
	}

	std::pair<double, double> getBroadSigmasForBin(const unsigned int& binNumber) {
		TVector2 tm = getBinCoordinates(binNumber);
		antok::RpdHelperHelper* rpdHelperHelper = antok::RpdHelperHelper::getInstance();
		return rpdHelperHelper->getBroadSigmas(tm.Y(), tm.X());
	}

	double getBroadSigmaContributionForBin(const unsigned int& binNumber) {
		TVector2 tm = getBinCoordinates(binNumber);
		antok::RpdHelperHelper* rpdHelperHelper = antok::RpdHelperHelper::getInstance();
		return rpdHelperHelper->getBroadSigmasScalingFactor(tm.Y(), tm.X());
	}

	static double getL(const TVector3& vertex, const TVector3& proton) {
		static const double TARGET_RADIUS = antok::Constants::TargetRadius();

		const double a = proton.X()*proton.X() + proton.Y()*proton.Y();
		const double b = 2 * (proton.X()*vertex.X() + proton.Y()*vertex.Y());
		double c = vertex.X()*vertex.X() + vertex.Y()*vertex.Y() - TARGET_RADIUS*TARGET_RADIUS;
		double n = antok::utils::getPositiveSolutionOfQuadraticEquation(a, b, c);

		TVector3 protonInTarget = n * proton;
		TVector3 exitPoint = vertex + protonInTarget;

		static const double TARGET_ENLARGEMENT = 3.;
		static const double TARGET_DOWNSTREAM_EDGE = antok::Constants::TargetDownstreamEdge() + TARGET_ENLARGEMENT;
		static const double TARGET_UPSTREAM_EDGE = antok::Constants::TargetUpstreamEdge() - TARGET_ENLARGEMENT;
		double multiplier = 1.;
		if(exitPoint.Z() > TARGET_DOWNSTREAM_EDGE) {
			multiplier = (TARGET_DOWNSTREAM_EDGE - vertex.Z()) / protonInTarget.Z();
		}
		if(exitPoint.Z() < TARGET_UPSTREAM_EDGE) {
			multiplier = (vertex.Z() - TARGET_UPSTREAM_EDGE) / protonInTarget.Z();
		}

		if(multiplier < 0.) {
			return 0.;
		}
		double retval = multiplier * protonInTarget.Mag();
		if(std::isnan(retval)) {
			throw(42);
		}
		return multiplier * protonInTarget.Mag();
	}

  private:

	static const unsigned int NBINS_T = 5;
	static const unsigned int NBINS_M = 5;

	static const double LIMITS_T_LOWER = 0.3;
	static const double LIMITS_T_UPPER = 1.;
	static const double LIMITS_M_LOWER = 1.;
	static const double LIMITS_M_UPPER = 4.2;

	std::vector<bin0r> _bins;
	std::vector<std::string> _binNames;
	std::vector<double> _tLimits;
	std::vector<double> _mLimits;

	TDirectory* _dir;

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

	uberBin0r bins(outFile);
//	bin0r bins(outFile);

	bins.prepareOutFile();

	outFile->cd();
	TH1D* predictedProtonMag = new TH1D("predictedProtonMag", "predictedProtonMag", 1000, 0, 50);
//	TH1D* lengthInTarget = new TH1D("lengthInTarget", "lengthInTarget", 1000, 0, 50);
	TH1D* xMass = new TH1D("xMass", "xMass", 1000, 0, 7);

	long nBins = inTree->GetEntries();

	unsigned int roundingNumber = int(std::pow(10., (unsigned int)(log10((double)nBins / 100.) + 0.5)) + 0.5);

	const double rpdPeriod = antok::Constants::RPDPeriod();

	for(long i = 0; i < nBins; ++i) {

		if(not (i % roundingNumber)) {
			std::cout<<"Event "<<i<<" of "<<nBins<<" ("<<std::setprecision(2)
			         <<(i/(double)nBins*100)<<"%)"<<std::endl;
		}

//		if(i > (double)nBins/40.) break;

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
		xMass->Fill(xVector.M());

		TVector3 p3Beam;
		TLorentzVector pBeam;
		p3Beam.SetXYZ(gradx, grady, 1.);
		pBeam = antok::getBeamEnergy(p3Beam, xVector);
		p3Beam = pBeam.Vect();

		TVector3 proton = p3Beam - xVector.Vect();
		double t = proton.Mag();
		predictedProtonMag->Fill(proton.Mag());
/*
		TVector3 vertex(vertexX, vertexY, vertexZ);
		double l = 0.;
		try {
			l = uberBin0r::getL(vertex, proton);
		} catch(const int& e) {
			std::cerr<<"Length in target could not be calculated."<<std::endl;
			if(e != 42) {
				std::cerr<<"In fitErrorFunctionForRPD: Something went very wrong when "
				         <<"calculating the proton exit point. Aborting..."<<std::endl;
				throw;
			}
			continue;
		}
		lengthInTarget->Fill(l);
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

		}

//		TH1D* hist = bins.getHist(planarVertex.Mag(), correctedVertexPhi, correctedProtonPhi);
		TH1D* hist = bins.getHist(planarVertex.Mag(), correctedVertexPhi, correctedProtonPhi,
		                          t, xVector.M());
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


	std::vector<std::pair<unsigned int, unsigned int> > failedFitBins;

	for(unsigned int i = 0; i < bins.getNBins(); ++i) {
		std::vector<TH1D*> histsInBin = bins.getHistsInBin(i);
		std::vector<TFitResultPtr> fitResults;
		for(unsigned int j = 0; j < histsInBin.size(); ++j) {

			if(not (currentFit % roundingNumber)) {
				std::cout<<"Fit "<<currentFit<<" of "<<nFits<<" ("<<std::setprecision(2)
				         <<(currentFit/(double)nFits*100)<<"%)"<<std::endl;
			}

			TH1D* fitHist = histsInBin[j];
			TF1* fitFunction = new TF1("fitFunction", "(-[0])*(TMath::Erf((x-[1])/[3])-TMath::Erf((x+[2])/[4])+[5]*(TMath::Erf((x-[1])/[6])-TMath::Erf((x+[2])/[7])))", -0.4, 0.4);
			fitFunction->SetParameter(0, fitHist->GetBinContent(fitHist->GetMaximumBin()) / 2.);        // overall strength
			TFitter::SetMaxIterations(20000);
			fitFunction->SetParLimits(0, fitHist->GetBinContent(fitHist->GetMaximumBin()) / 20., 10000);

			std::pair<double, double> limits = bins.getLimitsForFit(i, fitHist);
			std::pair<double, double> sharpSigmas = bins.getSharpSigmasForBin(i);
			std::pair<double, double> broadSigmas = bins.getBroadSigmasForBin(i);
			double broadSigmaScalingFactor = bins.getBroadSigmaContributionForBin(i);
/*
			std::cout<<"sigmas.first="<<sigmas.first<<std::endl;
			std::cout<<"sigmas.second="<<sigmas.second<<std::endl;
*/
			fitFunction->FixParameter(1, limits.second); // position both erf pairs right
			fitFunction->FixParameter(2, -limits.first); // position both erf pairs left

			fitFunction->FixParameter(3, sharpSigmas.second);       // sigma 1st erf pair right
			fitFunction->FixParameter(4, sharpSigmas.first);       // sigma 1st erf pair left

			fitFunction->FixParameter(5, broadSigmaScalingFactor);        // relative strength 2nd erf pair

			fitFunction->FixParameter(6, broadSigmas.second);         // sigma 2nd erf pair right
			fitFunction->FixParameter(7, broadSigmas.first);         // sigma 2nd erf pair left

/*			fitFunction->SetParameter(1, limits.second);
			fitFunction->SetParameter(2, -limits.first);

			fitFunction->SetParameter(3, 0.025);       // sigma 1st erf pair right
			fitFunction->SetParLimits(3, 0, 5);
			fitFunction->SetParameter(4, 0.025);       // sigma 1st erf pair left
			fitFunction->SetParLimits(4, 0, 5);

			fitFunction->SetParameter(5, 0.02);        // relative strength 2nd erf pair
			fitFunction->SetParLimits(5, 0., 0.5);
			fitFunction->SetParameter(6, 0.1);         // sigma 2nd erf pair right
			fitFunction->SetParLimits(6, 0, 10);
			fitFunction->SetParameter(7, 0.1);         // sigma 2nd erf pair left
			fitFunction->SetParLimits(7, 0, 10);
*/
//			TFitResultPtr result = fitHist->Fit(fitFunction, "RSEMIL");
			TFitResultPtr result = fitHist->Fit(fitFunction, "RSLQ");
			fitResults.push_back(result);
			int status = result->Status();
			if(status != 0) {
				std::cout<<"Fit failed..."<<std::endl;
				failedFitBins.push_back(std::pair<unsigned int, unsigned int>(i, j));
			}
			currentFit++;
		}
		bins.addFitResults(i, fitResults);
	}

	if(failedFitBins.size() > 0) {
		std::cout<<"Fit failed in bins: ["<<failedFitBins[0].first<<"."<<failedFitBins[0].second;
		for(unsigned int i = 1; i < failedFitBins.size(); ++i) {
			std::cout<<", "<<failedFitBins[i].first<<"."<<failedFitBins[i].second;
		}
		std::cout<<"]"<<std::endl;
	}

	bins.writeFitResultsAsGraphs();
//	bins.writeFitResultsToAscii(graphDataFileNamePrefix);
	bins.writeFitResultsToRoot();

	outFile->Write();
	outFile->Close();

}

int main(int argc, char** argv) {

	fitErrorFunctionForRPD(argv[1], argv[2], argv[3], argv[4]);

}
