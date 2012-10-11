#ifndef ANTOK_DATA_H
#define ANTOK_DATA_H

#include<iostream>
#include<map>
#include<string>
#include<sstream>

#include<Rtypes.h>

#include<constants.h>


class TLorentzVector;

namespace antok {

	struct Data {

		std::map<std::string, std::string> global_map;

		std::map<std::string, double> doubles;
		std::map<std::string, int> ints;
		std::map<std::string, Long64_t> long64_ts;

		std::map<std::string, TLorentzVector> lorentzVectors;

		bool insertDouble(std::string name) {
			if(doubles.count(name) > 0) {
				return false;
			}
			global_map[name] = "double";
			doubles[name] = -8888.8;
			return true;
		}

		bool insertInt(std::string name) {
			if(ints.count(name) > 0) {
				return false;
			}
			global_map[name] = "int";
			ints[name] = -8888;
			return true;
		}

		bool insertLong64_t(std::string name) {
			if(long64_ts.count(name) > 0) {
				return false;
			}
			global_map[name] = "Long64_t";
			long64_ts[name] = -8888;
			return true;
		}

		bool insertLorentzVector(std::string name) {
			if(lorentzVectors.count(name) > 0) {
				return false;
			}
			global_map[name] = "TLorentzVector";
			lorentzVectors[name] = *(new TLorentzVector);
			return true;
		}

		double* getDoubleAddr(std::string name) {
			if(doubles.count(name) < 1) {
				return 0;
			}
			return &doubles[name];
		}

		int* getIntAddr(std::string name) {
			if(ints.count(name) < 1) {
				return 0;
			}
			return &ints[name];
		}

		Long64_t* getLong64_tAddr(std::string name) {
			if(long64_ts.count(name) < 1) {
				return 0;
			}
			return &long64_ts[name];
		}

		TLorentzVector* getLorentzVectorAddr(std::string name) {
			if(lorentzVectors.count(name) < 1) {
				return 0;
			}
			return &lorentzVectors[name];
		}

		static std::string getVariableInsertionErrorMsg(std::vector<std::string> quantityNames, std::string quantityName = "") {
			std::stringstream msgStream;
			if(quantityNames.size() > 1) {
				msgStream<<"Could not insert variable \""<<quantityName<<"\" when registering calculation for quantities \"[";
				for(unsigned int i = 0; i < quantityNames.size()-1; ++i) {
					msgStream<<quantityNames[i]<<", ";
				}
				msgStream<<quantityNames[quantityNames.size()-1]<<"]\" (double entry?)."<<std::endl;
			} else {
				getVariableInsertionErrorMsg(quantityNames[0]);
			}
			return msgStream.str();
		};

		static std::string getVariableInsertionErrorMsg(std::string variableName) {

			std::stringstream msgStream;
			msgStream<<"Could not insert variable \""<<variableName<<"\" (double entry?)."<<std::endl;
			return msgStream.str();

		};

// ---------------------------------------------- OLD STUFF

		int Run;
		int TrigMask;
		Long64_t EvNbr;
		int SpillNbr;

		double X_primV;
		double Y_primV;
		double Z_primV;

		double gradx;
		double grady;

		double beam_time;

		std::vector<double> Mom_x;
		std::vector<double> Mom_y;
		std::vector<double> Mom_z;

		std::vector<double> z_max;

		std::vector<double> theta_RICH;
		std::vector<int> PID_RICH;

		double chi2PV;

		double RPD_Px;
		double RPD_Py;
		double RPD_Pz;
		double RPD_E;
		double RPD_Tz;
		double RPD_z;
		double RPD_beta;
		double RPD_Phi;
		double RPD_dEA;
		double RPD_dEB;
		int nbrRPDTracks;

		int cedarID_bayes;
		double cedarTheta_X;
		double cedarTheta_Y;
		double cedarProbK1;
		double cedarProbK2;
		double cedarProbK3;

		Data() {
			const unsigned int& N_PARTICLES = antok::Constants::n_particles();
			Mom_x.resize(N_PARTICLES);
			Mom_y.resize(N_PARTICLES);
			Mom_z.resize(N_PARTICLES);
			z_max.resize(N_PARTICLES);
			theta_RICH.resize(N_PARTICLES);
			PID_RICH.resize(N_PARTICLES);
		}

	};

}

#endif
