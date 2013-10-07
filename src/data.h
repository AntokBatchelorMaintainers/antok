#ifndef ANTOK_DATA_H
#define ANTOK_DATA_H

#include<map>
#include<string>
#include<utility>

#include<Rtypes.h>

class TLorentzVector;
class TVector3;

namespace antok {

	class Initializer;
	class Function;

	class Data {

		friend class Initializer;

	  public:

		template<typename T> bool insert(std::string name);

		template<typename T> T* getAddr(std::string name);

		std::string getType(std::string name);

		static std::string getVariableInsertionErrorMsg(std::vector<std::string> quantityNames,
		                                                std::string quantityName = "");
		static std::string getVariableInsertionErrorMsg(std::string variableName);

	  private:

		std::map<std::string, std::string> global_map;

		std::map<std::string, double> doubles;
		std::map<std::string, int> ints;
		std::map<std::string, Long64_t> long64_ts;

		std::map<std::string, std::vector<double> > doubleVectors;

		std::map<std::string, TLorentzVector> lorentzVectors;
		std::map<std::string, TVector3> vectors;


	};

}

#endif

