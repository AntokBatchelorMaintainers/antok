#ifndef ANTOK_DATA_H
#define ANTOK_DATA_H

#include<map>
#include<string>
#include<utility>

#include<Rtypes.h>

#include<generators_functions.h>

class TLorentzVector;

namespace antok {

	class Initializer;
	class Function;

	class Data {

		friend class Initializer;

	  public:

		bool insertDouble(std::string name);
		bool insertInt(std::string name);
		bool insertLong64_t(std::string name);
		bool insertLorentzVector(std::string name);

		double* getDoubleAddr(std::string name);
		int* getIntAddr(std::string name);
		Long64_t* getLong64_tAddr(std::string name);
		TLorentzVector* getLorentzVectorAddr(std::string name);

		std::string getType(std::string name);

		static std::string getVariableInsertionErrorMsg(std::vector<std::string> quantityNames,
		                                                std::string quantityName = "");
		static std::string getVariableInsertionErrorMsg(std::string variableName);

	  private:

		std::map<std::string, std::string> global_map;

		std::map<std::string, double> doubles;
		std::map<std::string, int> ints;
		std::map<std::string, Long64_t> long64_ts;

		std::map<std::string, TLorentzVector> lorentzVectors;


	};

}

#endif
