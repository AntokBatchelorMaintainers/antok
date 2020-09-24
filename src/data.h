#ifndef ANTOK_DATA_H
#define ANTOK_DATA_H

#include <map>
#include <set>
#include <string>

#include "Rtypes.h"

class TLorentzVector;
class TVector3;

namespace antok {

	class Initializer;
	class Function;

	class Data {

		friend class Initializer;

	public:

		//TODO add proper destructor that deletes all allocated memory

		template <typename T> bool insert             (const std::string& name);
		template <typename T> bool insertInputVariable(const std::string& name);

		//TODO add const version of getAddr?
		template <typename T> T* getAddr        (const std::string& name);
		std::string              getType        (const std::string& name) const;

		bool isInputVariable(const std::string& name) const;
		bool isVector       (const std::string& name) const;

		static std::string getVariableInsertionErrorMsg(const std::vector<std::string>& quantityNames,
		                                                const std::string&              quantityName = "");
		static std::string getVariableInsertionErrorMsg(const std::string& variableName);

	private:

		std::map<std::string, std::string> _global_map;
		std::set<std::string>              _inputVariables;  // bookkeeping which variable comes from the input tree

		// numbers
		std::map<std::string, int>      _ints;
		std::map<std::string, Long64_t> _long64_ts;
		std::map<std::string, double>   _doubles;

		// 3- and 4-vectors
		std::map<std::string, TVector3*>       _vector3s;
		std::map<std::string, TLorentzVector*> _lorentzVectors;

		// std::vectors
		std::map<std::string, std::vector<int>* >            _intVectors;
		std::map<std::string, std::vector<Long64_t>* >       _long64_tVectors;
		std::map<std::string, std::vector<double>* >         _doubleVectors;
		std::map<std::string, std::vector<TVector3>* >       _vector3Vectors;
		std::map<std::string, std::vector<TLorentzVector>* > _lorentzVectorVectors;


	};

}


template <typename T>
bool
antok::Data::insertInputVariable(const std::string& name )
{
	const bool ok = insert<T>(name);
	if (ok) {
		_inputVariables.insert(name);
	}
	return ok;
}


inline bool
antok::Data::isInputVariable(const std::string& name) const
{
	return _inputVariables.find(name) != _inputVariables.end();
}

#endif  // ANTOK_DATA_H
