#ifndef ANTOK_DATA_H
#define ANTOK_DATA_H

#include <map>
#include <set>
#include <string>

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVectorD.h"

namespace antok {

	class Initializer;
	class Function;

	class Data {

		friend class Initializer;

	public:

		virtual ~Data()
		{
			// delete all allocated memory
			for (const auto& v2Pair        : _vector2s)             delete v2Pair.second;
			for (const auto& v3Pair        : _vector3s)             delete v3Pair.second;
			for (const auto& LVPair        : _lorentzVectors)       delete LVPair.second;
			for (const auto& vDPair        : _vectorDs)             delete vDPair.second;
			for (const auto& intVecPair    : _intVectors)           delete intVecPair.second;
			for (const auto& l64tVecPair   : _long64_tVectors)      delete l64tVecPair.second;
			for (const auto& doubleVecPair : _doubleVectors)        delete doubleVecPair.second;
			for (const auto& vecDVecPair   : _vectorDoubleVectors)  delete vecDVecPair.second;
			for (const auto& vecIVecPair   : _vectorIntVectors)     delete vecIVecPair.second;
			for (const auto& v2VecPair     : _vector2Vectors)       delete v2VecPair.second;
			for (const auto& v3VecPair     : _vector3Vectors)       delete v3VecPair.second;
			for (const auto& vTDVecPair    : _vectorTDVectors)      delete vTDVecPair.second;
			for (const auto& LVVecPair     : _lorentzVectorVectors) delete LVVecPair.second;
		};

		template <typename T> bool insert             (const std::string& name);
		template <typename T> bool insertInputVariable(const std::string& name);

		//TODO add const version of getAddr?
		template <typename T> T* const getAddr(const std::string& name);
		std::string                    getType(const std::string& name) const;

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

		// 2-, 3-, and 4-vectors
		std::map<std::string, TVector2*>       _vector2s;
		std::map<std::string, TVector3*>       _vector3s;
		std::map<std::string, TLorentzVector*> _lorentzVectors;

		// n-dim vectors
		std::map<std::string, TVectorD*>       _vectorDs;

		// std::vectors of numbers
		std::map<std::string, std::vector<int>*>       _intVectors;
		std::map<std::string, std::vector<Long64_t>*>  _long64_tVectors;
		std::map<std::string, std::vector<double>*>    _doubleVectors;

		// std::vector of std::vectors
		std::map<std::string, std::vector<std::vector<int>>*>     _vectorIntVectors;
		std::map<std::string, std::vector<std::vector<double>>*>  _vectorDoubleVectors;

		// std::vectors of 2-, 3-, and 4-vectors
		std::map<std::string, std::vector<TVector2>*>       _vector2Vectors;
		std::map<std::string, std::vector<TVector3>*>       _vector3Vectors;
		std::map<std::string, std::vector<TLorentzVector>*> _lorentzVectorVectors;

		// std::vectors of n-dim vectors
		std::map<std::string, std::vector<TVectorD>*> _vectorTDVectors;

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
