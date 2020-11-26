#include <iostream>
#include <sstream>

#include "TLorentzVector.h"

#include "data.h"

namespace antok {

	template <typename T>
	bool
	Data::insert(const std::string& name)
	{
		std::cerr << "Could not insert variable of this type." << std::endl;
		return false;
	}


	template <typename T>
	T*
	Data::getAddr(const std::string& name)
	{
		std::cerr << "Could not get address for variable of this type." << std::endl;
		return nullptr;
	}


	template <>
	bool
	Data::insert<int>(const std::string& name)
	{
		if (_ints.count(name) > 0) {
			return false;
		}
		_global_map[name] = "int";
		_ints      [name] = -8888;  //TODO make magic int configurable
		return true;
	}


	template <>
	bool
	Data::insert<Long64_t>(const std::string& name)
	{
		if (_long64_ts.count(name) > 0) {
			return false;
		}
		_global_map[name] = "Long64_t";
		_long64_ts [name] = -8888;  //TODO make magic int configurable
		return true;
	}


	template <>
	bool
	Data::insert<double>(const std::string& name)
	{
		if (_doubles.count(name) > 0) {
			return false;
		}
		_global_map[name] = "double";
		_doubles   [name] = -8888.8;  //TODO use nan
		return true;
	}


	template <>
	bool
	Data::insert<TVector3>(const std::string& name)
	{
		if (_vector3s.count(name) > 0) {
			return false;
		}
		_global_map[name] = "TVector3";
		_vector3s  [name] = new TVector3;
		return true;
	}


	template <>
	bool
	Data::insert<TLorentzVector>(const std::string& name)
	{
		if (_lorentzVectors.count(name) > 0) {
			return false;
		}
		_global_map    [name] = "TLorentzVector";
		_lorentzVectors[name] = new TLorentzVector;
		return true;
	}


	template <>
	bool
	Data::insert<std::vector<int> >(const std::string& name)
	{
		if (_intVectors.count(name) > 0) {
			return false;
		}
		_global_map[name] = "std::vector<int>";
		_intVectors[name] = new std::vector<int>;
		return true;
	}


	template <>
	bool
	Data::insert<std::vector<Long64_t> >(const std::string& name)
	{
		if (_long64_tVectors.count(name) > 0) {
			return false;
		}
		_global_map     [name] = "std::vector<Long64_t>";
		_long64_tVectors[name] = new std::vector<Long64_t>;
		return true;
	}


	template <>
	bool
	Data::insert<std::vector<double> >(const std::string& name)
	{
		if (_doubleVectors.count(name) > 0) {
			return false;
		}
		_global_map   [name] = "std::vector<double>";
		_doubleVectors[name] = new std::vector<double>;
		return true;
	}


	template <>
	bool
	Data::insert<std::vector<TVector3> >(const std::string& name)
	{
		if (_vector3Vectors.count(name) > 0) {
			return false;
		}
		_global_map    [name] = "std::vector<TVector3>";
		_vector3Vectors[name] = new std::vector<TVector3>;
		return true;
	}


	template <>
	bool
	Data::insert<std::vector<TLorentzVector> >(const std::string& name)
	{
		if (_lorentzVectorVectors.count(name) > 0) {
			return false;
		}
		_global_map          [name] = "std::vector<TLorentzVector>";
		_lorentzVectorVectors[name] = new std::vector<TLorentzVector>;
		return true;
	}


	template <>
	int*
	Data::getAddr<int>(const std::string& name)
	{
		auto it = _ints.find(name);
		if (it == _ints.end()) {
			return nullptr;
		}
		return &it->second;
	}


	template <>
	Long64_t*
	Data::getAddr<Long64_t>(const std::string& name)
	{
		auto it = _long64_ts.find(name);
		if (it == _long64_ts.end()) {
			return nullptr;
		}
		return &it->second;
	}


	template <>
	double*
	Data::getAddr<double>(const std::string& name)
	{
		auto it = _doubles.find(name);
		if (it == _doubles.end()) {
			return nullptr;
		}
		return &it->second;
	}


	template <>
	TVector3*
	Data::getAddr<TVector3>(const std::string& name)
	{
		auto it = _vector3s.find(name);
		if (it == _vector3s.end()) {
			return nullptr;
		}
		return it->second;
	}


	template <>
	TLorentzVector*
	Data::getAddr<TLorentzVector>(const std::string& name)
	{
		auto it = _lorentzVectors.find(name);
		if (it == _lorentzVectors.end()) {
			return nullptr;
		}
		return it->second;
	}


	template <>
	std::vector<int>*
	Data::getAddr<std::vector<int> >(const std::string& name)
	{
		auto it = _intVectors.find(name);
		if (it == _intVectors.end()) {
			return nullptr;
		}
		return it->second;
	}


	template <>
	std::vector<Long64_t>*
	Data::getAddr<std::vector<Long64_t>>(const std::string& name)
	{
		auto it = _long64_tVectors.find(name);
		if (it == _long64_tVectors.end()) {
			return nullptr;
		}
		return it->second;
	}


	template <>
	std::vector<double>*
	Data::getAddr<std::vector<double> >(const std::string& name)
	{
		auto it = _doubleVectors.find(name);
		if (it == _doubleVectors.end()) {
			return nullptr;
		}
		return it->second;
	}


	template <>
	std::vector<TVector3>*
	Data::getAddr<std::vector<TVector3> >(const std::string& name)
	{
		auto it = _vector3Vectors.find(name);
		if (it == _vector3Vectors.end()) {
			return nullptr;
		}
		return it->second;
	}


	template <>
	std::vector<TLorentzVector>*
	Data::getAddr<std::vector<TLorentzVector> >(const std::string& name)
	{
		auto it = _lorentzVectorVectors.find(name);
		if (it == _lorentzVectorVectors.end()) {
			return nullptr;
		}
		return it->second;
	}


	std::string
	Data::getType(const std::string& name) const
	{
		auto it = _global_map.find(name);
		if (it == _global_map.end()) {
			return "";
		}
		return it->second;
	}


	bool
	Data::isVector(const std::string& name) const
	{
		const std::string  prefix   = "std::vector<";
		const std::string  postfix  = ">";
		const std::string& typeName = getType(name);
		if ((not typeName.compare(0, prefix.size(), prefix))
		    and not typeName.compare(typeName.length() - postfix.length(), postfix.length(), postfix)) {
			return true;
		}
		return false;
	}


	std::string
	Data::getVariableInsertionErrorMsg(const std::vector<std::string>& quantityNames,
	                                   const std::string&              quantityName)
	{
		std::stringstream msgStream;
		if (quantityNames.size() > 1) {
			msgStream << "Could not insert variable \"" << quantityName << "\" when registering calculation for quantities \"[";
			for (size_t i = 0; i < quantityNames.size() - 1; ++i) {
				msgStream << quantityNames[i] << ", ";
			}
			msgStream << quantityNames[quantityNames.size() - 1] << "]\" (double entry?)." << std::endl;
		} else {
			msgStream << getVariableInsertionErrorMsg(quantityNames[0]);
		}
		return msgStream.str();
	}


	std::string
	Data::getVariableInsertionErrorMsg(const std::string& variableName)
	{
		std::stringstream msgStream;
		msgStream << "Could not insert variable \"" << variableName << "\" (double entry?)." << std::endl;
		return msgStream.str();
	}

}  // antok namespace