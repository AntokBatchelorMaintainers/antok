#include<data.h>

#include<iostream>
#include<sstream>

#include<TLorentzVector.h>

bool antok::Data::insertDouble(std::string name) {
	if(doubles.count(name) > 0) {
		return false;
	}
	global_map[name] = "double";
	doubles[name] = -8888.8;
	return true;
}

bool antok::Data::insertInt(std::string name) {
	if(ints.count(name) > 0) {
		return false;
	}
	global_map[name] = "int";
	ints[name] = -8888;
	return true;
}

bool antok::Data::insertLong64_t(std::string name) {
	if(long64_ts.count(name) > 0) {
		return false;
	}
	global_map[name] = "Long64_t";
	long64_ts[name] = -8888;
	return true;
}

bool antok::Data::insertLorentzVector(std::string name) {
	if(lorentzVectors.count(name) > 0) {
		return false;
	}
	global_map[name] = "TLorentzVector";
	lorentzVectors[name] = *(new TLorentzVector);
	return true;
}

double* antok::Data::getDoubleAddr(std::string name) {
	if(doubles.count(name) < 1) {
		return 0;
	}
	return &doubles[name];
}

int* antok::Data::getIntAddr(std::string name) {
	if(ints.count(name) < 1) {
		return 0;
	}
	return &ints[name];
}

Long64_t* antok::Data::getLong64_tAddr(std::string name) {
	if(long64_ts.count(name) < 1) {
		return 0;
	}
	return &long64_ts[name];
}

TLorentzVector* antok::Data::getLorentzVectorAddr(std::string name) {
	if(lorentzVectors.count(name) < 1) {
		return 0;
	}
	return &lorentzVectors[name];
}

std::string antok::Data::getType(std::string name) {
	if(global_map.count(name) > 0) {
		return global_map[name];
	}
	return "";
}

std::string antok::Data::getVariableInsertionErrorMsg(std::vector<std::string> quantityNames,
                                                      std::string quantityName)
{
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

std::string antok::Data::getVariableInsertionErrorMsg(std::string variableName) {

	std::stringstream msgStream;
	msgStream<<"Could not insert variable \""<<variableName<<"\" (double entry?)."<<std::endl;
	return msgStream.str();

};

