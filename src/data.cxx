#include<data.h>

#include<iostream>
#include<sstream>

#include<TLorentzVector.h>

template<typename T>
bool antok::Data::insert(std::string name) {
	std::cerr<<"Could not insert variable of this type."<<std::endl;
	return false;
}

template<typename T>
T* antok::Data::getAddr(std::string name) {
	std::cerr<<"Could not get address for variable of this type."<<std::endl;
	return 0;
}

namespace antok {

	template<>
	bool antok::Data::insert<double>(std::string name) {
		if(doubles.count(name) > 0) {
			return false;
		}
		global_map[name] = "double";
		doubles[name] = -8888.8;
		return true;
	}

	template<>
	bool antok::Data::insert<int>(std::string name) {
		if(ints.count(name) > 0) {
			return false;
		}
	   global_map[name] = "int";
	   ints[name] = -8888;
	   return true;
	}

	template<>
	bool antok::Data::insert<Long64_t>(std::string name) {
		if(long64_ts.count(name) > 0) {
			return false;
		}
		global_map[name] = "Long64_t";
		long64_ts[name] = -8888;
		return true;
	}

    template<>
    bool antok::Data::insert<std::vector<int> >(std::string name) {
        if(intVectors.count(name) > 0) {
            return false;
        }
        global_map[name] = "std::vector<int>";
        intVectors[name] = (new std::vector<int>);
        return true;
    }

    template<>
	bool antok::Data::insert<std::vector<double> >(std::string name) {
		if(doubleVectors.count(name) > 0) {
			return false;
		}
		global_map[name] = "std::vector<double>";
		doubleVectors[name] = (new std::vector<double>);
		return true;
	}

	template<>
	bool antok::Data::insert<TLorentzVector>(std::string name) {
		if(lorentzVectors.count(name) > 0) {
			return false;
		}
		global_map[name] = "TLorentzVector";
		lorentzVectors[name] = (new TLorentzVector);
		return true;
	}

	template<>
	bool antok::Data::insert<std::vector<TLorentzVector> >(std::string name) {
		if(lorentzVectorVectors.count(name) > 0) {
			return false;
		}
		global_map[name] = "std::vector<TLorentzVector>";
		lorentzVectorVectors[name] = (new std::vector<TLorentzVector>);
		return true;
	}

	template<>
	bool antok::Data::insert<TVector3>(std::string name) {
		if(lorentzVectors.count(name) > 0) {
			return false;
		}
		global_map[name] = "TVector3";
		vectors[name] = *(new TVector3);
		return true;
	}

	template<>
	bool antok::Data::insert<std::vector<TVector3> >(std::string name) {
		if(vectorVectors.count(name) > 0) {
			return false;
		}
		global_map[name] = "std::vector<TVector3>";
		vectorVectors[name] = (new std::vector<TVector3>);
		return true;
	}

	template<>
	double* antok::Data::getAddr<double>(std::string name) {
		if(doubles.count(name) < 1) {
			return 0;
		}
		return &doubles[name];
	}

	template<>
	int* antok::Data::getAddr<int>(std::string name) {
		if(ints.count(name) < 1) {
			return 0;
		}
		return &ints[name];
	}

	template<>
	Long64_t* antok::Data::getAddr<Long64_t>(std::string name) {
		if(long64_ts.count(name) < 1) {
			return 0;
		}
		return &long64_ts[name];
	}

    template<>
    std::vector<int>* antok::Data::getAddr<std::vector<int> >(std::string name) {
        if(intVectors.count(name) < 1) {
            return 0;
        }
        return intVectors[name];
    }

	template<>
	std::vector<double>* antok::Data::getAddr<std::vector<double> >(std::string name) {
		if(doubleVectors.count(name) < 1) {
			return 0;
		}
		return doubleVectors[name];
	}

	template<>
	TLorentzVector* antok::Data::getAddr<TLorentzVector>(std::string name) {
		if(lorentzVectors.count(name) < 1) {
			return 0;
		}
		return lorentzVectors[name];
	}

	template<>
	std::vector<TLorentzVector>* antok::Data::getAddr<std::vector<TLorentzVector> >(std::string name) {
		if(lorentzVectorVectors.count(name) < 1) {
			return 0;
		}
		return lorentzVectorVectors[name];
	}

	template<>
	TVector3* antok::Data::getAddr<TVector3>(std::string name) {
		if(vectors.count(name) < 1) {
			return 0;
		}
		return &vectors[name];
	}

	template<>
	std::vector<TVector3>* antok::Data::getAddr<std::vector<TVector3> >(std::string name) {
		if(vectorVectors.count(name) < 1) {
			return 0;
		}
		return vectorVectors[name];
	}

}

std::string antok::Data::getType(std::string name) {
	if(global_map.count(name) > 0) {
		return global_map[name];
	}
	return "";
}

bool antok::Data::isVector(std::string name) {
	std::string prefix = "std::vector<";
	std::string postfix = ">";
	const std::string& typeName = getType(name);
	if((not typeName.compare(0, prefix.size(), prefix)) and
	    not typeName.compare(typeName.length() - postfix.length(), postfix.length(), postfix)) {
		return true;
	}
	return false;
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
		msgStream<<getVariableInsertionErrorMsg(quantityNames[0]);
	}
	return msgStream.str();
};

std::string antok::Data::getVariableInsertionErrorMsg(std::string variableName) {

	std::stringstream msgStream;
	msgStream<<"Could not insert variable \""<<variableName<<"\" (double entry?)."<<std::endl;
	return msgStream.str();

};

