#ifndef ANTOK_YAML_UTILS_HPP
#define ANTOK_YAML_UTILS_HPP

#include<iostream>

#include<yaml-cpp/yaml.h>

#include<data.h>
#include<object_manager.h>

namespace antok {

	namespace YAMLUtils {

		inline std::string getString(const YAML::Node& node) {
			try{
				return node.as<std::string>();
			} catch(const YAML::TypedBadConversion<std::string>& e) {
				return "";
			}
		}

		template<typename T>
		inline bool getValue(const YAML::Node& node, T* valPtr) {
			try{
				*valPtr = node.as<T>();
				return true;
			} catch(const YAML::TypedBadConversion<T>& e) {
				return false;
			}
		}

		template<typename T>
		inline T* getAddress(const YAML::Node& node) {
			T* retval = 0;
			try {
				T val = node.as<T>();
				retval = new T(val);
			} catch (const YAML::TypedBadConversion<T>& e) {
				// not bad yet, could be a variable name there
			}
			if(retval == 0) {
				antok::Data& data = antok::ObjectManager::instance()->getData();
				std::string name = antok::YAMLUtils::getString(node);
				if(name == "") {
					std::cerr<<"Entry has to be either a variable name or a convertible type."<<std::endl;
					return 0;
				}
				retval = data.getAddr<T>(name);
				if(retval == 0) {
					std::cerr<<"Variable \""<<name<<"\" not found in Data."<<std::endl;
					return 0;
				}
			}
			return retval;
		}

		template<>
		inline TLorentzVector* getAddress<TLorentzVector>(const YAML::Node& node) {
			TLorentzVector* retval = 0;
			antok::Data& data = antok::ObjectManager::instance()->getData();
			std::string name = antok::YAMLUtils::getString(node);
			if(name == "") {
				std::cerr<<"Entry has to be either a variable name or a convertible type."<<std::endl;
				return 0;
			}
			retval = data.getAddr<TLorentzVector>(name);
			if(retval == 0) {
				std::cerr<<"Variable \""<<name<<"\" not found in Data."<<std::endl;
				return 0;
			}
			return retval;
		}

		inline bool hasNodeKey(const YAML::Node& node, std::string key) {
			try {
				return node[key];
			} catch(const YAML::BadSubscript&) {
				return false;
			}
		}
		inline bool handleOnOffOption(std::string optionName, const YAML::Node& option, std::string location) {

			using antok::YAMLUtils::hasNodeKey;

			if (not hasNodeKey(option, optionName)) {
				std::cerr << "Warning: \"" << optionName << "\" not found in \"" << location << "\", switching it off" << std::endl;
			} else {
				std::string optionValue = antok::YAMLUtils::getString(option[optionName]);
				if (optionValue == "On") {
					return true;
				} else if (optionValue == "Off") {
					// returning at end of function
				} else {
					std::cerr << "Warning: \"" << location << "\"'s \"" << optionName << "\" is \"" << optionValue
					        << "\" instead of \"On\" or \"Off\", switching it off" << std::endl;
				}
			}
			return false;

		}

	}

}

#endif

