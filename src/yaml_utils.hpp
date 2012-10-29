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
			} catch(YAML::TypedBadConversion<std::string> e) {
				return "";
			}
		}

		template<typename T>
		inline bool getValue(const YAML::Node& node, T* valPtr) {
			try{
				*valPtr = node.as<T>();
				return true;
			} catch(YAML::TypedBadConversion<T> e) {
				return false;
			}
		}

		template<typename T>
		inline T* getAddress(const YAML::Node& node) {
			T* retval = 0;
			try {
				T val = node.as<T>();
				retval = new T(val);
			} catch (YAML::TypedBadConversion<T> e) {
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
			} catch(YAML::BadSubscript) {
				return false;
			}
		}

	}

}

#endif

