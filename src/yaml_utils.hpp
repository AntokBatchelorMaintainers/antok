#ifndef ANTOK_YAML_UTILS_HPP
#define ANTOK_YAML_UTILS_HPP

#include <iostream>

#include "yaml-cpp/yaml.h"

#include "data.h"
#include "object_manager.h"

namespace antok {

	namespace YAMLUtils {

		inline
		std::string
		getString(const YAML::Node& node)
		{
			try {
				return node.as<std::string>();
			} catch (const YAML::TypedBadConversion<std::string>& e) {
				return "";
			}
		}


		template<typename T>
		inline
		bool
		getValue(const YAML::Node& node,
		         T&                val)
		{
			try {
				val = node.as<T>();
				return true;
			} catch (const YAML::TypedBadConversion<T>& e) {
				return false;
			}
		}


		template<typename T>
		inline
		T*
		getAddress(const YAML::Node& node)
		{
			antok::Data& data = antok::ObjectManager::instance()->getData();
			T* retval = nullptr;
			try {
				// use antok::Data memory management to get new adress
				T val = node.as<T>();
				size_t nmb = 0;
				while (!data.insert<T>("__YAMLUtilsAddress" + std::to_string(nmb))) {
					if (nmb > 10000) { // TODO arbitrary limit
						std::cerr << "Could not allocate memory for YAMLUtils::getAddress.";
					}
					nmb++;
				}
				retval = data.getAddr<T>("__YAMLUtilsAddress" + std::to_string(nmb));
				*retval = val;
			} catch (const YAML::TypedBadConversion<T>& e) {
				// not bad yet, could be a variable name there
			}
			if (retval == nullptr) {
				const std::string name = antok::YAMLUtils::getString(node);
				if (name == "") {
					std::cerr << "Entry has to be either a variable name or a convertible type." << std::endl;
					return nullptr;
				}
				retval = data.getAddr<T>(name);
				if (retval == nullptr) {
					std::cerr << "Variable '"<< name << "' not found in Data." << std::endl;
					return nullptr;
				}
			}
			return retval;
		}


		template<>
		inline
		TLorentzVector*
		getAddress<TLorentzVector>(const YAML::Node& node)
		{
			TLorentzVector* retval = nullptr;
			antok::Data& data = antok::ObjectManager::instance()->getData();
			const std::string name = antok::YAMLUtils::getString(node);
			if (name == "") {
				std::cerr << "Entry has to be either a variable name or a convertible type." << std::endl;
				return nullptr;
			}
			retval = data.getAddr<TLorentzVector>(name);
			if (retval == nullptr) {
				std::cerr << "Variable '" << name << "' not found in Data." << std::endl;
				return nullptr;
			}
			return retval;
		}


		inline
		bool
		hasNodeKey(const YAML::Node&  node,
		           const std::string& key)
		{
			try {
				const YAML::Node tmpNode = node[key];
				return true;
			} catch(const YAML::BadSubscript& e) {
				return false;
			}
		}


		inline
		bool
		handleOnOffOption(const std::string& optionName,
		                  const YAML::Node&  option,
		                  const std::string& location)
		{
			using antok::YAMLUtils::hasNodeKey;
			if (not hasNodeKey(option, optionName)) {
				std::cerr << "Warning: option '" << optionName << "' not found at '" << location << "'. Switching it off." << std::endl;
			} else {
				const std::string optionValue = antok::YAMLUtils::getString(option[optionName]);
				if (optionValue == "On") {
					return true;
				} else if (optionValue == "Off") {
					// returning at end of function
				} else {
					std::cerr << "Warning: the value of option '" << optionName << "' at '" << location << "' is "
					          << "'" << optionValue << "' instead of 'On' or 'Off'. Switching it off." << std::endl;
				}
			}
			return false;

		}

	}

}

#endif  // ANTOK_YAML_UTILS_HPP
