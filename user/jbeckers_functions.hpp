#ifndef ANTOK_USER_JBECKERS_FUNCTIONS_HPP
#define ANTOK_USER_JBECKERS_FUNCTIONS_HPP

#include "TVector3.h"
#include "TLorentzVector.h"

namespace antok {
namespace user {
namespace jbeckers {
namespace functions {

	// Multiplies a variable by a constant scale (must be of type double):
	template <typename T>
	class Scale: public Function
	{

	public:

		Scale(const double& scaleFactor,
		      const T&      value,
		      T&            result)
			: _scaleFactor(scaleFactor),
			  _value      (value),
			  _result     (result)
		{ }

		virtual ~Scale() { }

		bool
		operator() ()
		{
			_result = _scaleFactor * _value;
			return true;
		}

	private:

		const double _scaleFactor;  // constant parameter, needs to be copied (-> no & !)
		const T&     _value;
		T&           _result;

	};


	template <>
	class Scale<std::vector<double>>: public Function
	{

	public:

		Scale(const double&               scaleFactor,
		      const std::vector<double>&  values,
		      std::vector<double>&        results)
			: _scaleFactor(scaleFactor),
			  _values     (values),
			  _results    (results)
		{ }

		virtual ~Scale() { }

		bool
		operator() ()
		{
			const size_t sizeVec = _values.size();
			_results.resize(sizeVec);
			for (size_t i = 0; i < sizeVec; ++i) {
				_results[i] = _scaleFactor * _values[i];
			}
			return true;
		}

	private:

		const double               _scaleFactor;  // constant parameter, needs to be copied (-> no & !)
		const std::vector<double>& _values;
		std::vector<double>&       _results;

	};

}
}
}
}

#endif  // ANTOK_USER_JBECKERS_FUNCTIONS_HPP