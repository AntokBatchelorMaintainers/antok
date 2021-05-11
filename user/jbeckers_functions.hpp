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

		Scale(const double& scale,
		      const T&      scaleOn,
		      T&            out)
			  : _scale(scale),
			    _scaleOn(scaleOn),
			    _out(out)
		{ }

		virtual ~Scale() { }

		bool
		operator() ()
		{
			_out = _scale * _scaleOn;
			return true;
		}

	private:

		const double _scale; // constant, needs to be copied (-> no & !)
		const T&     _scaleOn;
		T&           _out;

	};

	template <>
	class Scale<std::vector<double>>: public Function
	{

	public:

		Scale(const double&               scale,
		      const std::vector<double>&  scaleOn,
		      std::vector<double>&        out)
			  : _scale(scale),
			    _scaleOn(scaleOn),
			    _out(out)
		{ }

		virtual ~Scale() { }

		bool
		operator() ()
		{
			const size_t sizeVec = _scaleOn.size();
			_out.resize(sizeVec);
			for (size_t i = 0; i < sizeVec; ++i) {
				_out[i] = _scale * _scaleOn[i];
			}
			return true;
		}

	private:

		const double                   _scale; // constant, needs to be copied (-> no & !)
		const std::vector<double>&     _scaleOn;
		std::vector<double>&           _out;

	};

}
}
}
}

#endif