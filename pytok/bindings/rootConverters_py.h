#ifndef ANTOK_ROOTCONVERTERS_PY_H
#define ANTOK_ROOTCONVERTERS_PY_H

#include<boost/python.hpp>

class TLorentzRotation;
class TVector3;

namespace antok {

	namespace py {

		template<typename T>
		PyObject* convertToPy(const T& cxxObj);

		template<typename T>
		T convertFromPy(PyObject* pyObj);

		void exportRootConverters();

	}

}

#endif
