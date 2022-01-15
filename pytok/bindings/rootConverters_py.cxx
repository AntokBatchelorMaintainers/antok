#include "rootConverters_py.h"

#include<TLorentzVector.h>
#include<TPython.h>
#include<TVector3.h>

namespace bp = boost::python;

template<typename T>
PyObject* antok::py::convertToPy(const T& cxxObj) {
	T* newCxxObj = new T(cxxObj);
	return TPython::CPPInstance_FromVoidPtr(newCxxObj, newCxxObj->ClassName(), true);
};

template<typename T>
T antok::py::convertFromPy(PyObject* pyObj) {
	TObject* TObj = (TObject*)(TPython::CPPInstance_AsVoidPtr(pyObj));
	T cxxObj = dynamic_cast<T>(TObj);
	return cxxObj;
};

void antok::py::exportRootConverters() {

	bp::def("__RootConverters_convertToPy_TVector3", &antok::py::convertToPy<TVector3>);
	bp::def(
		"__RootConverters_convertFromPy_TVector3", &antok::py::convertFromPy<TVector3*>
		, bp::return_internal_reference<1>()
	);

	bp::def("__RootConverters_convertToPy_TLorentzVector", &antok::py::convertToPy<TLorentzVector>);
	bp::def(
		"__RootConverters_convertFromPy_TLorentzVector", &antok::py::convertFromPy<TLorentzVector*>
		, bp::return_internal_reference<1>()
	);

}
