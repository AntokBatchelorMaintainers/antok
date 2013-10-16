#include<basic_calcs_py.h>

#include<TVector3.h>
#include<TLorentzVector.h>

#include<rootConverters_py.h>

namespace bp = boost::python;

namespace {

	PyObject* __getBeamEnergy(PyObject* pyBeam, PyObject* pyLV_X) {
		TVector3* beam = antok::py::convertFromPy<TVector3*>(pyBeam);
		TLorentzVector* LV_X = antok::py::convertFromPy<TLorentzVector*>(pyLV_X);
		return antok::py::convertToPy<TLorentzVector>(antok::getBeamEnergy(*beam, *LV_X));
	}

	bp::tuple __getBoostToCenterOfMassSystem(PyObject* pyBeam)
	{
		const TLorentzVector* beam = antok::py::convertFromPy<TLorentzVector*>(pyBeam);
		double centerOfMassEnergy;
		TVector3 boostVector;
		antok::getBoostToCenterOfMassSystem(*beam, centerOfMassEnergy, boostVector);
		PyObject* pyBoostVector = antok::py::convertToPy<TVector3>(boostVector);
		return bp::make_tuple(centerOfMassEnergy, *pyBoostVector);
	}

}

void antok::py::exportBasicCalcs() {

	bp::def("getBeamEnergy", &__getBeamEnergy);
	bp::def("getBoostToCenterOfMassSystem", &__getBoostToCenterOfMassSystem);
	bp::def("getPositiveSolutionOfQuadraticEquation", &antok::utils::getPositiveSolutionOfQuadraticEquation);

}
