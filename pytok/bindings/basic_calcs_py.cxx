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

	bp::tuple __getRPDDeltaPhiResProjection(PyObject* pyBeam,
	                                             PyObject* pyProton,
	                                             PyObject* pyX)
	{
		const TLorentzVector* beam = antok::py::convertFromPy<TLorentzVector*>(pyBeam);
		const TLorentzVector* proton = antok::py::convertFromPy<TLorentzVector*>(pyProton);
		const TLorentzVector* X = antok::py::convertFromPy<TLorentzVector*>(pyX);
		double deltaPhi, res;
		antok::getRPDDeltaPhiResProjection(*beam, *proton, *X, deltaPhi, res);
		return bp::make_tuple(deltaPhi, res);
	}

	bp::tuple __getRPDDeltaPhiResRotation(PyObject* pyBeam,
	                                           PyObject* pyProton,
	                                           PyObject* pyX)
	{
		const TLorentzVector* beam = antok::py::convertFromPy<TLorentzVector*>(pyBeam);
		const TLorentzVector* proton = antok::py::convertFromPy<TLorentzVector*>(pyProton);
		const TLorentzVector* X = antok::py::convertFromPy<TLorentzVector*>(pyX);
		double deltaPhi, res;
		antok::getRPDDeltaPhiResRotation(*beam, *proton, *X, deltaPhi, res);
		return bp::make_tuple(deltaPhi, res);
	}


}

void antok::py::exportBasicCalcs() {

	bp::def("getBeamEnergy", &__getBeamEnergy);
	bp::def("getRPDDeltaPhiResProjection", &__getRPDDeltaPhiResProjection);
	bp::def("getRPDDeltaPhiResRotation", &__getRPDDeltaPhiResRotation);

}
