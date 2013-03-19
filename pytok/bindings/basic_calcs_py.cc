#include<basic_calcs_py.h>

#include<TVector3.h>
#include<TLorentzVector.h>

#include<rootConverters_py.h>

namespace bp = boost::python;

namespace {

	PyObject* __get_beam_energy(PyObject* pyBeam, PyObject* pyLV_X) {
		TVector3* beam = antok::py::convertFromPy<TVector3*>(pyBeam);
		TLorentzVector* LV_X = antok::py::convertFromPy<TLorentzVector*>(pyLV_X);
		return antok::py::convertToPy<TLorentzVector>(antok::get_beam_energy(*beam, *LV_X));
	}

	bp::tuple __get_RPD_delta_phi_res_projection(PyObject* pyBeam,
	                                             PyObject* pyProton,
	                                             PyObject* pyX)
	{
		const TLorentzVector* beam = antok::py::convertFromPy<TLorentzVector*>(pyBeam);
		const TLorentzVector* proton = antok::py::convertFromPy<TLorentzVector*>(pyProton);
		const TLorentzVector* X = antok::py::convertFromPy<TLorentzVector*>(pyX);
		double deltaPhi, res;
		antok::get_RPD_delta_phi_res_projection(*beam, *proton, *X, deltaPhi, res);
		return bp::make_tuple(deltaPhi, res);
	}

	bp::tuple __get_RPD_delta_phi_res_rotation(PyObject* pyBeam,
	                                           PyObject* pyProton,
	                                           PyObject* pyX)
	{
		const TLorentzVector* beam = antok::py::convertFromPy<TLorentzVector*>(pyBeam);
		const TLorentzVector* proton = antok::py::convertFromPy<TLorentzVector*>(pyProton);
		const TLorentzVector* X = antok::py::convertFromPy<TLorentzVector*>(pyX);
		double deltaPhi, res;
		antok::get_RPD_delta_phi_res_rotation(*beam, *proton, *X, deltaPhi, res);
		return bp::make_tuple(deltaPhi, res);
	}


}

void antok::py::exportBasicCalcs() {

	bp::def("get_beam_energy", &__get_beam_energy);
	bp::def("get_RPD_delta_phi_res_projection", &__get_RPD_delta_phi_res_projection);
	bp::def("get_RPD_delta_phi_res_rotation", &__get_RPD_delta_phi_res_rotation);

}
