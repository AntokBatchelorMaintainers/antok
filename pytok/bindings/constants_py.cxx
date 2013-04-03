#include<constants_py.h>

namespace bp = boost::python;

void antok::py::exportConstants() {

	bp::class_<antok::Constants>("Constants")

		.add_static_property("chargedKaonMass"
		                     , bp::make_function(&antok::Constants::chargedKaonMass
		                                         , bp::return_value_policy<bp::return_by_value>()))
		.add_static_property("chargedPionMass"
		                     , bp::make_function(&antok::Constants::chargedPionMass
		                                         , bp::return_value_policy<bp::return_by_value>()))

		.add_static_property("protonMass"
		                     , bp::make_function(&antok::Constants::protonMass
		                                         , bp::return_value_policy<bp::return_by_value>()))
		.add_static_property("nParticles"
		                     , bp::make_function(&antok::Constants::nParticles
		                                         , bp::return_value_policy<bp::return_by_value>()))
		.add_static_property("initialized"
		                     , bp::make_function(&antok::Constants::initialized
		                                         , bp::return_value_policy<bp::return_by_value>()));

}
