#include<constants_py.h>

namespace bp = boost::python;

void antok::py::exportConstants() {

	bp::class_<antok::Constants>("Constants")

		.add_static_property("charged_kaon_mass"
		                     , bp::make_function(&antok::Constants::charged_kaon_mass
		                                         , bp::return_value_policy<bp::return_by_value>()))
		.add_static_property("charged_pion_mass"
		                     , bp::make_function(&antok::Constants::charged_pion_mass
		                                         , bp::return_value_policy<bp::return_by_value>()))

		.add_static_property("proton_mass"
		                     , bp::make_function(&antok::Constants::proton_mass
		                                         , bp::return_value_policy<bp::return_by_value>()))
		.add_static_property("n_particles"
		                     , bp::make_function(&antok::Constants::n_particles
		                                         , bp::return_value_policy<bp::return_by_value>()))
		.add_static_property("initialized"
		                     , bp::make_function(&antok::Constants::initialized
		                                         , bp::return_value_policy<bp::return_by_value>()));

}
