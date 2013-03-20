#include<initializer_py.h>

namespace bp = boost::python;

void antok::py::exportInitializer() {

	bp::class_<antok::Initializer, boost::noncopyable>("Initializer", bp::no_init)
		.add_static_property("instance"
		                     , bp::make_function(&antok::Initializer::instance
	                         , bp::return_value_policy<bp::reference_existing_object>()))
		.def("readConfigFile", &antok::Initializer::readConfigFile)

		.def("initAll", &antok::Initializer::initAll)

		.def("initializeCutter", &antok::Initializer::initializeCutter)
		.def("initializeData", &antok::Initializer::initializeData)
		.def("initializeEvent", &antok::Initializer::initializeEvent)
		.def("initializePlotter", &antok::Initializer::initializePlotter);

}
