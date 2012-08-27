#include<basic_calcs.h>
#include<constants.hpp>
#include<event.h>

using hlib::PION_MASS;

hlib::Event* hlib::Event::_event = NULL;

hlib::Event* hlib::Event::instance() {
	if(_event == NULL) {
		_event = new hlib::Event();
	}
	return _event;
}


hlib::Event::Event() {

	_p.resize(hlib::N_PARTICLES);

};

void hlib::Event::update(const hlib::Data& data) {

	rawData = &data;

	_pSum.SetXYZT(0., 0., 0., 0.);
	for(unsigned int i = 0; i < _p.size(); ++i) {
		_p.at(i).SetXYZM(data.Mom_x.at(i), data.Mom_y.at(i), data.Mom_z.at(i), PION_MASS);
		_pSum += _p.at(i);
	}

	_pProton.SetPxPyPzE(data.RPD_Px, data.RPD_Py, data.RPD_Pz, data.RPD_E);

	_p3Beam.SetXYZ(data.gradx, data.grady, 1.);
	_pBeam = hlib::get_beam_energy(_p3Beam, _pSum);
	_p3Beam = _pBeam.Vect();

	_t = std::fabs((_pBeam - _pSum).Mag2());
	_tMin = std::fabs((std::pow(_pSum.M2() - _pBeam.M2(), 2)) / (4. * _p3Beam.Mag2()));
	_tPrime = _t - _tMin;

	hlib::get_RPD_delta_phi_res(_pBeam, _pProton, _pSum, _RpdDeltaPhi, _RpdPhiRes);

};

