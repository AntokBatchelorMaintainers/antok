#include<basic_calcs.h>
#include<constants.h>
#include<event.h>

using hlib::PION_MASS;

void hlib::Event::update(const hlib::Data& data) {

	p1.SetXYZM(data.Mom_x1, data.Mom_y1, data.Mom_z1, PION_MASS);
	p2.SetXYZM(data.Mom_x2, data.Mom_y2, data.Mom_z2, PION_MASS);
	p3.SetXYZM(data.Mom_x3, data.Mom_y3, data.Mom_z3, PION_MASS);
	p4.SetXYZM(data.Mom_x4, data.Mom_y4, data.Mom_z4, PION_MASS);
	p5.SetXYZM(data.Mom_x5, data.Mom_y5, data.Mom_z5, PION_MASS);

	pSum = p1+p2+p3+p4+p5;

	pProton.SetPxPyPzE(data.RPD_Px, data.RPD_Py, data.RPD_Pz, data.RPD_E);

	p3Beam.SetXYZ(data.gradx, data.grady, 1.);
	pBeam = hlib::get_beam_energy(p3Beam, pSum);
	p3Beam = pBeam.Vect();

	t = std::fabs((pBeam - pSum).Mag2());
	tMin = std::fabs((std::pow(pSum.M2() - pBeam.M2(), 2)) / (4. * p3Beam.Mag2()));
	tPrime = t - tMin;

	hlib::get_RPD_delta_phi_res(pBeam, pProton, pSum, RpdDeltaPhi, RpdPhiRes);

};

