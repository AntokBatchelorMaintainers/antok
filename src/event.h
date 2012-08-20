#ifndef HLIB_EVENT_H
#define HLIB_EVENT_H

#include<TLorentzVector.h>

#include<data.hpp>

namespace hlib {

	class Event {

	  public:

		Event() { };

		void update(const hlib::Data& data);

		const TLorentzVector& get_p1() const { return p1; };
		const TLorentzVector& get_p2() const { return p2; };
		const TLorentzVector& get_p3() const { return p3; };
		const TLorentzVector& get_p4() const { return p4; };
		const TLorentzVector& get_p5() const { return p5; };

		const TLorentzVector& get_pSum() const { return pSum; };
		const TLorentzVector& get_pProton() const { return pProton; };

		const TVector3& get_p3Beam() const { return p3Beam; };
		const TLorentzVector& get_pBeam() const { return pBeam; };

		const double& get_t() const { return t; };
		const double& get_tMin() const { return tMin; };
		const double& get_tPrime() const { return tPrime; };

		const double& get_RpdDeltaPhi() const { return RpdDeltaPhi; };
		const double& get_RpdPhiRes() const { return RpdPhiRes; };

	  private:

		TLorentzVector get_beam_energy(TVector3 p3_beam, const TLorentzVector& LV_X) const;

		TLorentzVector p1;
		TLorentzVector p2;
		TLorentzVector p3;
		TLorentzVector p4;
		TLorentzVector p5;

		TLorentzVector pSum;

		TLorentzVector pProton;

		TVector3 p3Beam;
		TLorentzVector pBeam;

		double t;
		double tMin;
		double tPrime;

		double RpdDeltaPhi;
		double RpdPhiRes;

		const hlib::Data* rawData;

	};

}

#endif
