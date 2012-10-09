#ifndef ANTOK_EVENT_H
#define ANTOK_EVENT_H

#include<vector>

#include<TLorentzVector.h>

#include<data.hpp>

namespace antok {

	class Function;

	class Event {

	  public:

		static Event* instance();

		void update(antok::Data& data);

		const TLorentzVector& get_p(unsigned int n) const { return _p.at(n); };

		const TLorentzVector& get_pSum() const { return _pSum; };
		const TLorentzVector& get_pProton() const { return _pProton; };

		const TVector3& get_p3Beam() const { return _p3Beam; };
		const TLorentzVector& get_pBeam() const { return _pBeam; };

		const double& get_t() const { return _t; };
		const double& get_tMin() const { return _tMin; };
		const double& get_tPrime() const { return _tPrime; };

		const double& get_RpdDeltaPhi() const { return _RpdDeltaPhi; };
		const double& get_RpdPhiRes() const { return _RpdPhiRes; };

		const double& get_RpdDeltaPhi_fhaas() const { return _RpdDeltaPhi_fhaas; };
		const double& get_RpdPhiRes_fhaas() const { return _RpdPhiRes_fhaas; };

		const antok::Data* rawData;

	  private:

		Event();

		static Event* _event;

		std::vector<antok::Function*> _functions;

		friend class antok::Initializer;

//		TLorentzVector get_beam_energy(TVector3 p3_beam, const TLorentzVector& LV_X) const;

		std::vector<TLorentzVector> _p;

		TLorentzVector _pSum;

		TLorentzVector _pProton;

		TVector3 _p3Beam;
		TLorentzVector _pBeam;

		double _t;
		double _tMin;
		double _tPrime;

		double _RpdDeltaPhi;
		double _RpdPhiRes;

		double _RpdDeltaPhi_fhaas;
		double _RpdPhiRes_fhaas;

	};

}

#endif
