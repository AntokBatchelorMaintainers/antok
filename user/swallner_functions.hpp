/*
 * swallner_functions.hpp
 *
 *  Created on: Dec 22, 2015
 *      Author: ga34liw
 */

#ifndef USER_SWALLNER_FUNCTIONS_HPP_
#define USER_SWALLNER_FUNCTIONS_HPP_

#include <cmath>
#include "swallner.h"


#include "TVector3.h"

namespace antok{
namespace user{
namespace stefan{
namespace functions{

	class CalcLTProjection: public Function {
	public:

		CalcLTProjection(TVector3 const* vector, TVector3 const* direction,
				double* projection_longitudinal, double* projection_transversal):
			vector_(vector),
			direction_(direction),
			projection_longitudinal_( projection_longitudinal ),
			projection_transversal_( projection_transversal ){}


		bool operator() (){
			TVector3 const normalized_direction = direction_->Unit();
			*projection_longitudinal_ = fabs( normalized_direction.Dot(*vector_) );
			*projection_transversal_ = sqrt( vector_->Mag2() - *projection_longitudinal_ * *projection_longitudinal_ );
			return true;
		}



	private:
		TVector3 const * const vector_;
		TVector3 const * const direction_;
		double* const projection_longitudinal_;
		double* const projection_transversal_;

	};

	class CalcArmenterosAlpha: public Function {
	public:

		CalcArmenterosAlpha(double const* longitudinal_mom_1, double const* longitudinal_mom_2,
							double* alpha ):
								longitudinal_mom_1_( longitudinal_mom_1 ),
								longitudinal_mom_2_( longitudinal_mom_2 ),
								alpha_( alpha )
								{}



		bool operator() (){
				*alpha_ = (*longitudinal_mom_1_ - *longitudinal_mom_2_ ) / (*longitudinal_mom_1_ + *longitudinal_mom_2_);
			return true;
		}



	private:
		double const* const longitudinal_mom_1_;
		double const* const longitudinal_mom_2_;
		double* const alpha_;

	};

	class CalcRICHPID: public Function {
	public:

		CalcRICHPID( double const* L_pion,
							 double const* L_kaon,
							 double const* L_proton,
							 double const* L_electron,
							 double const* L_muon,
							 double const* L_background,
		                     double* P_pion,
							 double* P_kaon,
							 double* P_proton,
							 double* P_electron,
							 double* P_muon,
							 double* P_background
							 ):
		                     L_pion_(         *L_pion),
							 L_kaon_(         *L_kaon ),
							 L_proton_(       *L_proton ),
							 L_electron_(     *L_electron ),
							 L_muon_(     *L_muon ),
							 L_background_(   *L_background ),
		                     P_pion_(         *P_pion ),
							 P_kaon_(         *P_kaon ),
							 P_proton_(       *P_proton ),
							 P_electron_(     *P_electron ),
							 P_muon_(     *P_muon ),
							 P_background_ (  *P_background )
		{}



		bool operator() (){
			P_pion_       = ( L_pion_       >= 0.0 )?  L_pion_ / fmax( L_kaon_, fmax(L_proton_, fmax( L_electron_, fmax(L_muon_, L_background_)))) : -1.0;
			P_kaon_       = ( L_kaon_       >= 0.0 )?  L_kaon_ / fmax( L_pion_, fmax(L_proton_, fmax( L_electron_, fmax(L_muon_, L_background_)))) : -1.0;
			P_proton_     = ( L_proton_     >= 0.0 )?  L_proton_ / fmax( L_kaon_, fmax(L_pion_, fmax( L_electron_, fmax(L_muon_, L_background_)))) : -1.0;
			P_electron_   = ( L_electron_   >= 0.0 )?  L_electron_ / fmax( L_kaon_, fmax(L_proton_, fmax( L_pion_, fmax(L_muon_, L_background_)))) : -1.0;
			P_muon_       = ( L_muon_       >= 0.0 )?  L_muon_/ fmax( L_kaon_, fmax(L_proton_, fmax( L_pion_, fmax(L_pion_ , L_background_))))     : -1.0;
			P_background_ = ( L_background_ >= 0.0 )?  L_background_/ fmax( L_kaon_, fmax(L_proton_, fmax( L_electron_, fmax(L_muon_, L_pion_))))  : -1.0;
			return true;
		}



	private:
		double const& L_pion_;
		double const& L_kaon_;
		double const& L_proton_;
		double const& L_electron_;
		double const& L_muon_;
		double const& L_background_;
		double& P_pion_;
		double& P_kaon_;
		double& P_proton_;
		double& P_electron_;
		double& P_muon_;
		double& P_background_;

	};
}
}
}
}




#endif /* USER_SWALLNER_FUNCTIONS_HPP_ */
