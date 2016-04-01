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

	class CalcRICHProbabilities: public Function {
	public:

		/***
		 * If pid should be determined,
		 * 	kaon: 0
		 * 	pion: 1
		 * 	proton: 2
		 * 	electron: 3
		 * 	muon: 4
		 * 	background: 5
		 */
		CalcRICHProbabilities( double const* L_pion,
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
			P_pion_       = ( L_pion_       > 0.0 )?  L_pion_ / fmax( L_kaon_, fmax(L_proton_, fmax( L_electron_, fmax(L_muon_, L_background_)))) : -1.0;
			P_kaon_       = ( L_kaon_       > 0.0 )?  L_kaon_ / fmax( L_pion_, fmax(L_proton_, fmax( L_electron_, fmax(L_muon_, L_background_)))) : -1.0;
			P_proton_     = ( L_proton_     > 0.0 )?  L_proton_ / fmax( L_kaon_, fmax(L_pion_, fmax( L_electron_, fmax(L_muon_, L_background_)))) : -1.0;
			P_electron_   = ( L_electron_   > 0.0 )?  L_electron_ / fmax( L_kaon_, fmax(L_proton_, fmax( L_pion_, fmax(L_muon_, L_background_)))) : -1.0;
			P_muon_       = ( L_muon_       > 0.0 )?  L_muon_/ fmax( L_kaon_, fmax(L_proton_, fmax( L_pion_, fmax(L_pion_ , L_background_))))     : -1.0;
			P_background_ = ( L_background_ > 0.0 )?  L_background_/ fmax( L_kaon_, fmax(L_proton_, fmax( L_electron_, fmax(L_muon_, L_pion_))))  : -1.0;
			return true;
		}



	protected:
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
	class CalcRICHPID: public CalcRICHProbabilities{
	public:

		/***
		 * If pid should be determined,
		 * 	kaon: 0
		 * 	pion: 1
		 * 	proton: 2
		 * 	electron: 3
		 * 	muon: 4
		 * 	background: 5
		 */
		CalcRICHPID(         double const* L_pion,
							 double const* L_kaon,
							 double const* L_proton,
							 double const* L_electron,
							 double const* L_muon,
							 double const* L_background,
							 TVector3 const* mom,
							 double const* P_ratio_cut,
							 double const* Mom_pion_min,
							 double const* Mom_pion_max,
							 double const* Mom_kaon_min,
							 double const* Mom_kaon_max,
							 double const* Mom_proton_min,
							 double const* Mom_proton_max,
							 double const* Mom_electron_min,
							 double const* Mom_electron_max,
							 double const* Mom_muon_min,
							 double const* Mom_muon_max,
		                     double* P_pion,
							 double* P_kaon,
							 double* P_proton,
							 double* P_electron,
							 double* P_muon,
							 double* P_background,
							 int* pid
							 ):
								 CalcRICHProbabilities(  L_pion,
							  L_kaon,
							  L_proton,
							  L_electron,
							  L_muon,
							  L_background,
		                     P_pion,
							 P_kaon,
							 P_proton,
							 P_electron,
							 P_muon,
							 P_background ),
							 mom_(*mom),
							 P_ratio_cut_(*P_ratio_cut),
							 Mom_pion_min_(*Mom_pion_min),
							 Mom_pion_max_(*Mom_pion_max),
							 Mom_kaon_min_(*Mom_kaon_min),
							 Mom_kaon_max_(*Mom_kaon_max),
							 Mom_proton_min_(*Mom_proton_min),
							 Mom_proton_max_(*Mom_proton_max),
							 Mom_electron_min_(*Mom_electron_min),
							 Mom_electron_max_(*Mom_electron_max),
							 Mom_muon_min_(*Mom_muon_min),
							 Mom_muon_max_(*Mom_muon_max),
							 pid_( *pid)
		{}



		bool operator() (){
			CalcRICHProbabilities::operator ()();
			const double mom = mom_.Mag();
			pid_ = -1;
			if(      P_pion_ > P_ratio_cut_ && mom > Mom_pion_min_ && mom < Mom_pion_max_)             pid_ = 0;
			else if( P_kaon_ > P_ratio_cut_ && mom > Mom_kaon_min_ && mom < Mom_kaon_max_)             pid_ = 1;
			else if( P_proton_ > P_ratio_cut_ && mom > Mom_proton_min_ && mom < Mom_proton_max_)       pid_ = 2;
			else if( P_electron_ > P_ratio_cut_ && mom > Mom_electron_min_ && mom < Mom_electron_max_) pid_ = 3;
			else if( P_muon_ > P_ratio_cut_ && mom > Mom_muon_min_ && mom < Mom_muon_max_)             pid_ = 4;
			else if( P_background_ > P_ratio_cut_ )                                                    pid_ = 5;
			return true;
		}



	private:
		TVector3 const& mom_;
		double const& P_ratio_cut_;
		double const& Mom_pion_min_;
		double const& Mom_pion_max_;
		double const& Mom_kaon_min_;
		double const& Mom_kaon_max_;
		double const& Mom_proton_min_;
		double const& Mom_proton_max_;
		double const& Mom_electron_min_;
		double const& Mom_electron_max_;
		double const& Mom_muon_min_;
		double const& Mom_muon_max_;
		int& pid_;

	};

	/***
	 * is_kp_pk = -2 -> no decision possible
	 * is_kp_pk = 1  -> no first particle is kaon, second particle is pion
	 * is_kp_pk = 2  -> no first particle is pion, second particle is kaon
	 */
	class DetermineKaonPionLV: public Function {
	public:

		DetermineKaonPionLV( TVector3* mom_1, int const* pid_1,
		                     TVector3* mom_2, int const* pid_2,
							 double const* mass_charged_kaon, double const* mass_charged_pion ,
							 TLorentzVector* kaon_lv, TLorentzVector* pion_lv, int* is_kp_pk, int* pid_kaon, int* pid_pion):
                                mom_1_(             *mom_1),
                                pid_1_(             *pid_1),
                                mom_2_(             *mom_2),
                                pid_2_(             *pid_2),
                                mass_charged_kaon_( *mass_charged_kaon),
                                mass_charged_pion_( *mass_charged_pion),
                                kaon_lv_(           *kaon_lv),
                                pion_lv_(           *pion_lv),
                                is_kp_pk_(          *is_kp_pk),
								pid_kaon_(          *pid_kaon),
								pid_pion_(          *pid_pion)
								{}



		bool operator() (){
			if(        ( pid_1_ == 1 && pid_2_ != 1) || ( pid_2_ == 0 && pid_1_ != 0 ) ){ // 1 = kaon, 2 = pion

				kaon_lv_ = TLorentzVector( mom_1_, sqrt( mass_charged_kaon_ * mass_charged_kaon_ + mom_1_.Mag2() ) );
				pion_lv_ = TLorentzVector( mom_2_, sqrt( mass_charged_pion_ * mass_charged_pion_ + mom_2_.Mag2() ) );
				is_kp_pk_ = 1;
				pid_kaon_ = pid_1_;
				pid_pion_ = pid_2_;
			} else if( ( pid_1_ == 0 && pid_2_ != 0) || ( pid_2_ == 1 && pid_1_ != 1 ) ){ // 1 = pion , 2 = kaon
				kaon_lv_ = TLorentzVector( mom_2_, sqrt( mass_charged_kaon_ * mass_charged_kaon_ + mom_2_.Mag2() ) );
				pion_lv_ = TLorentzVector( mom_1_, sqrt( mass_charged_pion_ * mass_charged_pion_ + mom_1_.Mag2() ) );
				is_kp_pk_ = 2;
				pid_kaon_ = pid_2_;
				pid_pion_ = pid_1_;
			} else {
				is_kp_pk_ = -1;
				kaon_lv_ = TLorentzVector( -4444.0, -4444.0, -4444.0, -4444.0 );
				pion_lv_ = TLorentzVector( -4444.0, -4444.0, -4444.0, -4444.0 );
			}
		}



	private:
		TVector3 const& mom_1_;
		int const& pid_1_;
		TVector3 const& mom_2_;
		int const& pid_2_;
		double const& mass_charged_kaon_;
		double const& mass_charged_pion_;
		TLorentzVector& kaon_lv_;
		TLorentzVector& pion_lv_;
		int& is_kp_pk_;
		int& pid_kaon_;
		int& pid_pion_;

	};
}
}
}
}




#endif /* USER_SWALLNER_FUNCTIONS_HPP_ */
