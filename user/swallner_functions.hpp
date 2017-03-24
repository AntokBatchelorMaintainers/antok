/*
 * swallner_functions.hpp
 *
 *  Created on: Dec 22, 2015
 *      Author: ga34liw
 */

#ifndef USER_SWALLNER_FUNCTIONS_HPP_
#define USER_SWALLNER_FUNCTIONS_HPP_

#include <cmath>
#include <limits>
#include <iostream>
//#include "swallner.h"


#include "TVector3.h"
#include "TLorentzVector.h"


namespace antok{
	class Function;
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
            const double L_max_nopion       = _getMaxL( L_kaon_, L_proton_, L_electron_, L_muon_,     L_background_ );
            const double L_max_nokaon       = _getMaxL( L_pion_, L_proton_, L_electron_, L_muon_,     L_background_ );
            const double L_max_noproton     = _getMaxL( L_pion_, L_kaon_,   L_electron_, L_muon_,     L_background_ );
            const double L_max_noelectron   = _getMaxL( L_pion_, L_kaon_,   L_proton_,   L_muon_,     L_background_ );
            const double L_max_nomuon       = _getMaxL( L_pion_, L_kaon_,   L_proton_,   L_electron_, L_background_ );
            const double L_max_nobackground = _getMaxL( L_pion_, L_kaon_,   L_proton_,   L_electron_, L_muon_ );
            P_pion_       = ( L_pion_        > 0.0 )? ( ( L_max_nopion       > 0.0 )? L_pion_       / L_max_nopion       : 10.0 ) : -1.0;
            P_kaon_       = ( L_kaon_        > 0.0 )? ( ( L_max_nokaon       > 0.0 )? L_kaon_       / L_max_nokaon       : 10.0 ) : -1.0;
            P_proton_     = ( L_proton_      > 0.0 )? ( ( L_max_noproton     > 0.0 )? L_proton_     / L_max_noproton     : 10.0 ) : -1.0;
            P_electron_   = ( L_electron_    > 0.0 )? ( ( L_max_noelectron   > 0.0 )? L_electron_   / L_max_noelectron   : 10.0 ) : -1.0;
            P_muon_       = ( L_muon_        > 0.0 )? ( ( L_max_nomuon       > 0.0 )? L_muon_       / L_max_nomuon       : 10.0 ) : -1.0;
            P_background_ = ( L_background_  > 0.0 )? ( ( L_max_nobackground > 0.0 )? L_background_ / L_max_nobackground : 10.0 ) : -1.0;
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

    private:
        /**
		 * @return Maximum of all likelihoods which are >= 0.0
		 */
		double _getMaxL( const double L1, const double L2, const double L3, const double L4, const double L5){
			return fmax( 0.0, fmax( L1, fmax( L2, fmax( L3, fmax( L4, L5 ) ) ) ) );
		}

	};
	class CalcRICHPID: public CalcRICHProbabilities{
	public:

		/***
		 * If pid should be determined,
		 * 	pion: 0
		 * 	kaon: 1
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
							 double const* mom_mag,
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
							 mom_(mom),
							 mom_mag_(mom_mag),
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
			const double mom = (mom_)? mom_->Mag() : *mom_mag_;
            if(      P_pion_       > P_ratio_cut_ && mom > Mom_pion_min_     && mom < Mom_pion_max_)     pid_ = 0;
            else if( P_kaon_       > P_ratio_cut_ && mom > Mom_kaon_min_     && mom < Mom_kaon_max_)     pid_ = 1;
            else if( P_proton_     > P_ratio_cut_ && mom > Mom_proton_min_   && mom < Mom_proton_max_)   pid_ = 2;
            else if( P_electron_   > P_ratio_cut_ && mom > Mom_electron_min_ && mom < Mom_electron_max_) pid_ = 3;
            else if( P_muon_       > P_ratio_cut_ && mom > Mom_muon_min_     && mom < Mom_muon_max_)     pid_ = 4;
            else if( P_background_ > P_ratio_cut_ )                                                      pid_ = 5;
            else                                                                                         pid_ = -1;
			return true;
		}



	private:
		TVector3 const* mom_;
		double const* mom_mag_;
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
	 * is_kp_pk = -1 -> no decision possible
	 * is_kp_pk = 1  -> no first particle is kaon, second particle is pion
	 * is_kp_pk = 2  -> no first particle is pion, second particle is kaon
	 */
	class DetermineKaonPionLV: public Function {
	public:

	DetermineKaonPionLV(TVector3* mom_1, int const* pid_1, TVector3* mom_2, int const* pid_2, double const* mass_charged_kaon, double const* mass_charged_pion,
	                    TLorentzVector* kaon_lv, TLorentzVector* pion_lv, int* is_kp_pk, int* pid_kaon, int* pid_pion, int* pid_kaon_mct, int* pid_pion_mct ,const int* method, const int* pid_1_mct, const int* pid_2_mct) :
			mom_1_(*mom_1),
			pid_1_(*pid_1),
			mom_2_(*mom_2),
			pid_2_(*pid_2),
			mass_charged_kaon_(*mass_charged_kaon),
			mass_charged_pion_(*mass_charged_pion),
			kaon_lv_(*kaon_lv),
			pion_lv_(*pion_lv),
			is_kp_pk_(*is_kp_pk),
			pid_kaon_(*pid_kaon),
			pid_pion_(*pid_pion),
			pid_kaon_mct_(pid_kaon_mct),
			pid_pion_mct_(pid_pion_mct),
			method_(*method),
			pid_1_mct_(pid_1_mct),
			pid_2_mct_(pid_2_mct){
	}



		bool operator()() {
		is_kp_pk_ = -1;
		kaon_lv_ = TLorentzVector(-4444.0, -4444.0, -4444.0, -4444.0);
		pion_lv_ = TLorentzVector(-4444.0, -4444.0, -4444.0, -4444.0);
		switch (method_) {
			case 0: // identify K- or pi-
			{
				if ((pid_1_ == 1 && pid_2_ != 1) || (pid_2_ == 0 && pid_1_ != 0)) { // 1 = kaon, 2 = pion
					kaon_lv_ = TLorentzVector(mom_1_, sqrt(mass_charged_kaon_ * mass_charged_kaon_ + mom_1_.Mag2()));
					pion_lv_ = TLorentzVector(mom_2_, sqrt(mass_charged_pion_ * mass_charged_pion_ + mom_2_.Mag2()));
					is_kp_pk_ = 1;
					pid_kaon_ = pid_1_;
					pid_pion_ = pid_2_;
					if(pid_kaon_mct_ != nullptr){
						*pid_kaon_mct_ = *pid_1_mct_;
						*pid_pion_mct_ = *pid_2_mct_;
					}
				} else if ((pid_1_ == 0 && pid_2_ != 0) || (pid_2_ == 1 && pid_1_ != 1)) { // 1 = pion , 2 = kaon
					kaon_lv_ = TLorentzVector(mom_2_,
					        sqrt(mass_charged_kaon_ * mass_charged_kaon_ + mom_2_.Mag2()));
					pion_lv_ = TLorentzVector(mom_1_,
					        sqrt(mass_charged_pion_ * mass_charged_pion_ + mom_1_.Mag2()));
					is_kp_pk_ = 2;
					pid_kaon_ = pid_2_;
					pid_pion_ = pid_1_;
					if(pid_kaon_mct_ != nullptr){
						*pid_kaon_mct_ = *pid_2_mct_;
						*pid_pion_mct_ = *pid_1_mct_;
					}
				}
				break;
			}
			case 1: // identify K-
			{
				if ((pid_1_ == 1 && pid_2_ != 1)) { // 1 = kaon, 2 = pion
					kaon_lv_ = TLorentzVector(mom_1_, sqrt(mass_charged_kaon_ * mass_charged_kaon_ + mom_1_.Mag2()));
					pion_lv_ = TLorentzVector(mom_2_, sqrt(mass_charged_pion_ * mass_charged_pion_ + mom_2_.Mag2()));
					is_kp_pk_ = 1;
					pid_kaon_ = pid_1_;
					pid_pion_ = pid_2_;
					if(pid_kaon_mct_ != nullptr){
						*pid_kaon_mct_ = *pid_1_mct_;
						*pid_pion_mct_ = *pid_2_mct_;
					}
				} else if (pid_2_ == 1 && pid_1_ != 1) { // 1 = pion , 2 = kaon
					kaon_lv_ = TLorentzVector(mom_2_, sqrt(mass_charged_kaon_ * mass_charged_kaon_ + mom_2_.Mag2()));
					pion_lv_ = TLorentzVector(mom_1_, sqrt(mass_charged_pion_ * mass_charged_pion_ + mom_1_.Mag2()));
					is_kp_pk_ = 2;
					pid_kaon_ = pid_2_;
					pid_pion_ = pid_1_;
					if(pid_kaon_mct_ != nullptr){
						*pid_kaon_mct_ = *pid_2_mct_;
						*pid_pion_mct_ = *pid_1_mct_;
					}
				}
			break;
			}
			case 2: // identify pi-
			{
				if (pid_2_ == 0 && pid_1_ != 0) { // 1 = kaon, 2 = pion
					kaon_lv_ = TLorentzVector(mom_1_, sqrt(mass_charged_kaon_ * mass_charged_kaon_ + mom_1_.Mag2()));
					pion_lv_ = TLorentzVector(mom_2_, sqrt(mass_charged_pion_ * mass_charged_pion_ + mom_2_.Mag2()));
					is_kp_pk_ = 1;
					pid_kaon_ = pid_1_;
					pid_pion_ = pid_2_;
					if(pid_kaon_mct_ != nullptr){
						*pid_kaon_mct_ = *pid_1_mct_;
						*pid_pion_mct_ = *pid_2_mct_;
					}
				} else if (pid_1_ == 0 && pid_2_ != 0) { // 1 = pion , 2 = kaon
					kaon_lv_ = TLorentzVector(mom_2_,
					        sqrt(mass_charged_kaon_ * mass_charged_kaon_ + mom_2_.Mag2()));
					pion_lv_ = TLorentzVector(mom_1_,
					        sqrt(mass_charged_pion_ * mass_charged_pion_ + mom_1_.Mag2()));
					is_kp_pk_ = 2;
					pid_kaon_ = pid_2_;
					pid_pion_ = pid_1_;
					if(pid_kaon_mct_ != nullptr){
						*pid_kaon_mct_ = *pid_2_mct_;
						*pid_pion_mct_ = *pid_1_mct_;
					}
				}
				break;
			}
			default:
				throw "Method not implemented for DetermineKaonPionLV.";
		}
		return true;
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
		int* const pid_kaon_mct_;
		int* const pid_pion_mct_;
		const int& method_;
		const int* const pid_1_mct_;
		const int* const pid_2_mct_;

	};


	/***
	 * is_kp_pk = -1 -> no decision possible
	 * is_kp_pk = 1  -> no first particle is kaon, second particle is pion
	 * is_kp_pk = 2  -> no first particle is pion, second particle is kaon
	 */
	class DetermineKaonPionLVLikelihood: public Function {
	public:

	DetermineKaonPionLVLikelihood(TVector3* mom_1, double const* l_pion_1, double const* l_kaon_1, int const* pid_1,
	                              TVector3* mom_2, double const* l_pion_2, double const* l_kaon_2, int const* pid_2,
	                              double const* mass_charged_kaon, double const* mass_charged_pion,
	                              double const* threshold_kpi, double const* threshold_pik,
	                              TLorentzVector* kaon_lv, TLorentzVector* pion_lv, int* is_kp_pk, int* pid_kaon, int* pid_pion) :
			mom_1_(*mom_1),
			l_pion_1_(*l_pion_1),
			l_kaon_1_(*l_kaon_1),
			pid_1_(*pid_1),
			mom_2_(*mom_2),
			l_pion_2_(*l_pion_2),
			l_kaon_2_(*l_kaon_2),
			pid_2_(*pid_2),
			mass_charged_kaon_(*mass_charged_kaon),
			mass_charged_pion_(*mass_charged_pion),
			threshold_kpi_(*threshold_kpi),
			threshold_pik_(*threshold_pik),
			kaon_lv_(*kaon_lv),
			pion_lv_(*pion_lv),
			is_kp_pk_(*is_kp_pk),
			pid_kaon_(*pid_kaon),
			pid_pion_(*pid_pion){
	}



	bool operator()() {
		is_kp_pk_ = -1;
		kaon_lv_ = TLorentzVector(-4444.0, -4444.0, -4444.0, -4444.0);
		pion_lv_ = TLorentzVector(-4444.0, -4444.0, -4444.0, -4444.0);

		const double l_kpi = l_kaon_1_ * l_pion_2_;
		const double l_pik = l_pion_1_ * l_kaon_2_;
		const double log_l_ratio_kpi_piK = log(l_kpi / l_pik);

		if (log_l_ratio_kpi_piK > threshold_kpi_) {
			is_kp_pk_ = 1;
			kaon_lv_ = TLorentzVector(mom_1_, sqrt(mass_charged_kaon_ * mass_charged_kaon_ + mom_1_.Mag2()));
			pion_lv_ = TLorentzVector(mom_2_, sqrt(mass_charged_pion_ * mass_charged_pion_ + mom_2_.Mag2()));
			pid_kaon_ = pid_1_;
			pid_pion_ = pid_2_;
		} else if (log_l_ratio_kpi_piK < threshold_pik_) {
			is_kp_pk_ = 2;
			kaon_lv_ = TLorentzVector(mom_2_, sqrt(mass_charged_kaon_ * mass_charged_kaon_ + mom_2_.Mag2()));
			pion_lv_ = TLorentzVector(mom_1_, sqrt(mass_charged_pion_ * mass_charged_pion_ + mom_1_.Mag2()));
			pid_kaon_ = pid_2_;
			pid_pion_ = pid_1_;
		}
		return true;
	}



	private:
		TVector3 const& mom_1_;
		double const& l_pion_1_;
		double const& l_kaon_1_;
		int const& pid_1_;
		TVector3 const& mom_2_;
		double const& l_pion_2_;
		double const& l_kaon_2_;
		int const& pid_2_;
		double const& mass_charged_kaon_;
		double const& mass_charged_pion_;
		double const& threshold_kpi_;
		double const& threshold_pik_;
		TLorentzVector& kaon_lv_;
		TLorentzVector& pion_lv_;
		int& is_kp_pk_;
		int& pid_kaon_;
		int& pid_pion_;

	};

	class CalcCEDARPID: public Function {
	public:

		/***
		 * If pid should be determined,
		 *  noPID: -1
		 * 	kaon: 0
		 * 	pion: 1
		 * 	proton: 2
		 *
		 * This method takes into account the likelihood of each CEDAR individually and finally combines both CEDARs.
		 *
		 * The Thresholds are given for the difference of the log10 of the likelihood to be a kaon minus the log10 of the likelihood to be a pion
		 *   - the difference must be larger than the kaon threshold to be identified as a kaon
		 *   - the difference must be smaller than the pion threshold to be identified as a pion
		 * A PID for that CEDAR is assigned if exactly one of the two conditions is fullfilled
		 * A overall PID is assigned if both CEDARs give the same PID or one of the CEDARs has no PID
		 *
		 *
		 */

		enum PIDs {pidNo = -1, pidKaon = 0, pidPion = 1, pidProton = 2};
		CalcCEDARPID(double const* L_pion_CEDAR1,
		            double const* L_kaon_CEDAR1,
		            double const* L_proton_CEDAR1,
		            int const* n_hits_CEDAR1,
		            std::vector<double> const& thresholds_kaon_DeltaLogLike_CEDAR1,
		            std::vector<double> const& thresholds_pion_DeltaLogLike_CEDAR1,
		            double const* L_pion_CEDAR2,
		            double const* L_kaon_CEDAR2,
		            double const* L_proton_CEDAR2,
		            int const* n_hits_CEDAR2,
		            std::vector<double> const& thresholds_kaon_DeltaLogLike_CEDAR2,
		            std::vector<double> const& thresholds_pion_DeltaLogLike_CEDAR2,
		            int* CEDARPid_CEDAR1,
		            int* CEDARPid_CEDAR2,
		            int* CEDARPid,
		            double* LLdiff_CEDAR1,
		            double* LLdiff_CEDAR2
		            ):
				L_pion_CEDAR1_(*L_pion_CEDAR1),
				L_kaon_CEDAR1_(*L_kaon_CEDAR1),
				L_proton_CEDAR1_(*L_proton_CEDAR1),
				n_hits_CEDAR1_(*n_hits_CEDAR1),
				thresholds_kaon_DeltaLogLike_CEDAR1_(thresholds_kaon_DeltaLogLike_CEDAR1),
				thresholds_pion_DeltaLogLike_CEDAR1_(thresholds_pion_DeltaLogLike_CEDAR1),
				CEDARPid_CEDAR1_(*CEDARPid_CEDAR1),
				L_pion_CEDAR2_(*L_pion_CEDAR2),
				L_kaon_CEDAR2_(*L_kaon_CEDAR2),
				L_proton_CEDAR2_(*L_proton_CEDAR2),
				n_hits_CEDAR2_(*n_hits_CEDAR2),
				thresholds_kaon_DeltaLogLike_CEDAR2_(thresholds_kaon_DeltaLogLike_CEDAR2),
				thresholds_pion_DeltaLogLike_CEDAR2_(thresholds_pion_DeltaLogLike_CEDAR2),
				CEDARPid_CEDAR2_(*CEDARPid_CEDAR2),
				CEDARPid_(*CEDARPid),
				LLdiff_CEDAR1_( *((LLdiff_CEDAR1)? LLdiff_CEDAR1 : new double())),
				LLdiff_CEDAR2_( *((LLdiff_CEDAR2)? LLdiff_CEDAR2 : new double()))
		{
		}

		bool operator()() {

			CEDARPid_CEDAR1_ = getPIDForCEDAR(L_kaon_CEDAR1_, L_pion_CEDAR1_, n_hits_CEDAR1_, thresholds_kaon_DeltaLogLike_CEDAR1_, thresholds_pion_DeltaLogLike_CEDAR1_, LLdiff_CEDAR1_);
			CEDARPid_CEDAR2_ = getPIDForCEDAR(L_kaon_CEDAR2_, L_pion_CEDAR2_, n_hits_CEDAR2_, thresholds_kaon_DeltaLogLike_CEDAR2_, thresholds_pion_DeltaLogLike_CEDAR2_, LLdiff_CEDAR2_);


			if (CEDARPid_CEDAR1_ == CEDARPid_CEDAR2_ )      CEDARPid_ = CEDARPid_CEDAR1_;
			else if (CEDARPid_CEDAR2_ == pidNo)             CEDARPid_ = CEDARPid_CEDAR1_;
			else if (CEDARPid_CEDAR1_ == pidNo)             CEDARPid_ = CEDARPid_CEDAR2_;
			else                                            CEDARPid_ = pidNo;
			return true;
		}

	protected:
		double const& L_pion_CEDAR1_;
		double const& L_kaon_CEDAR1_;
		double const& L_proton_CEDAR1_;
		const int& n_hits_CEDAR1_;
		const std::vector<double> thresholds_kaon_DeltaLogLike_CEDAR1_;
		const std::vector<double> thresholds_pion_DeltaLogLike_CEDAR1_;
		int& CEDARPid_CEDAR1_;
		double const& L_pion_CEDAR2_;
		double const& L_kaon_CEDAR2_;
		double const& L_proton_CEDAR2_;
		const int& n_hits_CEDAR2_;
		const std::vector<double> thresholds_kaon_DeltaLogLike_CEDAR2_;
		const std::vector<double> thresholds_pion_DeltaLogLike_CEDAR2_;
		int& CEDARPid_CEDAR2_;
		int& CEDARPid_;
		double& LLdiff_CEDAR1_;
		double& LLdiff_CEDAR2_;


	private:

		int getPIDForCEDAR(const double L_kaon, const double L_pion, const int n_hits,
		                   const std::vector<double>& thresholds_kaon, const std::vector<double>& thresholds_pion,
		                   double& LogLikeDiff) {
			int pid = pidNo;
			LogLikeDiff = std::numeric_limits<double>::quiet_NaN();
			if(L_kaon <= 0.0 and L_pion > 0.0){
				pid = pidPion;
				LogLikeDiff = -15.0;
			} else if (L_pion <= 0.0 and L_kaon > 0.0){
				pid = pidKaon;
				LogLikeDiff = 15.0;
			} else if (L_kaon > 0.0 and L_pion > 0.0) {
				LogLikeDiff = log10(L_kaon) - log10(L_pion);
				const bool is_kaon = LogLikeDiff > thresholds_kaon[n_hits];
				const bool is_pion = LogLikeDiff < thresholds_pion[n_hits];
				if        (is_kaon and not is_pion) {
					pid = pidKaon;
				} else if (is_pion and not is_kaon) {
					pid = pidPion;
				} // else stay with no pid
			}

			return pid;
		}
	};


	class CalcCEDARPIDOneL: public Function {
	public:

		/***
		 * If pid should be determined,
		 *  noPID: -1
		 * 	kaon: 0
		 * 	pion: 1
		 * 	proton: 2
		 *
		 * The Thresholds are given for the difference of the log10 of the likelihood to be a kaon minus the log10 of the likelihood to be a pion
		 *   - the difference must be larger than the kaon threshold to be identified as a kaon
		 *   - the difference must be smaller than the pion threshold to be identified as a pion
		 *
		 *
		 */

		enum PIDs {pidNo = -1, pidKaon = 0, pidPion = 1, pidProton = 2};
		CalcCEDARPIDOneL(double const* L_pion_CEDAR1,
		            double const* L_kaon_CEDAR1,
		            double const* L_proton_CEDAR1,
		            double const* L_pion_CEDAR2,
		            double const* L_kaon_CEDAR2,
		            double const* L_proton_CEDAR2,
		            double const* threshold_kaon_DeltaLogLike,
		            double const* threshold_pion_DeltaLogLike,
		            int* CEDARPid,
		            double* LLdiff_CEDAR1,
		            double* LLdiff_CEDAR2,
		            double* LLdiff_CEDARs
		            ):
				L_pion_CEDAR1_(*L_pion_CEDAR1),
				L_kaon_CEDAR1_(*L_kaon_CEDAR1),
				L_proton_CEDAR1_(*L_proton_CEDAR1),
				L_pion_CEDAR2_(*L_pion_CEDAR2),
				L_kaon_CEDAR2_(*L_kaon_CEDAR2),
				L_proton_CEDAR2_(*L_proton_CEDAR2),
				threshold_kaon_DeltaLogLike_(*threshold_kaon_DeltaLogLike),
				threshold_pion_DeltaLogLike_(*threshold_pion_DeltaLogLike),
				CEDARPid_(*CEDARPid),
				LLdiff_CEDAR1_( *((LLdiff_CEDAR1)? LLdiff_CEDAR1 : new double())),
				LLdiff_CEDAR2_( *((LLdiff_CEDAR2)? LLdiff_CEDAR2 : new double())),
				LLdiff_CEDARs_( *((LLdiff_CEDARs)? LLdiff_CEDARs : new double()))
		{
		}

		bool operator()() {

			LLdiff_CEDAR1_ = getLogLikeDiff(L_kaon_CEDAR1_, L_pion_CEDAR1_);
			LLdiff_CEDAR2_ = getLogLikeDiff(L_kaon_CEDAR2_, L_pion_CEDAR2_);
			LLdiff_CEDARs_ = getLogLikeDiff(L_kaon_CEDAR2_*L_kaon_CEDAR1_, L_pion_CEDAR2_*L_pion_CEDAR1_);


			if      (LLdiff_CEDARs_ > threshold_kaon_DeltaLogLike_ )  CEDARPid_ = pidKaon;
			else if (LLdiff_CEDARs_ < threshold_pion_DeltaLogLike_ )  CEDARPid_ = pidPion;
			else                                                      CEDARPid_ = pidNo;
			return true;
		}

	protected:
		double const& L_pion_CEDAR1_;
		double const& L_kaon_CEDAR1_;
		double const& L_proton_CEDAR1_;
		double const& L_pion_CEDAR2_;
		double const& L_kaon_CEDAR2_;
		double const& L_proton_CEDAR2_;
		double const& threshold_kaon_DeltaLogLike_;
		double const& threshold_pion_DeltaLogLike_;
		int& CEDARPid_;
		double& LLdiff_CEDAR1_;
		double& LLdiff_CEDAR2_;
		double& LLdiff_CEDARs_;


	private:

		double getLogLikeDiff(const double L_kaon, const double L_pion) {
			double LogLikeDiff = std::numeric_limits<double>::quiet_NaN();
			if(L_kaon <= 0.0 and L_pion > 0.0){
				LogLikeDiff = -15.0;
			} else if (L_pion <= 0.0 and L_kaon > 0.0){
				LogLikeDiff = 15.0;
			} else if (L_kaon > 0.0 and L_pion > 0.0) {
				LogLikeDiff = log10(L_kaon) - log10(L_pion);
			}

			return LogLikeDiff;
		}
	};


	class CalcAngles3P: public Function {
	public:

		/**
		 * Calculates the Gottfried-Jackson and Helicity frame angles from the four-momenta (without symmetrization).
		 * Does only for for 3-particle final states
		 *
		 * @param lv11Addr Four momentum of the batchelor
		 * @param lv21Addr Four momentum of the first isobar daughter
		 * @param lv22Addr Four momentum of the second isobar daughter
		 * @param lvBeam Four momentum of the beam particle
		 * @param targetMassAddr Target mass
		 * @param GJ_costhetaAddr Gottfried-Jackson frame costheta angle
		 * @param GJ_phiAddr Gottfried-Jackson frame phi angle
		 * @param HF_costhetaAddr Helicity frame costheta angle of the isobar decay
		 * @param HF_phiAddr Helicity frame  phi angle of the isobar decay
		 */
		CalcAngles3P(const TLorentzVector* lv11Addr, const TLorentzVector* lv21Addr, const TLorentzVector* lv22Addr,
		             const TLorentzVector* lvBeamAddr, const double* targetMassAddr,
		             double* GJ_costhetaAddr, double* GJ_phiAddr, double* HF_costhetaAddr, double* HF_phiAddr):
			lv11_(*lv11Addr),
			lv21_(*lv21Addr),
			lv22_(*lv22Addr),
			lvBeam_(*lvBeamAddr),
			targetMass_(*targetMassAddr),
			GJ_costheta_(*GJ_costhetaAddr),
			GJ_phi_(*GJ_phiAddr),
			HF_costheta_(*HF_costhetaAddr),
			HF_phi_(*HF_phiAddr)
			{}



		bool operator()() {

			const TLorentzVector lvIsobar = lv21_ + lv22_; // LV of isobar
			const TLorentzVector lvX = lv11_ + lv21_ + lv22_; // LV of X
			const TLorentzVector lvTarget(0,0,0,targetMass_);

			const TVector3 boostLab2X(-lvX.BoostVector());

			// boost in X rest frame
			TLorentzVector lvX_X(lvX);                          lvX_X.Boost(boostLab2X);
			TLorentzVector lvBeam_X(lvBeam_);                   lvBeam_X.Boost(boostLab2X);
			TLorentzVector lvTarget_X(lvTarget);                lvTarget_X.Boost(boostLab2X);
			TLorentzVector lvRecoil_X(lvBeam_ - lvX+ lvTarget);	lvRecoil_X.Boost(boostLab2X);
			TLorentzVector lvIsobar_X(lvIsobar);                lvIsobar_X.Boost(boostLab2X);
			TLorentzVector lv11_X(lv11_);                       lv11_X.Boost(boostLab2X);
			TLorentzVector lv21_X(lv21_);                       lv21_X.Boost(boostLab2X);
			TLorentzVector lv22_X(lv22_);                       lv22_X.Boost(boostLab2X);



			// get z-, y-,x- axis of GJ frame in X rest frame
			const TVector3 gjAxisZ(lvBeam_X.Vect().Unit());
			const TVector3 gjAxisY(lvTarget_X.Vect().Cross(lvRecoil_X.Vect()).Unit());
			const TVector3 gjAxisX(gjAxisY.Cross(gjAxisZ).Unit());

			// get rotation from GJ to lab frame and vice versa
			TRotation GJ2X;
			GJ2X = GJ2X.RotateAxes(gjAxisX, gjAxisY, gjAxisZ);
			TRotation X2GJ(GJ2X);
			X2GJ.Invert();
			assert(!X2GJ.IsIdentity());

			// rotate in GJ frame
			TLorentzVector lvX_GJ(lvX_X);            lvX_GJ *= X2GJ;
			TLorentzVector lvBeam_GJ(lvBeam_X);      lvBeam_GJ *= X2GJ;
			TLorentzVector lvTarget_GJ(lvTarget_X);  lvTarget_GJ *= X2GJ;
			TLorentzVector lvRecoil_GJ(lvRecoil_X);  lvRecoil_GJ *= X2GJ;
			TLorentzVector lvIsobar_GJ(lvIsobar_X);  lvIsobar_GJ *= X2GJ;
			TLorentzVector lv11_GJ(lv11_X);          lv11_GJ *= X2GJ;
			TLorentzVector lv21_GJ(lv21_X);          lv21_GJ *= X2GJ;
			TLorentzVector lv22_GJ(lv22_X);          lv22_GJ *= X2GJ;

			// calculate theta, phi in GJ frame
			GJ_costheta_ = lvIsobar_GJ.CosTheta();
			GJ_phi_ = lvIsobar_GJ.Phi();

			// boost in isobar rest frame
			const TVector3 boostX2I(-lvIsobar_X.BoostVector());
			TLorentzVector lvIsobar_I(lvIsobar_X);     lvIsobar_I.Boost(boostX2I);
			TLorentzVector lv21_I(lv21_X);             lv21_I.Boost(boostX2I);
			TLorentzVector lv22_I(lv22_X);             lv22_I.Boost(boostX2I);

			// get z-, y-, x-axis in helicity frame
			const TVector3 axisZ_Hel(lvIsobar_X.Vect().Unit());
			const TVector3 axisY_Hel(lvBeam_X.Vect().Cross(axisZ_Hel).Unit());
			const TVector3 axisX_Hel(axisY_Hel.Cross(axisZ_Hel));

			// get rotation in helicity frame
			TRotation Hel2I;
			Hel2I = Hel2I.RotateAxes(axisX_Hel, axisY_Hel, axisZ_Hel);
			TRotation I2Hel(Hel2I);
			I2Hel.Invert();
			assert(!I2Hel.IsIdentity());

			TLorentzVector lvIsobar_Hel(lvIsobar_I);      lvIsobar_Hel *= I2Hel;
			TLorentzVector lv21_Hel(lv21_I);              lv21_Hel *= I2Hel;
			TLorentzVector lv22_Hel(lv22_I);              lv22_Hel *= I2Hel;

			HF_costheta_ = lv21_Hel.CosTheta();
			HF_phi_ = lv21_Hel.Phi();

			return true;
		}

	protected:
		const TLorentzVector& lv11_;
		const TLorentzVector& lv21_;
		const TLorentzVector& lv22_;
		const TLorentzVector& lvBeam_;
		const double& targetMass_;
		double& GJ_costheta_;
		double& GJ_phi_;
		double& HF_costheta_;
		double& HF_phi_;


	private:

		double getLogLikeDiff(const double L_kaon, const double L_pion) {
			double LogLikeDiff = std::numeric_limits<double>::quiet_NaN();
			if(L_kaon <= 0.0 and L_pion > 0.0){
				LogLikeDiff = -15.0;
			} else if (L_pion <= 0.0 and L_kaon > 0.0){
				LogLikeDiff = 15.0;
			} else if (L_kaon > 0.0 and L_pion > 0.0) {
				LogLikeDiff = log10(L_kaon) - log10(L_pion);
			}

			return LogLikeDiff;
		}
	};
}
}
}
}



#endif /* USER_SWALLNER_FUNCTIONS_HPP_ */
