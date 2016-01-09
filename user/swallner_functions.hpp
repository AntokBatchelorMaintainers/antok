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

}
}
}
}




#endif /* USER_SWALLNER_FUNCTIONS_HPP_ */
