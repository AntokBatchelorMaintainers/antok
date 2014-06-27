#ifndef ANTOK_USER_CDREIS_FUNCTIONS_HPP
#define ANTOK_USER_CDREIS_FUNCTIONS_HPP

#include<iostream>
#include<vector>


#include<TLorentzVector.h>
#include<TRotation.h>
#include<TLorentzRotation.h>

namespace antok {

	namespace user {

		namespace cdreis {

			namespace functions {

				class GetRecoilLorentzVec : public Function
				{
					public:
						GetRecoilLorentzVec(TLorentzVector* BeamLorentzVec, TLorentzVector* XLorentzVec, double* RecoilMass, TLorentzVector* RecoilLorentzVec)
						                   : _BeamLorentzVec(BeamLorentzVec), _XLorentzVec(XLorentzVec), _RecoilMass(RecoilMass),
						                     _RecoilLorentzVec(RecoilLorentzVec) { }

						virtual ~GetRecoilLorentzVec() { }

						bool operator() () {
							_RecoilLorentzVec->SetVectM(((_BeamLorentzVec->Vect()) - (_XLorentzVec->Vect())), *_RecoilMass);
							return true;
						}

					private:
						TLorentzVector *_BeamLorentzVec, *_XLorentzVec;
						double *_RecoilMass;
						TLorentzVector* _RecoilLorentzVec;
				};

			}

		}

	}

}

#endif
