#ifndef ANTOK_USER_HUBERS_FUNCTIONS_HPP
#define ANTOK_USER_HUBERS_FUNCTIONS_HPP

#include<iostream>
#include<vector>
#include<NNpoly.h>

#include<TLorentzVector.h>

namespace antok {

	namespace user {

		namespace hubers {

			namespace functions {

				class GetNeuronalBeamEnergy : public Function
				{
					public:
						GetNeuronalBeamEnergy(double* xAddr,
								double* yAddr,
								double* dxAddr,
								double* dyAddr,
								double* eAddr)
							: _xAddr(xAddr),
							_yAddr(yAddr),
							_dxAddr(dxAddr),
							_dyAddr(dyAddr),
							_eAddr(eAddr) { }

						virtual ~GetNeuronalBeamEnergy() { }

						bool operator() () {
 							double xarr[4]={*_xAddr,*_yAddr,*_dxAddr,*_dyAddr};
 							*_eAddr=					NNpoly::Ebeam(xarr,NNpoly::params);
							return true;
						}

					private:

						double* _xAddr;
						double* _yAddr;
						double* _dxAddr;
						double* _dyAddr;
						double* _eAddr;

				};






			}

		}

	}

}


#endif
