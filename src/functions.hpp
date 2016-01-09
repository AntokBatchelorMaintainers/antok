#ifndef ANTOK_FUNCTIONS_HPP
#define ANTOK_FUNCTIONS_HPP

#include<iostream>
#include<vector>

#include<TLorentzVector.h>

#include<basic_calcs.h>

namespace antok {

	class Function
	{

	  public:

		virtual bool operator() () = 0;
		virtual ~Function() { }

	};

	namespace functions {

		template<typename T>
		class Sum : public Function
		{

		  public:

			Sum(std::vector<T*> inputAddrsSummands, std::vector<T*> inputAddrsSubtrahends, T* outAddr) {
				_inputAddrsSummands = inputAddrsSummands;
				_inputAddrsSubtrahends = inputAddrsSubtrahends;
				_outAddr = outAddr;
				if( (_inputAddrsSummands.size() < 1) and (_inputAddrsSubtrahends.size() < 1) ) {
					std::cerr<<"Got empty address vector as input for a sum."<<std::endl;
					throw 1;
				}
			};

			virtual ~Sum() { }

			bool operator() () {
				(*_outAddr) = T();
				for(unsigned int i = 0; i < _inputAddrsSummands.size(); ++i) {
					*_outAddr = *_outAddr + *(_inputAddrsSummands[i]);
				}
				for(unsigned int i = 0; i < _inputAddrsSubtrahends.size(); ++i) {
					*_outAddr = *_outAddr - *(_inputAddrsSubtrahends[i]);
				}
				return true;
			};


		  private:

			std::vector<T*> _inputAddrsSummands;
			std::vector<T*> _inputAddrsSubtrahends;
			T* _outAddr;

		};



		class Mass: public Function
		{

		  public:

			Mass(TLorentzVector* inputAddr, double* outAddr)
				: _inputAddr(inputAddr),
				  _outAddr(outAddr) { }

			virtual ~Mass() { }

			bool operator() () {
				(*_outAddr) = _inputAddr->M();
				return true;
			}

		  private:

			TLorentzVector* _inputAddr;
			double* _outAddr;

		};

		class GetLorentzVec: public Function
		{

		  public:

			GetLorentzVec(double* xAddr, double* yAddr, double* zAddr, double* mAddr, TLorentzVector* outAddr, int pType)
				: _xAddr(xAddr),
				  _yAddr(yAddr),
				  _zAddr(zAddr),
				  _mAddr(mAddr),
				  _outAddr(outAddr),
				  _pType(pType) { }

			GetLorentzVec(TVector3* vec3Addr, double* mAddr, TLorentzVector* outAddr, int pType)
				: _vec3Addr(vec3Addr),
				  _mAddr(mAddr),
				  _outAddr(outAddr),
				  _pType(pType) { }


			virtual ~GetLorentzVec() { }

			bool operator() () {
				switch( _pType )
				{
					case 0:
						_outAddr->SetXYZM(*_xAddr, *_yAddr, *_zAddr, *_mAddr);
						break;
					case 1:
						_outAddr->SetPxPyPzE(*_xAddr, *_yAddr, *_zAddr, *_mAddr);
						break;
					case 2:
						_outAddr->SetXYZM(_vec3Addr->X(),_vec3Addr->Y(),_vec3Addr->Z(),  *_mAddr);
						break;
					case 3:
						_outAddr->SetPxPyPzE(_vec3Addr->X(),_vec3Addr->Y(),_vec3Addr->Z(), *_mAddr);
						break;
				}
				return true;
			}

		  private:

			double* _xAddr;
			double* _yAddr;
			double* _zAddr;
			TVector3* _vec3Addr;
			double* _mAddr;
			TLorentzVector* _outAddr;
			int _pType;

		};

		class GetBeamLorentzVec : public Function
		{

		  public:

			GetBeamLorentzVec(double* gradxAddr, double* gradyAddr, TLorentzVector* xLorentzVec, TLorentzVector* outAddr)
				: _gradxAddr(gradxAddr),
				  _gradyAddr(gradyAddr),
				  _xLorentzVec(xLorentzVec),
				  _outAddr(outAddr) { }

			virtual ~GetBeamLorentzVec() { }

			bool operator() () {
				TVector3 p3Beam((*_gradxAddr), (*_gradyAddr), 1.);
				(*_outAddr) = antok::getBeamEnergy(p3Beam, (*_xLorentzVec));
				return true;
			}

		  private:

			double* _gradxAddr;
			double* _gradyAddr;
			TLorentzVector* _xLorentzVec;
			TLorentzVector* _outAddr;

		};

		class GetTs : public Function
		{

		  public:

			GetTs(TLorentzVector* xLorentzVec, TLorentzVector* beamLorentzVec, double* tAddr, double* tMinAddr, double* tPrimeAddr)
				: _xLorentzVec(xLorentzVec),
				  _beamLorentzVec(beamLorentzVec),
				  _tAddr(tAddr),
				  _tMinAddr(tMinAddr),
				  _tPrimeAddr(tPrimeAddr) { }

			virtual ~GetTs() { }

			bool operator() () {
			    (*_tAddr) = std::fabs(((*_beamLorentzVec) - (*_xLorentzVec)).Mag2());
			    (*_tMinAddr) = std::fabs((std::pow((*_xLorentzVec).M2() - (*_beamLorentzVec).M2(), 2)) / (4. * (_beamLorentzVec->Vect()).Mag2()));
			    (*_tPrimeAddr) = (*_tAddr) - (*_tMinAddr);
				return true;
			}

		  private:

			TLorentzVector* _xLorentzVec;
			TLorentzVector* _beamLorentzVec;
			double* _tAddr;
			double* _tMinAddr;
			double* _tPrimeAddr;

		};

		class Sum2 : public Function
		{

		  public:

			Sum2(std::vector<double*> inAddrs, double* outAddr)
				: _inAddrs(inAddrs),
				  _outAddr(outAddr) { }

			virtual ~Sum2() { }

			bool operator() () {
				(*_outAddr) = 0.;
				for(unsigned int i = 0; i < _inAddrs.size(); ++i) {
					(*_outAddr) += ((*_inAddrs[i]) * (*_inAddrs[i]));
				}
				(*_outAddr) = std::sqrt(*_outAddr);
				return true;
			}

		  private:

			std::vector<double*> _inAddrs;
			double* _outAddr;

		};

		class Diff : public Function
		{

		  public:

			Diff(double* inAddr1, double* inAddr2, double* outAddr)
				: _inAddr1(inAddr1),
				  _inAddr2(inAddr2),
				  _outAddr(outAddr) { }

			virtual ~Diff() { }

			bool operator() () {
				(*_outAddr) = (*_inAddr1) - (*_inAddr2);
				return true;
			}

		  private:

			double* _inAddr1;
			double* _inAddr2;
			double* _outAddr;

		};

		template < typename T >
		class Quotient: public Function
		{

		  public:

			Quotient(T* inAddr1, T* inAddr2, T* outAddr)
				: _inAddr1(inAddr1),
				  _inAddr2(inAddr2),
				  _outAddr(outAddr) { }

			virtual ~Quotient() { }

			bool operator() () {
				(*_outAddr) = (*_inAddr1) / (*_inAddr2);
				return true;
			}

		  private:

			T* _inAddr1;
			T* _inAddr2;
			T* _outAddr;

		};

		template < typename T >
		class Mul: public Function
		{

		  public:

			Mul(T* inAddr1, T* inAddr2, T* outAddr)
				: _inAddr1(inAddr1),
				  _inAddr2(inAddr2),
				  _outAddr(outAddr) { }

			virtual ~Mul() { }

			bool operator() () {
				(*_outAddr) = (*_inAddr1) * (*_inAddr2);
				return true;
			}

		  private:

			T* _inAddr1;
			T* _inAddr2;
			T* _outAddr;

		};


		class Abs : public Function
		{

		  public:

			Abs(double* inAddr, double* outAddr)
				: _inAddr(inAddr),
				  _outAddr(outAddr) { }

			virtual ~Abs() { }

			bool operator() () {
				(*_outAddr) = std::fabs(*_inAddr);
				return true;
			}

		  private:

			double* _inAddr;
			double* _outAddr;

		};

		class Energy : public Function
		{

		  public:

			Energy(TLorentzVector* inAddr, double* outAddr)
				: _inAddr(inAddr),
				  _outAddr(outAddr) { }

			virtual ~Energy() { }

			bool operator() () {
				(*_outAddr) = _inAddr->E();
				return true;
			}

		  private:

			TLorentzVector* _inAddr;
			double* _outAddr;

		};

		class RadToDegree : public Function
		{

		  public:

			RadToDegree(double* inAddr, double* outAddr)
				: _inAddr(inAddr),
				  _outAddr(outAddr) { }

			virtual ~RadToDegree() { }

			bool operator() () {
				(*_outAddr) = ((*_inAddr) / TMath::Pi()) * 180.;
				return true;
			}

		  private:

			double* _inAddr;
			double* _outAddr;

		};

		class ConvertIntToDouble : public Function
		{

		  public:

			ConvertIntToDouble(int* inAddr, double* outAddr)
				: _inAddr(inAddr),
				  _outAddr(outAddr) { }

			virtual ~ConvertIntToDouble() { }

			bool operator() () {
				(*_outAddr) = (*_inAddr);
				return true;
			}

		  private:

			int* _inAddr;
			double* _outAddr;

		};

		class GetGradXGradY : public Function
		{

		  public:

			GetGradXGradY(TLorentzVector* lorentzVecAddr, double* xGradAddr, double* yGradAddr)
				: _lorentzVecAddr(lorentzVecAddr),
				  _xGradAddr(xGradAddr),
				  _yGradAddr(yGradAddr) { }

			virtual ~GetGradXGradY() { }

			bool operator() () {
				TVector3 vect = _lorentzVecAddr->Vect();
				double z = vect.Z();
				(*_xGradAddr) = vect.X() / z;
				(*_yGradAddr) = vect.Y() / z;
				return true;
			}

		  private:

			TLorentzVector* _lorentzVecAddr;
			double* _xGradAddr;
			double* _yGradAddr;

		};

		class GetLorentzVectorAttributes : public Function
		{

		  public:

			GetLorentzVectorAttributes(TLorentzVector* lorentzVecAddr, double* xAddr,
			                                                           double* yAddr,
			                                                           double* zAddr,
			                                                           double* phiAddr,
			                                                           double* thetaAddr)
				: _lorentzVecAddr(lorentzVecAddr),
				  _xAddr(xAddr),
				  _yAddr(yAddr),
				  _zAddr(zAddr),
				  _phiAddr(phiAddr),
				  _thetaAddr(thetaAddr) { }

			virtual ~GetLorentzVectorAttributes() { }

			bool operator() () {
				(*_xAddr) = _lorentzVecAddr->X();
				(*_yAddr) = _lorentzVecAddr->Y();
				(*_zAddr) = _lorentzVecAddr->Z();
				(*_phiAddr) = _lorentzVecAddr->Phi();
				(*_thetaAddr) = _lorentzVecAddr->Theta();
				return true;
			}

		  private:

			TLorentzVector* _lorentzVecAddr;
			double* _xAddr;
			double* _yAddr;
			double* _zAddr;
			double* _phiAddr;
			double* _thetaAddr;

		};

		class GetTVector3: public Function
		{

		  public:

			GetTVector3(double* xAddr, double* yAddr, double* zAddr, TVector3* outAddr)
				: _xAddr(xAddr),
				  _yAddr(yAddr),
				  _zAddr(zAddr),
				  _outAddr(outAddr) { }

			virtual ~GetTVector3() { }

			bool operator() () {
				_outAddr->SetXYZ(*_xAddr, *_yAddr, *_zAddr);
				return true;
			}

		  private:

			double* _xAddr;
			double* _yAddr;
			double* _zAddr;
			TVector3* _outAddr;

		};

	}

}

#endif

