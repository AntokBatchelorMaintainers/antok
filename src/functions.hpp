#ifndef ANTOK_FUNCTIONS_HPP
#define ANTOK_FUNCTIONS_HPP

#include<iostream>
#include<vector>

#include<TLorentzVector.h>

#include<basic_calcs.h>
#include<constants.h>

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
                  _vec3Addr( nullptr ),
				  _mAddr(mAddr),
				  _outAddr(outAddr),
				  _pType(pType) { }

			GetLorentzVec(TVector3* vec3Addr, double* mAddr, TLorentzVector* outAddr, int pType)
				: _xAddr(NULL),
				  _yAddr(NULL),
				  _zAddr(NULL),
				  _vec3Addr(vec3Addr),
				  _mAddr(mAddr),
				  _outAddr(outAddr),
				  _pType(pType) {}

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

			GetBeamLorentzVec(double* gradxAddr, double* gradyAddr, TLorentzVector* xLorentzVec, TLorentzVector* outAddr, const double* mass_beam, const double* massTarget)
				: _gradx(*gradxAddr),
				  _grady(*gradyAddr),
				  _xLorentzVec(*xLorentzVec),
				  _massBeamAddr( mass_beam ),
				  _out(*outAddr),
				  _massTargetAddr(massTarget){
				if(_massTargetAddr == nullptr){
					std::cout << "INFO: No target mass is given to reconstruct the beam energy -> using proton mass" << std::endl;
					_massTargetAddr = &antok::Constants::protonMass();
				}
				if(_massBeamAddr== nullptr){
					std::cout << "INFO: No beam mass is given to reconstruct the beam energy -> using pion mass" << std::endl;
					_massBeamAddr= &antok::Constants::chargedPionMass();
				}
			}

			virtual ~GetBeamLorentzVec() { }

			bool operator() () {
				TVector3 p3Beam(_gradx, _grady, 1.);
				_out = antok::getBeamEnergy(p3Beam, _xLorentzVec, (*_massBeamAddr), (*_massTargetAddr));
				return true;
			}

		  private:

			const double& _gradx;
			const double& _grady;
			const TLorentzVector& _xLorentzVec;
			const double* _massBeamAddr;
			TLorentzVector& _out;
			const double* _massTargetAddr;

		};

		class GetTs : public Function
		{

		  public:

			GetTs(TLorentzVector* xLorentzVec, TLorentzVector* beamLorentzVec, double* tAddr, double* tMinAddr, double* tPrimeAddr, const double* targetMassAddr)
				: _xLorentzVec(*xLorentzVec),
				  _beamLorentzVec(*beamLorentzVec),
				  _t(*tAddr),
				  _tMin(*tMinAddr),
				  _tPrime(*tPrimeAddr),
				  _targetMassAddr(targetMassAddr){
				if(_targetMassAddr == nullptr){
					std::cout << "INFO: No target mass given to calculate t' -> Using approximation." << std::endl;
				}
			}

			virtual ~GetTs() { }

			bool operator() () {
			    _t = (_beamLorentzVec - _xLorentzVec).Mag2();


			    if(_targetMassAddr != NULL){
			    	const double EB_Lab = _beamLorentzVec.E();
			    	// calculate center of mass system quantities to calculate tmin
					const double mT2 = (*_targetMassAddr) * (*_targetMassAddr); // target mass
					const double mX2 = _xLorentzVec.M2(); // X-state mass
					const double mB2 = _beamLorentzVec.M2(); // beam mass
					const double s = mB2 + mT2 + 2.0*EB_Lab*sqrt(mT2);
					const double EB_CM = (s - mT2 + mB2)/(2.0*sqrt(s));
					const double EX_CM = (s + mX2 - mT2)/(2.0*sqrt(s));
					const double pB_CM = sqrt(EB_CM*EB_CM - mB2);
					const double pX_CM = sqrt(EX_CM*EX_CM - mX2);

					_tMin = -mB2 - mX2 + 2.0*( EB_CM*EX_CM - pB_CM*pX_CM );
			    } else {
			    	// use approximation
			    	_tMin = std::fabs((std::pow(_xLorentzVec.M2() - _beamLorentzVec.M2(), 2)) / (4. * _beamLorentzVec.Vect().Mag2()));
			    }

			    _tPrime = std::fabs(_t) - _tMin;
				return true;
			}

		  private:

			const TLorentzVector& _xLorentzVec;
			const TLorentzVector& _beamLorentzVec;
			double& _t;
			double& _tMin;
			double& _tPrime;
			const double* _targetMassAddr;

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



		template< typename T>
		class Abs : public Function
		{

		  public:

			Abs(T* inAddr, double* outAddr)
				: _inAddr(inAddr),
				  _outAddr(outAddr) { }

			virtual ~Abs() { }

			bool operator() () {
				(*_outAddr) = __abs(*_inAddr);
				return true;
			}

		  private:

			T const*const _inAddr;
			double* _outAddr;

			double __abs (const TVector3 & x){ return x.Mag(); }
			template <typename Q, class Dummy=int>
			double __abs(Q const& x){ return std::fabs( x ); }

		};

		template< typename T>
		class Log : public Function
		{

		  public:

			Log(T* inAddr, double* outAddr)
				: _inAddr(inAddr),
				  _outAddr(outAddr) { }

			virtual ~Log() { }

			bool operator() () {
				(*_outAddr) = std::log(*_inAddr);
				return true;
			}

		  private:

			T const*const _inAddr;
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
		class GetTVector3FromTLorenzVector: public Function
		{

		  public:

			GetTVector3FromTLorenzVector(const TLorentzVector* lvAddr, TVector3* vAddr)
				: _lv(*lvAddr),
				  _v(*vAddr){}

			virtual ~GetTVector3FromTLorenzVector() { }

			bool operator() () {
				_v = _lv.Vect();
				return true;
			}

		  private:

			const TLorentzVector& _lv;
			TVector3& _v;

		};

	}

}

#endif

